function draw_cross_section(me, Se, M, G, W, i3, ib, opt) %#ok<INUSL>
%DRAW_CROSS_SECTION  Cartes de la machine STRUCTUREES facon ANSYS.
%
%   mec.draw_cross_section(me,Se,M,G,W,i3,ib,opt) produit, a partir d'une
%   resolution de champ (mec.mesh_refined + mec.solve_mesh, entrefer Fourier),
%   DEUX figures au rendu comparable aux cartes ANSYS Maxwell :
%
%     1. CARTE |B| STRUCTUREE : cellules de fer colorees (jet, 0..Bmax),
%        encoches remplies (B~0) avec leurs VRAIS contours (stator tc6 :
%        isthme + biseau + trapeze + fond arrondi ; rotor : goutte a
%        EPAULEMENT, cercle superieur tangent a l'ouverture), arbre et
%        frontieres — et DISTRIBUTION DU CHAMP par LIGNES DE FLUX :
%        contours de la fonction de flux psi reconstruite par integration
%        des flux de branches (equivalent MEC des lignes iso-A d'ANSYS).
%        La fermeture de l'integration (garantie par la conservation du
%        flux, KCL) est verifiee et affichee.
%     2. CARTE DE DENSITE DE COURANT J : barres rotoriques au courant
%        instantane J = i_bar/Abar, encoches stator par COUCHE (bobinage
%        double couche reel reconstruit de M.wind), echelle symetrique.
%
%   opt.Bmax (2.2), opt.Jmax (8e6), opt.nlines (24), opt.prefix ('').
%
%   LANGUE. Les libelles affiches sont surchargeables, de sorte que la
%   fonction de trace ne fige aucune langue :
%       opt.titleB / opt.titleJ : titres des deux figures
%       opt.labB   / opt.labJ   : libelles d'echelle de couleur
%       opt.nameB  / opt.nameJ  : noms de fenetre
%   RUN_ARTICLE passe l'ANGLAIS ; les defauts ci-dessous sont un repli
%   francais. Ils sont en ASCII pur a dessein : un accent dans un litteral
%   .m depend de l'encodage du fichier, et rien ne garantit qu'il arrive
%   intact dans un PNG destine a la publication.
%   Ces libelles n'ont aucun lien avec les noms de fichier choisis par
%   l'appelant.
%
%   Voir aussi : RUN_ARTICLE (section IV), mec.mesh_refined.

if nargin<8, opt=struct(); end
Bmax=getdef(opt,'Bmax',2.2);
Jmax=getdef(opt,'Jmax',8e6);
nlin=getdef(opt,'nlines',24);
pfx =getdef(opt,'prefix','');
ttlB=getdef(opt,'titleB','|B| et lignes de flux (geometrie reelle des encoches)');
ttlJ=getdef(opt,'titleJ','Densite de courant instantanee (barres + 2 couches stator)');
labB=getdef(opt,'labB','|B| [T]');
labJ=getdef(opt,'labJ','J [A/m^2]');
namB=getdef(opt,'nameB','Carte |B| structuree');
namJ=getdef(opt,'nameJ','Carte J');

Ms=me.Ms; Mr=me.Mr; Ls=me.Ls; Lr=me.Lr;
rs=me.rs_edges(:).'; rr=me.rr_edges(:).';
ths=me.th_s(:).'; dths=me.dth_s(:).';
thr=me.th_r(:).'; dthr=me.dth_r(:).';
Ng0=Ms*Ls+Mr*Lr;
rot=thr(1)-dthr(1)/2;                        % position rotorique (theta_rot)
% Centres angulaires des encoches COHERENTS avec le maillage : dans chaque
% pas, conform_ring place le bloc DENT d'abord puis le bloc ENCOCHE
% -> centre d'encoche = debut du pas + (1+ft)/2 * pas  (ft = bt/tau).
% (Sans cela, les contours dessines sont decales de ~(ft-? )*pas vs les
% cellules et les lignes de flux : c'est le "decalage" visible.)
fts=G.bts/G.taus; ftr=G.btr/G.taur;
thS=@(s2) 2*pi*((s2-1)+(1+fts)/2)/M.Ns;
thR=@(m2) rot + 2*pi*((m2-1)+(1+ftr)/2)/M.Nr;

% ============ 1. flux de branches ranges en grilles ============
nb=me.gapfirst-1; a=me.a(1:nb); b=me.b(1:nb); Phi=Se.Phi(1:nb);
PhirS=zeros(Ms,Ls-1); PhitS=zeros(Ms,Ls);    % radial la->la+1 ; tangentiel c->c+1
PhirR=zeros(Mr,Lr-1); PhitR=zeros(Mr,Lr);
PsurS=zeros(Ms,1);    PsurR=zeros(Mr,1);     % cellule(couche 1) -> surface
for k=1:nb
    na=a(k); nb2=b(k);
    if nb2>Ng0
        if nb2<=Ng0+Ms, PsurS(nb2-Ng0)=Phi(k);
        else,           PsurR(nb2-Ng0-Ms)=Phi(k); end
    elseif na<=Ms*Ls
        c=mod(na-1,Ms)+1; la=ceil(na/Ms);
        c2=mod(nb2-1,Ms)+1;
        if c==c2, PhirS(c,la)=Phi(k); else, PhitS(c,la)=Phi(k); end
    else
        na_=na-Ms*Ls; nb_=nb2-Ms*Ls;
        c=mod(na_-1,Mr)+1; la=ceil(na_/Mr);
        c2=mod(nb_-1,Mr)+1;
        if c==c2, PhirR(c,la)=Phi(k); else, PhitR(c,la)=Phi(k); end
    end
end

% ============ 2. fonction de flux psi (lignes de flux) ============
%  psi aux COINS de cellules ; flux (+r) d'une face horizontale =
%  psi(c+1)-psi(c) ; flux (+theta) d'une face verticale = psi(la)-psi(la+1).
%  L'integrabilite est garantie par KCL ; la fermeture (wrap) est verifiee.
[PSIs,errS]=stream_block(-PsurS,PhirS,PhitS,+1);   % stator : +r vers l'exterieur
[PSIr,errR]=stream_block(+PsurR,PhirR,PhitR,-1);   % rotor  : couches vers -r
PSIr=PSIr-mean(PSIr(1:end-1,1))+mean(PSIs(1:end-1,1));  % continuite via l'entrefer
fprintf('   [draw_cross_section] fermeture de psi : stator %.2e, rotor %.2e Wb\n',errS,errR);

% grilles CONFORMES par couche (repli : grille unique repetee)
if isfield(me,'THS')
    THS=me.THS; DTHS=me.DTHS; THR=me.THR; DTHR=me.DTHR;
else
    THS=repmat(ths,Ls,1); DTHS=repmat(dths,Ls,1);
    THR=repmat(thr,Lr,1); DTHR=repmat(dthr,Lr,1);
end

% ============ FIGURE 1 : |B| structuree + lignes de flux ============
figure('Name',[pfx namB],'Position',[60 40 880 800]); hold on;
Bn=cellB(me,Se,THS,THR);
paint_cells(Bn,THS,DTHS,THR,DTHR,rs,rr,Ms,Mr,Ls,Lr);
% encoches et barres remplies (B ~ 0) avec contours reels
for s2=1:M.Ns
    P=stator_poly(M,thS(s2));
    fill(P(1,:),P(2,:),0,'EdgeColor',[.75 .75 .75],'LineWidth',0.4);
end
for m2=1:M.Nr
    P=rotor_poly(M,thR(m2));
    fill(P(1,:),P(2,:),0,'EdgeColor',[.75 .75 .75],'LineWidth',0.4);
end
fill_disk(M.Dsh/2,[0.80 0.87 0.97]);
colormap(jet(11)); clim([0 Bmax]);
cb=colorbar; cb.Label.String=labB;
lv=linspace(min([PSIs(:);PSIr(:)]),max([PSIs(:);PSIr(:)]),nlin);
contour_block(PSIs,THS,DTHS,rs,lv);
contour_block(PSIr,THR,DTHR,rr,lv);
draw_circles(M,G);
axis equal off;
title([pfx ttlB]);

% ============ FIGURE 2 : densite de courant J ============
figure('Name',[pfx namJ],'Position',[110 60 880 800]); hold on;
fill_annulus(G.Rs,M.Dso/2,[0.92 0.92 0.94]);
fill_annulus(M.Dsh/2,G.Rr,[0.92 0.92 0.94]);
fill_disk(M.Dsh/2,[0.80 0.87 0.97]);
for m2=1:M.Nr
    P=rotor_poly(M,thR(m2));
    fill(P(1,:),P(2,:),ib(m2)/G.Abar,'EdgeColor',[.3 .3 .3],'LineWidth',0.3);
end
[Jd,Ju]=stator_layerJ(M,G,i3);
for s2=1:M.Ns
    [Pd,Pu]=stator_poly_layers(M,thS(s2));
    fill(Pd(1,:),Pd(2,:),Jd(s2),'EdgeColor',[.3 .3 .3],'LineWidth',0.3);
    fill(Pu(1,:),Pu(2,:),Ju(s2),'EdgeColor',[.3 .3 .3],'LineWidth',0.3);
end
colormap(jet(11)); clim([-Jmax Jmax]);
cb=colorbar; cb.Label.String=labJ;
draw_circles(M,G);
axis equal off;
title([pfx ttlJ]);
end

% ======================================================================
function v=getdef(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end

function [PSI,err]=stream_block(Fr1,Phir,Phit,srad)
%  PSI aux coins (K+1 x L+1). Fr1 = flux +r des faces de SURFACE (rangee 1) ;
%  Phir (K x L-1) = branches radiales la->la+1 (flux +r si srad=+1) ;
%  Phit (K x L)   = branches tangentielles c->c+1 (flux +theta).
    [K,Lm1]=size(Phir); L=Lm1+1;
    PSI=zeros(K+1,L+1);
    PSI(2:K+1,1)=cumsum(Fr1(:));                       % rangee de surface
    for la=1:L                                         % colonne 1 (arete theta_1)
        PSI(1,la+1)=PSI(1,la)-Phit(K,la);
    end
    for la=2:L                                         % rangees internes
        PSI(2:K+1,la)=PSI(1,la)+cumsum(srad*Phir(:,la-1));
    end
    PSI(2:K+1,L+1)=PSI(1,L+1);                         % frontiere externe : flux nul
    err=max(abs(PSI(K+1,1:L)-PSI(1,1:L)));             % fermeture (wrap)
end

function Bn=cellB(me,Se,THS,THR)
    nb=me.gapfirst-1; a=me.a(1:nb); b=me.b(1:nb);
    Bsig=zeros(nb,1); ok=me.A(1:nb)>0; Bsig(ok)=Se.Phi(ok)./me.A(ok);
    Ms=me.Ms; Mr=me.Mr; Ls=me.Ls; Lr=me.Lr;
    rs_c=0.5*(me.rs_edges(1:end-1)+me.rs_edges(2:end));
    rr_c=0.5*(me.rr_edges(1:end-1)+me.rr_edges(2:end));
    x=zeros(me.Nnodes,1); y=x;
    for la=1:Ls, id=(la-1)*Ms+(1:Ms); x(id)=rs_c(la)*cos(THS(la,:)); y(id)=rs_c(la)*sin(THS(la,:)); end
    for la=1:Lr, id=Ms*Ls+(la-1)*Mr+(1:Mr); x(id)=rr_c(la)*cos(THR(la,:)); y(id)=rr_c(la)*sin(THR(la,:)); end
    dx=x(b)-x(a); dy=y(b)-y(a); dl=hypot(dx,dy); dl(dl==0)=1;
    ux=dx./dl; uy=dy./dl;
    rn=hypot(x(a),y(a)); rn(rn==0)=1;
    isRad=abs(ux.*x(a)./rn+uy.*y(a)./rn)>0.7;
    N=me.Nnodes;
    acc=@(m)deal(accumarray([a(m);b(m)],[Bsig(m).*ux(m);Bsig(m).*ux(m)],[N 1]),...
                 accumarray([a(m);b(m)],[Bsig(m).*uy(m);Bsig(m).*uy(m)],[N 1]),...
                 accumarray([a(m);b(m)],1,[N 1]));
    [vxR,vyR,cR]=acc(isRad); [vxT,vyT,cT]=acc(~isRad);
    cR(cR==0)=1; cT(cT==0)=1;
    Bn=hypot(vxR./cR+vxT./cT,vyR./cR+vyT./cT);
end

function paint_cells(Bn,THS,DTHS,THR,DTHR,rs,rr,Ms,Mr,Ls,Lr)
%  peint TOUTES les cellules (fer ET air : l'air est ~0 T -> bleu sombre),
%  avec la grille CONFORME de chaque couche -> aucun vide residuel.
    [XX,YY,CC]=deal([]);
    for la=1:Ls
        [X,Y]=cellpoly(THS(la,:),DTHS(la,:),rs(la),rs(la+1));
        id=(la-1)*Ms+(1:Ms);
        XX=[XX X]; YY=[YY Y]; CC=[CC Bn(id).']; %#ok<AGROW>
    end
    for la=1:Lr
        [X,Y]=cellpoly(THR(la,:),DTHR(la,:),rr(la),rr(la+1));
        id=Ms*Ls+(la-1)*Mr+(1:Mr);
        XX=[XX X]; YY=[YY Y]; CC=[CC Bn(id).']; %#ok<AGROW>
    end
    patch(XX,YY,CC,'EdgeColor','none');
end

function [X,Y]=cellpoly(th,dth,r1,r2)
    t1=th-dth/2; t2=th+dth/2;
    X=[r1*cos(t1);r1*cos(t2);r2*cos(t2);r2*cos(t1)];
    Y=[r1*sin(t1);r1*sin(t2);r2*sin(t2);r2*sin(t1)];
end

function contour_block(PSI,TH,DTH,redge,lv)
%  contours de psi sur la grille CURVILIGNE des coins (aretes par couche)
    L=size(TH,1); K=size(TH,2);
    X=zeros(L+1,K+1); Y=X;
    for la=1:L+1
        g=min(la,L);
        e=[TH(g,:)-DTH(g,:)/2, TH(g,1)-DTH(g,1)/2+2*pi];
        X(la,:)=redge(la)*cos(e); Y(la,:)=redge(la)*sin(e);
    end
    contour(X,Y,PSI.',lv,'k','LineWidth',0.5);
end

function draw_circles(M,G)
    tt=linspace(0,2*pi,361);
    plot((M.Dso/2)*cos(tt),(M.Dso/2)*sin(tt),'k-','LineWidth',1.0);
    plot((M.Dsh/2)*cos(tt),(M.Dsh/2)*sin(tt),'k-','LineWidth',0.8);
    plot(G.Rs*cos(tt),G.Rs*sin(tt),'-','Color',[.5 .5 .5],'LineWidth',0.3);
    plot(G.Rr*cos(tt),G.Rr*sin(tt),'-','Color',[.5 .5 .5],'LineWidth',0.3);
end

function fill_annulus(r1,r2,col)
    tt=linspace(0,2*pi,181);
    fill([r1*cos(tt), r2*cos(fliplr(tt))],[r1*sin(tt), r2*sin(fliplr(tt))],...
        col,'EdgeColor','none');
end
function fill_disk(r,col)
    tt=linspace(0,2*pi,181);
    fill(r*cos(tt),r*sin(tt),col,'EdgeColor','none');
end

function P=stator_poly(M,thc)
%  profil tc6 (ouverture, biseau, trapeze, fond a coins arrondis rayon Rs)
    Rb=M.Ds/2; Rc=M.Rs;
    y1=M.hs0; y2=M.hs0+M.hs1; y3=y2+M.hs2; y4=y3+Rc;
    xh=[M.bs0/2, M.bs0/2, M.bs1/2, M.bs2/2];
    yh=[0,       y1,      y2,      y3     ];
    tt=linspace(0,pi/2,7);
    xc=(M.bs2/2-Rc)+Rc*cos(tt); yc=y3+Rc*sin(tt);
    xr=[xh, xc, 0];  yr=[yh, yc, y4];
    x=[xr, -fliplr(xr)]; y=[yr, fliplr(yr)];
    r=Rb+y; th=thc+x./r;
    P=[r.*cos(th); r.*sin(th)];
end

function [Pd,Pu]=stator_poly_layers(M,thc)
%  encoche scindee en 2 couches radiales egales du CORPS (double couche) :
%  Pu = couche cote entrefer, Pd = couche de fond
    Rb=M.Ds/2;
    y2=M.hs0+M.hs1; y3=y2+M.hs2; ym=0.5*(y2+y3);
    bm=M.bs1+(M.bs2-M.bs1)*(ym-y2)/(y3-y2);
    Pu=locpoly(Rb,thc,[ M.bs1/2,  bm/2, -bm/2, -M.bs1/2],[y2, ym, ym, y2]);
    Pd=locpoly(Rb,thc,[ bm/2, M.bs2/2, -M.bs2/2, -bm/2],[ym, y3, y3, ym]);
end
function P=locpoly(Rb,thc,xs,ys)
    r=Rb+ys; th=thc+xs./r;
    P=[r.*cos(th); r.*sin(th)];
end

function P=rotor_poly(M,thc)
%  goutte a EPAULEMENT (figure ANSYS) : cercle superieur tangent au bas de
%  l'ouverture (centre a Br1/2), flancs tangents, cercle inferieur Br2.
    Rr=M.Ds/2-M.g;
    r1=M.br1/2; r2=M.br2/2;
    d0=M.hr0+M.hr01; c1=d0+r1; c2=c1+M.hr1;
    sa=min(max((r1-r2)/M.hr1,-1),1); al=asin(sa);
    ttop=linspace(0,pi/2-al,9);
    xt=r1*sin(ttop); yt=c1-r1*cos(ttop);
    tbot=linspace(pi/2-al,pi,11);
    xb=r2*sin(tbot); yb=c2-r2*cos(tbot);
    xr=[M.br0/2, M.br0/2, xt, xb];
    yr=[0,       d0,      yt, yb];
    x=[xr, -fliplr(xr)]; y=[yr, fliplr(yr)];
    r=Rr-y; th=thc+x./r;
    P=[r.*cos(th); r.*sin(th)];
end

function [Jd,Ju]=stator_layerJ(M,G,i3)
%  densite de courant instantanee des 2 couches (bobinage reel M.wind)
    Ju=zeros(M.Ns,1); Jd=Ju;
    Ic=i3(:)/M.a;                              % courant par voie parallele
    Ahalf=G.Aslot_s/2;
    wN=@(x)mod(x-1,M.Ns)+1;
    for ph=1:3
        sg=1; if ph==3 && isfield(M.wind,'Cinv') && M.wind.Cinv, sg=-1; end
        pD=wN(M.wind.pDown+M.wind.shift(ph)); sD=sg*M.wind.sDown;
        pU=wN(M.wind.pUp  +M.wind.shift(ph)); sU=sg*M.wind.sUp;
        for k2=1:numel(pD), Jd(pD(k2))=Jd(pD(k2))+sD(k2)*M.Ntc*Ic(ph)/Ahalf; end
        for k2=1:numel(pU), Ju(pU(k2))=Ju(pU(k2))+sU(k2)*M.Ntc*Ic(ph)/Ahalf; end
    end
end
