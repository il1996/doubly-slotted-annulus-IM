function mesh = mesh_refined(M, G, nc, i3, theta_rot, nr, ibar)
%MESH_REFINED  Grille polaire 2D conforme raffinable (mailleur, A2/A3).
%
%   mesh = mec.mesh_refined(M,G,nc,i3,theta_rot,nr,ibar) construit un réseau
%   de réluctances polaire à DEUX résolutions réglables :
%     nc : colonnes tangentielles par dent ET par encoche (conforme) ;
%     nr : couches RADIALES par région (dents et culasses).
%   Chaque cellule est un élément à 4 branches (radial haut/bas, tangentiel
%   gauche/droite), fer ou air selon la région (dent/encoche/culasse) —
%   l'élément de Perho/Silva. L'entrefer est traité par couplage physique
%   entre les cellules de bec stator et rotor (permet la rotation).
%
%   Entrées :
%     nc,nr      : résolutions tangentielle / radiale (>=1)
%     i3         : courants stator instantanés [iA;iB;iC] (FMM de bobinage)
%     theta_rot  : position rotorique [rad] (0 défaut)
%     nr         : couches radiales par région (2 défaut)
%     ibar       : courants de barre instantanés (Nr x 1) [A] (charge ; [] à vide)
%
%   Sortie : mesh (listes de branches pour mec.solve_mesh) + métadonnées
%   (nœuds de becs, positions d'entrefer, dth) pour Bg et couple.
%
%   Voir aussi : mec.solve_mesh, RUN_MESH, RUN_RIPPLE.

mu0 = 4*pi*1e-7;
if nargin<5 || isempty(theta_rot), theta_rot=0; end
if nargin<6 || isempty(nr), nr=2; end
if nargin<7, ibar=[]; end
Ns=M.Ns; Nr=M.Nr; L=M.L; kFe=M.kFe; g=M.g;
Rs=G.Rs; Rr=G.Rr; Rg=(Rs+Rr)/2;

ft_s=G.bts/G.taus;  ftgap_s=(G.taus-M.bs0)/G.taus;
ft_r=G.btr/G.taur;  ftgap_r=(G.taur-M.br0)/G.taur;
nst=nr; nys=max(1,ceil(nr/2));      % couches dents / culasse stator
nrt=nr; nyr=max(1,ceil(nr/2));      % couches dents / culasse rotor
Ls=nst+nys; Lr=nrt+nyr;

% --- Grilles angulaires conformes ---
[ths,dths,toothS,pidxS]=conform_ring(Ns,ft_s,nc,0);
[thr,dthr,toothR,pidxR]=conform_ring(Nr,ft_r,nc,theta_rot);
Ms=numel(ths); Mr=numel(thr);

% --- Rayons de couches (lay=1 = coté entrefer) ---
% stator : dents [Rs, Rs+hs] puis culasse [Rs+hs, Dso/2]
rs_edges=[Rs + (0:nst)/nst*G.hs, Rs+G.hs + (1:nys)/nys*G.hys];
rs_c=0.5*(rs_edges(1:end-1)+rs_edges(2:end));          % rayons centres (Ls)
rs_th=diff(rs_edges);                                   % épaisseurs (Ls)
% rotor : dents [Rr-hr, Rr] puis culasse vers l'arbre
rr_edges=[Rr - (0:nrt)/nrt*G.hr, Rr-G.hr - (1:nyr)/nyr*G.hyr];
rr_c=0.5*(rr_edges(1:end-1)+rr_edges(2:end));
rr_th=abs(diff(rr_edges));

% --- MAILLAGE CONFORME PAR COUCHE ------------------------------------
%  La repartition dent/encoche de CHAQUE couche radiale suit la largeur
%  REELLE du profil d'encoche a ce rayon (stator tc6 : isthme/biseau/
%  trapeze/fond ; rotor : goutte a epaulement) au lieu d'une fraction
%  constante : les cellules epousent les contours (plus de vides entre
%  maillage et geometrie) ET les sections de dents deviennent exactes en
%  profondeur (dent stator PARALLELE reproduite par construction).
bslot_s=@(y) (y<M.hs0).*M.bs0 ...
  +(y>=M.hs0 & y<M.hs0+M.hs1).*(M.bs0+(M.bs1-M.bs0).*(y-M.hs0)./max(M.hs1,1e-9)) ...
  +(y>=M.hs0+M.hs1 & y<M.hs0+M.hs1+M.hs2).*(M.bs1+(M.bs2-M.bs1).*(y-M.hs0-M.hs1)./M.hs2) ...
  +(y>=M.hs0+M.hs1+M.hs2).*M.bs2;
r1b=M.br1/2; r2b=M.br2/2; d0b=M.hr0+M.hr01; c1b=d0b+r1b; c2b=c1b+M.hr1;
bslot_r=@(y) (y<d0b).*M.br0 ...
  +(y>=d0b & y<c1b).*(2*sqrt(max(r1b^2-(c1b-y).^2,0))) ...
  +(y>=c1b & y<c2b).*(2*r1b+(2*r2b-2*r1b).*(y-c1b)./M.hr1) ...
  +(y>=c2b).*(2*sqrt(max(r2b^2-(min(y,c2b+r2b)-c2b).^2,0)));
THS=repmat(ths(:).',Ls,1); DTHS=repmat(dths(:).',Ls,1);
for la=1:nst
    tau_l=2*pi*rs_c(la)/Ns;
    ftl=min(max(1-bslot_s(rs_c(la)-Rs)/tau_l,0.15),0.92);
    [t_,d_]=conform_ring(Ns,ftl,nc,0);
    THS(la,:)=t_(:).'; DTHS(la,:)=d_(:).';
end
for la=nst+1:Ls, THS(la,:)=THS(nst,:); DTHS(la,:)=DTHS(nst,:); end
THR=repmat(thr(:).',Lr,1); DTHR=repmat(dthr(:).',Lr,1);
for la=1:nrt
    tau_l=2*pi*rr_c(la)/Nr;
    ftl=min(max(1-bslot_r(Rr-rr_c(la))/tau_l,0.15),0.92);
    [t_,d_]=conform_ring(Nr,ftl,nc,theta_rot);
    THR(la,:)=t_(:).'; DTHR(la,:)=d_(:).';
end
for la=nrt+1:Lr, THR(la,:)=THR(nrt,:); DTHR(la,:)=DTHR(nrt,:); end
ths=THS(1,:).'; dths=DTHS(1,:).';        % couche d'entrefer (gap, exports)
thr=THR(1,:).'; dthr=DTHR(1,:).';
ft1_s=sum(DTHS(1,toothS))/(2*pi);        % fractions de dent en couche 1
ft1_r=sum(DTHR(1,toothR))/(2*pi);

% --- Indices de nœuds : S(col,lay), R(col,lay) ---
Sid=@(c,la)(la-1)*Ms+c;
Rid=@(c,la)Ms*Ls+(la-1)*Mr+c;
mesh.Nnodes=Ms*Ls+Mr*Lr; mesh.ref=Rid(1,Lr);

% --- FMM stator par colonne (conducteurs dans colonnes d'encoche) ---
W=mec.winding(M); icond=i3(:)/M.a;
Cc=zeros(Ms,3); slotcolS=~toothS;
for s=1:Ns
    cols=find(pidxS==s & slotcolS); if isempty(cols), continue; end
    for ph=1:3, Cc(cols,ph)=W.C(s,ph)/numel(cols); end
end
FcolS=cumsum(Cc*icond); FcolS=FcolS-mean(FcolS);
% --- FMM rotor par colonne (barres dans colonnes d'encoche rotor) ---
FcolR=zeros(Mr,1);
if ~isempty(ibar)
    Ar=zeros(Mr,1); slotcolR=~toothR;
    for m=1:Nr
        cols=find(pidxR==m & slotcolR); if isempty(cols), continue; end
        Ar(cols)=ibar(m)/numel(cols);
    end
    FcolR=cumsum(Ar); FcolR=FcolR-mean(FcolR);
end

% ============ Assemblage des branches ============
% préallocation généreuse
maxB=2*(Ms*Ls+Mr*Lr)+8*Ms;
a=zeros(maxB,1); b=a; iron=a; ll=a; AA=a; PP=a; EE=a; nb=0;
    function push(na,nb_,ir,len,area,perm,emf)
        nb=nb+1; a(nb)=na; b(nb)=nb_; iron(nb)=ir;
        ll(nb)=len; AA(nb)=area; PP(nb)=perm; EE(nb)=emf;
    end
wSm=@(k)mod(k-1,Ms)+1; wRm=@(m)mod(m-1,Mr)+1;

% ---- STATOR ----
for c=1:Ms
    isT=toothS(c);
    % radial (lay -> lay+1)
    for la=1:Ls-1
        rmid=0.5*(rs_c(la)+rs_c(la+1)); len=abs(rs_c(la)-rs_c(la+1));
        toothReg=(la<=nst-1) || (la==nst); % branche dans/vers région dents
        inYoke=(la>=nst+1);
        dthm=0.5*(DTHS(la,c)+DTHS(la+1,c));      % largeur conforme moyenne
        area=dthm*rmid*kFe*L;                    % section pour flux radial
        if inYoke || (toothReg && isT)
            emf=0;
            if toothReg && isT && la==1, emf=FcolS(c); end   % FMM injectée
            push(Sid(c,la),Sid(c,la+1),1,len,area,0,emf);
        else
            P=mu0*(dthm*rmid*L)/len;              % chemin d'encoche (air)
            push(Sid(c,la),Sid(c,la+1),0,0,0,P,0);
        end
    end
    % tangentiel (c -> c+1) par couche
    c2=wSm(c+1);
    for la=1:Ls
        rlay=rs_c(la); thick=rs_th(la);
        arc=0.5*(DTHS(la,c)+DTHS(la,c2))*rlay;
        area=thick*kFe*L;
        inYoke=(la>=nst+1);
        sameT=toothS(c)&&toothS(c2)&&(pidxS(c)==pidxS(c2));
        if inYoke || sameT
            push(Sid(c,la),Sid(c2,la),1,arc,area,0,0);
        else
            Pt=mu0*(thick*L)/arc;                 % fuite tangentielle (air)
            push(Sid(c,la),Sid(c2,la),0,0,0,Pt,0);
        end
    end
end

% ---- ROTOR ----
for c=1:Mr
    isT=toothR(c);
    for la=1:Lr-1
        rmid=0.5*(rr_c(la)+rr_c(la+1)); len=abs(rr_c(la)-rr_c(la+1));
        toothReg=(la<=nrt);
        inYoke=(la>=nrt+1);
        dthm=0.5*(DTHR(la,c)+DTHR(la+1,c));      % largeur conforme moyenne
        area=dthm*rmid*kFe*L;
        if inYoke || (toothReg && isT)
            emf=0;
            if toothReg && isT && la==1, emf=FcolR(c); end
            push(Rid(c,la),Rid(c,la+1),1,len,area,0,emf);
        else
            P=mu0*(dthm*rmid*L)/len;
            push(Rid(c,la),Rid(c,la+1),0,0,0,P,0);
        end
    end
    c2=wRm(c+1);
    for la=1:Lr
        rlay=rr_c(la); thick=rr_th(la);
        arc=0.5*(DTHR(la,c)+DTHR(la,c2))*rlay;
        area=thick*kFe*L;
        inYoke=(la>=nrt+1);
        sameT=toothR(c)&&toothR(c2)&&(pidxR(c)==pidxR(c2));
        if inYoke || sameT
            push(Rid(c,la),Rid(c2,la),1,arc,area,0,0);
        else
            Pt=mu0*(thick*L)/arc;
            push(Rid(c,la),Rid(c2,la),0,0,0,Pt,0);
        end
    end
end

% ---- ENTREFER : bec stator (lay1, dents) <-> bec rotor (lay1, dents) ----
%  Perméance C1 ANALYTIQUE : noyau gaussien de largeur PHYSIQUE sigma0
%  (indépendante du maillage — correction M1/A1). On stocke dLambda/dtheta
%  par branche (dérivée analytique, dsigma/dtheta = Rg) pour le couple par
%  travaux virtuels analytiques (A5), sans re-maillage.
Cper=2*pi*Rg; xs=ths*Rg; xr=thr*Rg; wxs=dths*Rg;
% Largeur de couplage C1 = largeur PHYSIQUE de la perméance dent-à-dent,
% IDENTIFIÉE PAR FEM (banc à deux dents, mec.fem_airgap_ident / RUN_FEM_IDENT,
% à la Gyselinck [29]) : sigma0 = 0.54*tau. C'est LE paramètre qui gouverne
% l'ondulation de couple (A1). NB : à ce sigma0 physique, l'ondulation MEC
% reste sur-estimée (~41 %) — résidu structurel du réseau dent-à-dent
% monocouche (défaut M6), non supprimable par la largeur d'entrefer seule.
sigma0 = 0.54 * 0.5*(G.taus + G.taur);     % FEM-identifié (0.54*tau)
if isfield(M,'opt') && isfield(M.opt,'gap_sigma0') && ~isempty(M.opt.gap_sigma0)
    sigma0 = M.opt.gap_sigma0;             % override (étude de sensibilité)
end
win = 3*sigma0;
rtooth=find(toothR); xr_t=xr(rtooth);
gdLall=zeros(maxB,1);                        % dLambda/dtheta par branche
Ng=0;
if isfield(M,'opt')&&isfield(M.opt,'gap_layers')&&~isempty(M.opt.gap_layers)
    Ng=M.opt.gap_layers;
end
gap0=nb;

%  DEFAUT = OPERATEUR HARMONIQUE. L'entrefer par branches discretes (branche
%  ci-dessous) est un modele de COMPARAISON : il ferme la couronne par un
%  noyau gaussien de largeur sigma0 ajustee, la ou l'operateur la resout.
%  Il faut donc opt.gap_fourier = false pour le retrouver, et non l'inverse.
gapF = true;
if isfield(M,'opt')&&isfield(M.opt,'gap_fourier')&&~isempty(M.opt.gap_fourier)
    gapF = logical(M.opt.gap_fourier);
end
if gapF
% ---------- P1 : ENTREFER HARMONIQUE (couronne de Laplace, leve M6) ----------
%  Des noeuds de SURFACE sont crees a r=Rs (stator) et r=Rr (rotor), relies
%  aux noeuds de couche 1 par la demi-branche radiale de leur cellule (fer
%  non lineaire sous les dents, AIR sous les ouvertures d'encoche — cette
%  demi-hauteur d'air est essentielle : sans elle l'interieur d'encoche
%  serait colle a la surface d'entrefer). Le couplage entre les deux
%  surfaces est la matrice de permeance harmonique EXACTE de la couronne
%  d'air (mec.airgap_fourier) : Carter, frange, composante tangentielle Bt
%  et attenuation des harmoniques de denture EMERGENT de la solution — ni
%  kC ni sigma0. Le couple s'obtient par MST harmonique analytique
%  (mesh.gapF.torque) ; il n'y a plus de branches d'entrefer discretes.
Ng0 = Ms*Ls + Mr*Lr;
mesh.Nnodes = Ng0 + Ms + Mr;
for c=1:Ms
    len=rs_th(1)/2; rmid=Rs+rs_th(1)/4;
    if toothS(c)
        push(Sid(c,1),Ng0+c,1,len,dths(c)*rmid*kFe*L,0,0);
    else
        push(Sid(c,1),Ng0+c,0,0,0,mu0*(dths(c)*rmid*L)/len,0);
    end
end
for c=1:Mr
    len=rr_th(1)/2; rmid=Rr-rr_th(1)/4;
    if toothR(c)
        push(Rid(c,1),Ng0+Ms+c,1,len,dthr(c)*rmid*kFe*L,0,0);
    else
        push(Rid(c,1),Ng0+Ms+c,0,0,0,mu0*(dthr(c)*rmid*L)/len,0);
    end
end
gap0=nb;                                  % pas de branches d'entrefer discretes
NhF=100;
if isfield(M.opt,'gap_fourier_Nh')&&~isempty(M.opt.gap_fourier_Nh)
    NhF=M.opt.gap_fourier_Nh;
end
% Grille de SURFACE : les arcs de colonnes sont remappes aux largeurs
% PHYSIQUES de la surface d'entrefer — face de dent (avec becs) = taus-bs0,
% ouverture = bs0 — et non aux largeurs du CORPS d'encoche (bts/taus) qui
% servent au maillage interne. Sans ce remappage, la surface serait « trop
% encochee » (ouverture ~45 % du pas au lieu de ~18 %) : Bg1 chuterait de
% ~25 % et la denture serait sur-representee. Etirement symetrique autour du
% centre de chaque bloc (dents/ouverture) -> blocs exactement contigus.
[thsF,dthsF]=surface_grid(ths,dths,toothS,ftgap_s/ft1_s,(1-ftgap_s)/(1-ft1_s));
[thrF,dthrF]=surface_grid(thr,dthr,toothR,ftgap_r/ft1_r,(1-ftgap_r)/(1-ft1_r));
AF=mec.airgap_fourier(thsF,dthsF,thrF,dthrF,Rs,Rr,L,NhF);
AF.th_s_face=thsF; AF.dth_s_face=dthsF; AF.th_r_face=thrF; AF.dth_r_face=dthrF;
ids=[Ng0+(1:Ms), Ng0+Ms+(1:Mr)];
[iiY,jjY]=ndgrid(ids,ids);
mesh.Yx=sparse(iiY(:),jjY(:),AF.Y(:),mesh.Nnodes,mesh.Nnodes);
AF.ids=ids; mesh.gapF=AF;

elseif Ng==0
% ---------- ENTREFER MONOCOUCHE (défaut) : dent stator <-> dent rotor ----------
for c=1:Ms
    if ~toothS(c), continue; end
    d=mod(xr_t-xs(c)+Cper/2,Cper)-Cper/2;    % sigma signé
    near=find(abs(d)<=win); if isempty(near), continue; end
    sig=d(near);
    k=exp(-(sig/sigma0).^2);                  % noyau gaussien C1
    S=sum(k); if S<=0, continue; end
    dk=k.*(-2*sig/sigma0^2)*Rg;               % dk/dtheta
    dS=sum(dk);
    P0=mu0*L*wxs(c)/g*(ftgap_s/ft1_s);        % perméance totale de colonne
    Lam=P0*k/S;                               % perméances (somme = P0)
    dLam=P0*(dk*S - k*dS)/S^2;                % dLambda/dtheta (analytique)
    for t=1:numel(near)
        push(Sid(c,1),Rid(rtooth(near(t)),1),0,0,0,Lam(t),0);
        gdLall(nb)=dLam(t);
    end
end
else
% ---------- ENTREFER MULTICOUCHE : anneau mi-entrefer + étalement tangentiel ----
%  Le flux traverse un demi-entrefer stator->anneau (fixe) puis un demi-entrefer
%  anneau->rotor (glissant). L'anneau porte des branches tangentielles d'air
%  (hauteur physique h_eff) qui LISSENT la distribution de flux avant le
%  couplage glissant -> teste si l'étalement de gap réduit l'ondulation (M6).
Ng0 = Ms*Ls + Mr*Lr;  MGid=@(k)Ng0+k;         % noeuds d'anneau mi-entrefer
mesh.Nnodes = Ng0 + Ms;
heff = g;                                     % hauteur d'étalement (physique ~ g)
if isfield(M,'opt')&&isfield(M.opt,'gap_heff')&&~isempty(M.opt.gap_heff)
    heff=M.opt.gap_heff;
end
gh = g/2;                                     % demi-entrefer
% branches tangentielles de l'anneau (air) : étalement
for c=1:Ms
    c2=wSm(c+1);
    arc=0.5*(dths(c)+dths(c2))*Rg;
    push(MGid(c),MGid(c2),0,0,0,mu0*heff*L/arc,0);
end
% stator dent -> anneau (demi-entrefer, FIXE, dLambda/dtheta=0)
for c=1:Ms
    if ~toothS(c), continue; end
    Phalf=mu0*L*wxs(c)/gh*(ftgap_s/ft1_s);
    push(Sid(c,1),MGid(c),0,0,0,Phalf,0);
end
% anneau -> dent rotor (demi-entrefer, GLISSANT, dLambda/dtheta analytique)
for c=1:Ms
    if ~toothS(c), continue; end
    d=mod(xr_t-xs(c)+Cper/2,Cper)-Cper/2;
    near=find(abs(d)<=win); if isempty(near), continue; end
    sig=d(near); k=exp(-(sig/sigma0).^2); S=sum(k); if S<=0, continue; end
    dk=k.*(-2*sig/sigma0^2)*Rg; dS=sum(dk);
    P0=mu0*L*wxs(c)/gh*(ftgap_s/ft1_s);       % demi-entrefer
    Lam=P0*k/S; dLam=P0*(dk*S-k*dS)/S^2;
    for t=1:numel(near)
        push(MGid(c),Rid(rtooth(near(t)),1),0,0,0,Lam(t),0);
        gdLall(nb)=dLam(t);
    end
end
end

% troncature
a=a(1:nb); b=b(1:nb); iron=logical(iron(1:nb));
ll=ll(1:nb); AA=AA(1:nb); PP=PP(1:nb); EE=EE(1:nb);
mesh.a=a; mesh.b=b; mesh.iron=iron; mesh.l=ll; mesh.A=AA; mesh.P=PP; mesh.E=EE;

% métadonnées
mesh.Ms=Ms; mesh.Mr=Mr; mesh.nc=nc; mesh.nr=nr; mesh.Ls=Ls; mesh.Lr=Lr;
mesh.th_s=ths; mesh.dth_s=dths; mesh.toothS=toothS;
mesh.Sid=Sid; mesh.Rid=Rid; mesh.Rg=Rg; mesh.p=M.p; mesh.L=L; mesh.Ns=Ns;
% métadonnées géométriques pour post-traitement graphique (cartes de B)
mesh.th_r=thr; mesh.dth_r=dthr; mesh.toothR=toothR;
mesh.rs_edges=rs_edges; mesh.rr_edges=rr_edges;
% grilles CONFORMES par couche (largeurs reelles du profil d'encoche)
mesh.THS=THS; mesh.DTHS=DTHS; mesh.THR=THR; mesh.DTHR=DTHR;
mesh.gapfirst=gap0+1;                 % 1er indice de branche d'entrefer
mesh.gap_dLdth=gdLall(gap0+1:nb);     % dLambda/dtheta des branches d'entrefer
end

% ======================================================================
function [th,dth,tooth,pidx]=conform_ring(N,ft,nc,shift)
dpitch=2*pi/N; wt=ft*dpitch/nc; ws=(1-ft)*dpitch/nc;
K=2*nc*N; th=zeros(K,1); dth=th; tooth=false(K,1); pidx=zeros(K,1);
edge=0; idx=0;
for pp=1:N
    for c=1:nc, idx=idx+1; th(idx)=edge+wt/2; dth(idx)=wt; tooth(idx)=true;  pidx(idx)=pp; edge=edge+wt; end
    for c=1:nc, idx=idx+1; th(idx)=edge+ws/2; dth(idx)=ws; tooth(idx)=false; pidx(idx)=pp; edge=edge+ws; end
end
th=th+shift;
end

function [thF,dthF]=surface_grid(th,dth,tooth,aT,aO)
%SURFACE_GRID  Remappe les arcs de colonnes aux largeurs de SURFACE (becs).
%   Chaque bloc contigu (dents ou ouverture) est dilate du facteur aT/aO
%   autour de SON centre : la conservation par pas (ft*aT+(1-ft)*aO=1)
%   garantit des blocs contigus sans recouvrement.
th=th(:); dth=dth(:); tooth=logical(tooth(:));
thF=th; dthF=dth;
K=numel(th); i=1;
while i<=K
    j=i; while j<K && tooth(j+1)==tooth(i), j=j+1; end
    al=aO; if tooth(i), al=aT; end
    e0=th(i)-dth(i)/2; e1=th(j)+dth(j)/2; m=0.5*(e0+e1);
    thF(i:j)=m+(th(i:j)-m)*al;
    dthF(i:j)=dth(i:j)*al;
    i=j+1;
end
end
