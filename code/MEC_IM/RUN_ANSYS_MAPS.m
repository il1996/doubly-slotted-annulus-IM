%% RUN_ANSYS_MAPS  -  Cartes de champ et distribution du flux : MEC vs ANSYS
%
%  Comparaison QUANTITATIVE des champs entre le reseau MEC raffine et l'EF
%  (ANSYS Maxwell 2D), a partir des profils de mi-entrefer exportes par ANSYS
%  (Calculator Expressions Plot 1.tab : Distance, A [Wb/m], |B|, Br, Bt sur
%  1000 points, snapshot t = 2 s) pour les DEUX essais (a vide, en charge) :
%
%    1. profil d'induction radiale Br(theta) au mi-entrefer, superpose ;
%    2. spectre spatial de Br (fondamental + harmoniques de denture) ;
%    3. distribution du flux : A(theta) ANSYS vs flux cumule MEC ;
%    4. cartes 2D de densite de flux |B| du MEC (a comparer aux cartes
%       ANSYS des .docx) + carte vectorielle (direction du flux) ;
%    5. tableau des sondes locales ANSYS (m1..m8) vs B regionaux MEC.
%
%  METHODE D'APPARIEMENT. Le point de fonctionnement MEC est evalue au
%  glissement DEDUIT de la vitesse du snapshot ANSYS (1500,2 tr/min a vide ;
%  1469,8 tr/min en charge). L'origine angulaire d'ANSYS etant arbitraire
%  (position rotorique du snapshot), les profils MEC sont recales en phase
%  sur le FONDAMENTAL (rotation de Dphi/p) : la comparaison porte sur les
%  amplitudes et le contenu harmonique, pas sur l'origine des angles.
%  Le MEC ne calcule PAS de composante tangentielle dans l'entrefer
%  (branches purement radiales) : Bt_ANSYS est quantifie comme limite (M6).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
OUT=fileparts(mfilename('fullpath'));

fprintf('=== Cartes de champ MEC vs ANSYS (18,5 kW) ===\n\n');
M=mec.machine_18_5kW(); ctx=mec.build_context(M);
G=ctx.G; W=ctx.W; BH=ctx.BH; p=M.p; Nr=M.Nr; Ns=M.Ns; mu0=4*pi*1e-7;

%% 1. Profils ANSYS de mi-entrefer
P0=read_profile(fullfile(ROOT,'transitoire','a vide','Calculator Expressions Plot 1.tab'));
P1=read_profile(fullfile(ROOT,'transitoire','en charge','Calculator Expressions Plot 1.tab'));
Rmid=P0.C/(2*pi);
fprintf('Profils ANSYS : %d points, rayon d''echantillonnage %.2f mm (mi-entrefer MEC %.2f mm)\n',...
    numel(P0.th), Rmid*1e3, 0.5*(G.Rs+G.Rr)*1e3);

%% 2. Points de fonctionnement apparies (vitesses des snapshots)
n0=1500.212916; s0=max((M.ns*60-n0)/(M.ns*60),1e-4);   % a vide (s<=0 -> 1e-4)
n1=1469.772776; s1=(M.ns*60-n1)/(M.ns*60);             % en charge
r0=mec.equivalent_circuit(ctx,s0,ctx.Xm0);
r1=mec.equivalent_circuit(ctx,s1,ctx.Xm0);
fprintf('Points MEC apparies : a vide s=%.4g (I1=%.2f A) ; en charge s=%.4f (I1=%.2f A, I2=%.2f A)\n\n',...
    s0,r0.I1,s1,r1.I1,r1.I2);

%% 3. Champs MEC (maillage nc=6, nr=3)
nc=6; nr=3;
[me0,Se0]=solve_field(M,G,W,BH,ctx,r0,nc,nr);
[me1,Se1]=solve_field(M,G,W,BH,ctx,r1,nc,nr);
F0=gap_profile(me0,Se0,M);
F1=gap_profile(me1,Se1,M);

%% 4. Recalage de phase sur le fondamental + comparaison Br(theta)
[F0,c0m]=align_profile(F0,P0,p);
[F1,c1m]=align_profile(F1,P1,p);
c0a=fourier_uni(P0.Br,P0.th,p); c1a=fourier_uni(P1.Br,P1.th,p);

fprintf('--- Fondamental d''entrefer Bg1 (mi-entrefer) ---\n');
fprintf('%12s | %8s %8s %8s\n','essai','MEC [T]','EF [T]','ecart');
fprintf('%12s | %8.3f %8.3f %+7.1f%%\n','a vide',abs(c0m),abs(c0a),100*(abs(c0m)-abs(c0a))/abs(c0a));
fprintf('%12s | %8.3f %8.3f %+7.1f%%\n','en charge',abs(c1m),abs(c1a),100*(abs(c1m)-abs(c1a))/abs(c1a));

% flux par pole (via la distribution de A) et composante tangentielle
Phip_a0=(max(P0.A)-min(P0.A))*M.L;  Phip_m0=(max(F0.Acum)-min(F0.Acum))*M.L;
Phip_a1=(max(P1.A)-min(P1.A))*M.L;  Phip_m1=(max(F1.Acum)-min(F1.Acum))*M.L;
fprintf('\n--- Flux par pole (crete de la distribution A) ---\n');
fprintf('a vide    : MEC %.2f mWb   EF %.2f mWb   (%+.1f %%)\n',Phip_m0*1e3,Phip_a0*1e3,100*(Phip_m0-Phip_a0)/Phip_a0);
fprintf('en charge : MEC %.2f mWb   EF %.2f mWb   (%+.1f %%)\n',Phip_m1*1e3,Phip_a1*1e3,100*(Phip_m1-Phip_a1)/Phip_a1);
fprintf('\n--- Composante tangentielle au mi-entrefer (limite M6 du MEC) ---\n');
fprintf('Bt RMS ANSYS : %.3f T (a vide), %.3f T (en charge) ; MEC : 0 (branches radiales)\n\n',...
    rms(P0.Bt),rms(P1.Bt));

%% 5. Figure : Br(theta) superposes
figure('Name','Br entrefer MEC vs ANSYS','Position',[60 60 1150 640]);
cases={ {P0,F0,'a vide (s\approx0)'} , {P1,F1,sprintf('en charge (s=%.4f)',s1)} };
for k=1:2
    Pk=cases{k}{1}; Fk=cases{k}{2}; lab=cases{k}{3};
    subplot(2,2,2*k-1);
    plot(Pk.th*180/pi,Pk.Br,'r-','LineWidth',0.9); hold on;
    stairs(Fk.th_plot*180/pi,Fk.Bcol,'b-','LineWidth',1.1);
    grid on; xlim([0 360]); xlabel('\theta mecanique [deg]'); ylabel('B_r [T]');
    title(['B_r(\theta) mi-entrefer — ',lab]); legend('EF (ANSYS)','MEC','Location','southwest');
    subplot(2,2,2*k);
    plot(Pk.th*180/pi,Pk.Br,'r-','LineWidth',1.1); hold on;
    stairs(Fk.th_plot*180/pi,Fk.Bcol,'b-','LineWidth',1.3);
    grid on; xlim([0 90]); xlabel('\theta mecanique [deg]'); ylabel('B_r [T]');
    title('zoom sur un pole'); legend('EF','MEC','Location','southwest');
end
saveas(gcf,fullfile(OUT,'entrefer_Br_MEC_vs_ANSYS.png'));

%% 6. Figure : spectres spatiaux de Br (ordres electriques nu = k_mec/p)
numax=29;  nu=1:2:numax;                     % harmoniques impairs
Am0=arrayfun(@(n)abs(fourier_nonuni(F0.Bcol,F0.th,F0.dth,n*p)),nu);
Aa0=arrayfun(@(n)abs(fourier_uni(P0.Br,P0.th,n*p)),nu);
Am1=arrayfun(@(n)abs(fourier_nonuni(F1.Bcol,F1.th,F1.dth,n*p)),nu);
Aa1=arrayfun(@(n)abs(fourier_uni(P1.Br,P1.th,n*p)),nu);
figure('Name','Spectres Br','Position',[80 80 1150 420]);
tt={'a vide','en charge'};
for k=1:2
    if k==1, Aa=Aa0; Am=Am0; else, Aa=Aa1; Am=Am1; end
    subplot(1,2,k);
    bar(nu,[Aa(:) Am(:)],1.0,'grouped'); grid on;
    set(gca,'YScale','log'); ylim([1e-4 2]);
    xlabel('ordre harmonique \nu (electrique)'); ylabel('|B_{r,\nu}| [T]');
    title(['spectre spatial de B_r — ',tt{k}]);
    legend('EF (ANSYS)','MEC'); xticks([1 5 7 11 13 17 19 23 25 29]);
end
saveas(gcf,fullfile(OUT,'entrefer_spectre_MEC_vs_ANSYS.png'));

fprintf('--- Harmoniques de denture (ordres electriques 23/25 = Ns/p -/+ 1) ---\n');
fprintf('%10s | %10s %10s | %10s %10s\n','essai','EF nu=23','MEC nu=23','EF nu=25','MEC nu=25');
fprintf('%10s | %10.3f %10.3f | %10.3f %10.3f\n','a vide',Aa0(nu==23),Am0(nu==23),Aa0(nu==25),Am0(nu==25));
fprintf('%10s | %10.3f %10.3f | %10.3f %10.3f\n\n','en charge',Aa1(nu==23),Am1(nu==23),Aa1(nu==25),Am1(nu==25));

%% 7. Figure : distribution du flux A(theta)
figure('Name','Distribution du flux','Position',[100 100 1150 420]);
for k=1:2
    Pk=cases{k}{1}; Fk=cases{k}{2}; lab=cases{k}{3};
    subplot(1,2,k);
    plot(Pk.th*180/pi,(Pk.A-mean(Pk.A))*1e3,'r-','LineWidth',1.2); hold on;
    plot(Fk.th_plot*180/pi,(Fk.Acum-mean(Fk.Acum))*1e3,'b-','LineWidth',1.2);
    grid on; xlim([0 360]); xlabel('\theta mecanique [deg]');
    ylabel('A = \int B_r R d\theta  [mWb/m]');
    title(['distribution du flux — ',lab]); legend('EF (ANSYS)','MEC','Location','southwest');
end
saveas(gcf,fullfile(OUT,'flux_distribution_MEC_vs_ANSYS.png'));

%% 8. Cartes 2D MEC : |B| + direction du flux
map_field(me0,Se0,M,G,2.0,sprintf('MEC — |B| a vide (s=%.4g)',s0));
saveas(gcf,fullfile(OUT,'carte_B_MEC_avide.png'));
map_field(me1,Se1,M,G,2.2,sprintf('MEC — |B| en charge (s=%.4f)',s1));
saveas(gcf,fullfile(OUT,'carte_B_MEC_charge.png'));

%% 9. Sondes locales ANSYS vs B regionaux MEC
[Bts0,Bys0,Btr0,Byr0]=regional_max(me0,Se0);
[Bts1,Bys1,Btr1,Byr1]=regional_max(me1,Se1);
fprintf('--- Sondes locales ANSYS (.docx) vs maxima regionaux MEC [T] ---\n');
fprintf('%22s | %8s %8s %7s\n','region (essai)','ANSYS','MEC max','ecart');
prb('culasse stator (vide)',1.837,Bys0); prb('dent stator (vide)',1.633,Bts0);
prb('dent rotor (vide)',1.929,Btr0);     prb('culasse rotor* (vide)',1.674,Byr0);
prb('culasse stator (chg)',1.896,Bys1);  prb('dent stator (chg)',1.684,Bts1);
prb('dent rotor (chg)',1.870,Btr1);      prb('culasse rotor (chg)',1.577,Byr1);
fprintf('  (*sonde m4 ANSYS au sommet de culasse rotor ; B regional MEC = max sur la region)\n');
fprintf('\nFigures : entrefer_Br/_spectre/flux_distribution_MEC_vs_ANSYS.png, carte_B_MEC_avide/charge.png\n');
fprintf('=== Termine ===\n');

%% ================= fonctions locales =================
function P=read_profile(f)
    fid=fopen(f); hdr=fgetl(fid); fclose(fid);
    tk=regexp(hdr,'"([^"]+)"','tokens'); names=cellfun(@(c)c{1},tk,'UniformOutput',false);
    D=readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    gc=@(pat)D(:,find(contains(names,pat),1));
    dist=gc('Distance'); A=gc('Flux_Lines'); Br=gc('Br'); Bt=gc('Bt'); Bm=gc('Mag_B');
    P.C=dist(end)*1e-3;                               % circonference [m] (dernier pt = 2*pi)
    n=numel(dist)-1;                                  % dernier point = duplication de 0
    P.th=2*pi*dist(1:n)*1e-3/P.C; P.A=A(1:n); P.Br=Br(1:n); P.Bt=Bt(1:n); P.Bmag=Bm(1:n);
end

function [me,Se]=solve_field(M,G,W,BH,ctx,r,nc,nr)
    p=M.p; Nr=M.Nr;
    psi1=angle(r.I1c); psi2=angle(r.I2c);
    i3=sqrt(2)*r.I1*[cos(psi1);cos(psi1-2*pi/3);cos(psi1+2*pi/3)];
    Ibar=r.I2*(2*M.m*W.kw1*W.Nph)/Nr;
    thbar=2*pi*(0:Nr-1)'/Nr;
    ib=sqrt(2)*Ibar*cos(p*thbar-psi2);
    me=mec.mesh_refined(M,G,nc,i3,0,nr,ib);
    Se=mec.solve_mesh(me,BH,M.opt);
end

function F=gap_profile(me,Se,M)
    gi=me.gapfirst:numel(me.a);
    phig=accumarray(me.a(gi),Se.Phi(gi),[me.Ms 1]);   % flux d'entrefer par colonne stator
    F.th=me.th_s(:); F.dth=me.dth_s(:); F.Rg=me.Rg;
    F.Bcol=phig./(F.dth*me.Rg*M.L);
    F.Acum=cumsum(phig)/M.L;                          % distribution du flux [Wb/m]
end

function c=fourier_uni(y,th,k)
    n=numel(y); c=(2/n)*sum(y(:).*exp(-1i*k*th(:)));
end
function c=fourier_nonuni(y,th,dth,k)
    c=(1/pi)*sum(y(:).*exp(-1i*k*th(:)).*dth(:));
end

function [F,cm]=align_profile(F,P,p)
    ca=fourier_uni(P.Br,P.th,p);
    cm=fourier_nonuni(F.Bcol,F.th,F.dth,p);
    dphi=(angle(cm)-angle(ca))/p;                     % rotation qui recale le fondamental
    th2=mod(F.th+dphi,2*pi);
    [th2,ix]=sort(th2);
    F.th_plot=th2; F.Bcol=F.Bcol(ix); F.dth=F.dth(ix); F.th=F.th(ix);
    F.Acum=cumsum(F.Bcol.*F.dth)*F.Rg;                % [Wb/m], recalcule apres tri
    cm=fourier_nonuni(F.Bcol,F.th_plot,F.dth,p);
end

function map_field(me,Se,M,G,cmax,ttl)
    % carte |B| par cellule + direction du flux (quiver), style ANSYS
    Ms=me.Ms; Mr=me.Mr; Ls=me.Ls; Lr=me.Lr;
    rs=me.rs_edges; rr=me.rr_edges;
    nb=me.gapfirst-1;                                 % branches internes (hors entrefer)
    a=me.a(1:nb); b=me.b(1:nb);
    Bsig=zeros(nb,1); ok=me.A(1:nb)>0;
    Bsig(ok)=Se.Phi(ok)./me.A(ok);
    % coordonnees des noeuds
    [xn,yn]=node_xy(me,M,G);
    % accumulation vectorielle par noeud, SEPAREE radial/tangentiel (sinon le
    % moyennage sur les 4 branches divise chaque composante par 2)
    dx=xn(b)-xn(a); dy=yn(b)-yn(a); dl=hypot(dx,dy); dl(dl==0)=1;
    ux=dx./dl; uy=dy./dl;
    rn=hypot(xn(a),yn(a)); rn(rn==0)=1;
    crad=abs(ux.*xn(a)./rn + uy.*yn(a)./rn);          % |u . r_chapeau|
    isRad=crad>0.7;
    N=me.Nnodes;
    acc=@(m)deal( ...
        accumarray([a(m);b(m)],[Bsig(m).*ux(m);Bsig(m).*ux(m)],[N 1]), ...
        accumarray([a(m);b(m)],[Bsig(m).*uy(m);Bsig(m).*uy(m)],[N 1]), ...
        accumarray([a(m);b(m)],1,[N 1]));
    [vxR,vyR,cR]=acc(isRad);  [vxT,vyT,cT]=acc(~isRad);
    cR(cR==0)=1; cT(cT==0)=1;
    vx=vxR./cR+vxT./cT; vy=vyR./cR+vyT./cT; Bn=hypot(vx,vy);
    % polygones de cellules
    [XX,YY,CC]=deal([]);
    for la=1:Ls
        r1=rs(la); r2=rs(la+1);
        [X,Y]=cell_poly(me.th_s,me.dth_s,r1,r2);
        XX=[XX X]; YY=[YY Y]; CC=[CC Bn((la-1)*Ms+(1:Ms))']; %#ok<AGROW>
    end
    for la=1:Lr
        r1=rr(la); r2=rr(la+1);
        [X,Y]=cell_poly(me.th_r,me.dth_r,r1,r2);
        CC=[CC Bn(Ms*Ls+(la-1)*Mr+(1:Mr))']; XX=[XX X]; YY=[YY Y]; %#ok<AGROW>
    end
    figure('Name',ttl,'Position',[80 80 760 700]);
    patch(XX,YY,CC,'EdgeColor','none'); axis equal off; hold on;
    colormap(jet(11)); clim([0 cmax]); cb=colorbar; cb.Label.String='|B| [T]';
    % direction du flux (sous-echantillonnee)
    idx=(1:4:me.Nnodes)';
    q=0.005./max(Bn(idx),0.5);
    quiver(xn(idx),yn(idx),vx(idx).*q,vy(idx).*q,0,'k','LineWidth',0.4);
    title(ttl,'Interpreter','none');
end

function [x,y]=node_xy(me,M,G) %#ok<INUSD>
    Ms=me.Ms; Mr=me.Mr; Ls=me.Ls; Lr=me.Lr;
    rs_c=0.5*(me.rs_edges(1:end-1)+me.rs_edges(2:end));
    rr_c=0.5*(me.rr_edges(1:end-1)+me.rr_edges(2:end));
    x=zeros(me.Nnodes,1); y=x;
    for la=1:Ls
        id=(la-1)*Ms+(1:Ms);
        x(id)=rs_c(la)*cos(me.th_s); y(id)=rs_c(la)*sin(me.th_s);
    end
    for la=1:Lr
        id=Ms*Ls+(la-1)*Mr+(1:Mr);
        x(id)=rr_c(la)*cos(me.th_r); y(id)=rr_c(la)*sin(me.th_r);
    end
end

function [X,Y]=cell_poly(th,dth,r1,r2)
    th=th(:).'; dth=dth(:).';
    t1=th-dth/2; t2=th+dth/2;
    X=[r1*cos(t1); r1*cos(t2); r2*cos(t2); r2*cos(t1)];
    Y=[r1*sin(t1); r1*sin(t2); r2*sin(t2); r2*sin(t1)];
end

function [Bts,Bys,Btr,Byr]=regional_max(me,Se)
    % maxima de |B| de branche par region de fer
    nb=me.gapfirst-1; isFe=logical(me.iron(1:nb));
    a=me.a(1:nb); B=abs(Se.B(1:nb));
    Ms=me.Ms; Ls=me.Ls; Mr=me.Mr;
    nst=me.nr; nys=Ls-nst;  nrt=me.nr; nyr=me.Lr-nrt; %#ok<NASGU>
    laS=ceil(a/Ms); laS(a>Ms*Ls)=0;                   % couche stator (0 si rotor)
    aR=a-Ms*Ls; laR=ceil(aR/Mr); laR(a<=Ms*Ls)=0;
    Bts=max(B(isFe & laS>=1 & laS<=nst));
    Bys=max(B(isFe & laS>nst));
    Btr=max(B(isFe & laR>=1 & laR<=nrt));
    Byr=max(B(isFe & laR>nrt));
end

function prb(lab,ba,bm)
    fprintf('%22s | %8.3f %8.3f %+6.1f%%\n',lab,ba,bm,100*(bm-ba)/ba);
end
