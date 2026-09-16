%% RUN_Z8_CAVITY_LOAD - carte de cavite rotorique EN CHARGE : effet du
%  courant de barre dans l'element de cavite sur le point nominal.
%
%  La carte Q_r des cavites rotoriques est calculee a courant de barre nul
%  (parois a un seul potentiel). En charge, la barre porte un courant I_z :
%  les deux dents adjacentes different de I_z (Ampere) et le champ de la
%  cavite est la somme du champ de Laplace (Q_r) et du champ particulier du
%  courant uniforme dans la barre. L'element recoit un vecteur de flux
%  source par ampere (cavity_src_nO16.mat, cavity_graded.py), ajoute au
%  second membre de la condensation (eq. 7) et au bilan de flux du reseau
%  (mec.solve_network, AG.f), sans modifier Q.
%
%  TROIS VARIANTES au point nominal s = 0,0188, fermeture cavites (33,16),
%  vrillage harmonique neutralise :
%    A  chaine du manuscrit (Table 8) : X_m et R_fe mis a jour sur le champ
%       A VIDE du reseau au courant magnetisant |I_m| (mec.magnetizing) ;
%       le courant de barre n'y entre pas, la source est sans effet.
%    B  point fixe sur le champ EN CHARGE (FMM stator I1 + onde de barres I2,
%       comme la chaine de champ de RUN_Z2/Fig. 6), sans source de cavite :
%       X_m = E1(champ charge)/|I1 - I2|, R_fe des pertes fer du champ charge.
%    C  comme B, avec la source de cavite au courant de barre du reseau
%       (I_z tire de l'onde de FMM de barres, point fixe).
%  L'effet de la carte de cavite en charge est C - B ; B - A mesure, a part,
%  ce que l'excitation en charge du reseau change a la branche magnetisante.
%  Sortie : Z8_cavity_load_out.txt, Z8_cavity_load.mat
clear; clc; t0=tic;
if isfile('Z8_cavity_load_out.txt'), delete('Z8_cavity_load_out.txt'); end
diary('Z8_cavity_load_out.txt'); diary on;
ROOT=fullfile('..','..','reference','ANSYS_18_5kW');
Nh=8192; s_ch=0.0188; Im_ref=0.2; nT=33; nO=16;
M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
ctx=mec.build_context(M); G=ctx.G; W=ctx.W; ctx.M.opt.skew_harm=0;
p=M.p; Rm=0.5*(G.Rs+G.Rr); Nr=M.Nr; Ns=M.Ns;
cav=load(sprintf('cavity_nO%d.mat',nO)); src=load(sprintf('cavity_src_nO%d.mat',nO)); cav.sr=src.sr;
A1=mec.airgap_dtn_tooth_cav(M,G,0,nT,nO,Nh,'p1',cav); ctx.AG=A1;
RmS=mec.magnetizing(ctx,Im_ref); ctx.Xm0=RmS.Xm;
fprintf('=== Z8 : carte de cavite rotorique en charge (courant de barre dans l''element) ===\n\n');
fprintf('  fermeture cavites (%d,%d), N_h = %d, base p1, skew_harm = 0, s = %.4f | X_m0 = %.3f ohm (%.0f s)\n',nT,nO,Nh,s_ch,ctx.Xm0,toc(t0));
fprintf('  vecteur source (Wb/A) : colonnes %s ; paroi %.3e ; somme %.1e\n',mat2str(src.sr(1:nO).',3),src.sr(nO+1),sum(src.sr));
fprintf('  part du champ particulier |Phi_J|max / |Q ramp|max = %.1e (parois quasi paralleles : la rampe domine)\n', ...
    max(abs(src.PhiJ))/max(abs(src.sr-src.PhiJ)));

%% ---- A : chaine du manuscrit ----------------------------------------------
rA=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
PUB=[115.53 18.675 432.7 70.029]; got=[rA.Tem rA.I1 rA.Pcu_r ctx.Xm0]; okA=all(abs(got-PUB)./PUB<5e-4);
fprintf('\n  A (chaine, Table 8) : T %.4f | I1 %.4f | Pcu_r %.2f | X_m0 %.3f  -> garde Table 8 %s\n',got,tern(okA,'PASSEE','ECHOUEE'));

%% ---- convention de signe : potentiels de dents rotoriques contre l'onde de FMM de barres ----
[Fs,Fr]=mmf_from_phasors(rA.I1c,rA.I2c,M,W);
SL=mec.solve_network(ctx.net,G,ctx.BH,A1,Fs,Fr,M.opt);
dU=SL.Ur([2:end 1])-SL.Ur; dF=Fr([2:end 1])-Fr;
slope=(dF.'*dU)/(dF.'*dF); resid=norm(dU-slope*dF)/norm(dU);
fprintf('\n  convention : U_dent(i+1)-U_dent(i) = %.4f x [F_r(i+1)-F_r(i)] + residu %.1f %% (chute de fer du reseau)\n',slope,100*resid);
sgn=sign(slope);
fprintf('  Ampere dans la cavite : U_droite - U_gauche = -I_z  =>  I_z = %+d x [F_r(i+1)-F_r(i)] ; |I_z| max %.1f A, rms %.1f A\n', ...
    -sgn,max(abs(dF)),sqrt(mean(dF.^2)));
fprintf('  (rappel : courant de barre du circuit I_bar1 = %.1f A rms ; onde echantillonnee crete %.1f A)\n',rA.Ibar1,max(abs(dF)));

%% ---- B et C : point fixe sur le champ en charge ---------------------------
P1=read_profile(fullfile(ROOT,'transitoire','en charge','Calculator Expressions Plot 1.tab'));
RES=struct(); modes={'load','load+src'}; lab={'B (champ charge, sans source)','C (champ charge + source de cavite)'};
for k=1:2
    [r,info]=ec_loadfield(ctx,s_ch,ctx.Xm0,modes{k},sgn);
    RES.(sprintf('m%d',k))=r; RES.(sprintf('i%d',k))=info;
    fprintf('\n  %s : %d iterations (dX %.1e)\n',lab{k},r.iterFP,info.dX);
    fprintf('    T %.4f N.m | I1 %.4f A | I2 %.4f A | Pcu_r %.2f W | Pfe %.2f W | cosphi %.4f | eta %.4f\n',r.Tem,r.I1,r.I2,r.Pcu_r,r.Pfe,r.cosphi,r.eta);
    fprintf('    X_m %.4f ohm | E1 circuit %.3f V | E1 champ %.3f V | |I_m| %.4f | |I1-I2| %.4f | Bg1 (flux/dent) %.4f T\n', ...
        r.Xm,r.E1,info.E1f,r.Im,info.Iem,info.Bg1);
    %  champ au mi-entrefer, en charge, recale sur la phase du fondamental (comme Z4)
    AG=info.AG; thq=P1.th(:).';
    [Br,Bt]=AG.field(info.Usurf,Rm,thq); Br=Br(:).'; Bt=Bt(:).';
    cf=@(y,th,kk)(2/numel(y))*sum(y(:).*exp(-1j*kk*th(:)));
    ca=cf(P1.Br,P1.th,p); phA=atan2(-imag(ca),real(ca)); cm=cf(Br,thq,p); phM=atan2(-imag(cm),real(cm));
    dphi=(phM-phA)/p; [Br2,Bt2]=AG.field(info.Usurf,Rm,thq+dphi); Br2=Br2(:).'; Bt2=Bt2(:).';
    rmse=sqrt(mean((Br2-P1.Br(:).').^2)); rref=sqrt(mean(P1.Br.^2));
    fprintf('    champ mi-entrefer en charge : Bg1 %.4f T (EF %.4f) | Bt rms %.4f (EF %.4f) | RMSE Br %.4f T = %.1f %% du rms EF\n', ...
        abs(cm),abs(ca),sqrt(mean(Bt2.^2)),sqrt(mean(P1.Bt.^2)),rmse,100*rmse/rref);
    RES.(sprintf('f%d',k))=[abs(cm) sqrt(mean(Bt2.^2)) rmse 100*rmse/rref];
    if k==2
        fprintf('    potentiel d''ouverture induit par la source : max |uO| %.2f A (potentiels de dents rotor : max %.1f A) ; |f| max %.2e Wb\n', ...
            max(abs(AG.uO)),max(abs(info.Usurf(Ns+1:end))),max(abs(AG.f)));
    end
end
rB=RES.m1; rC=RES.m2;
fprintf('\n  ---- SYNTHESE (point nominal, cavites (33,16)) ----\n');
fprintf('  %-34s %10s %10s %10s %10s %10s\n','variante','T (N.m)','I1 (A)','Pcu_r (W)','X_m (ohm)','E1 (V)');
fprintf('  %-34s %10.3f %10.4f %10.2f %10.3f %10.2f\n','A chaine (Table 8)',rA.Tem,rA.I1,rA.Pcu_r,rA.Xm,rA.E1);
fprintf('  %-34s %10.3f %10.4f %10.2f %10.3f %10.2f\n','B champ charge sans source',rB.Tem,rB.I1,rB.Pcu_r,rB.Xm,rB.E1);
fprintf('  %-34s %10.3f %10.4f %10.2f %10.3f %10.2f\n','C champ charge + source',rC.Tem,rC.I1,rC.Pcu_r,rC.Xm,rC.E1);
fprintf('  %-34s %10s %10s %10s\n','reference EF','121.63','19.72','488.5');
fprintf('  effet de la source (C - B)        : T %+.3f %% | I1 %+.3f %% | Pcu_r %+.3f %% | X_m %+.3f %%\n', ...
    100*(rC.Tem-rB.Tem)/rB.Tem,100*(rC.I1-rB.I1)/rB.I1,100*(rC.Pcu_r-rB.Pcu_r)/rB.Pcu_r,100*(rC.Xm-rB.Xm)/rB.Xm);
fprintf('  excitation en charge (B - A)      : T %+.3f %% | I1 %+.3f %% | Pcu_r %+.3f %% | X_m %+.3f %%\n', ...
    100*(rB.Tem-rA.Tem)/rA.Tem,100*(rB.I1-rA.I1)/rA.I1,100*(rB.Pcu_r-rA.Pcu_r)/rA.Pcu_r,100*(rB.Xm-rA.Xm)/rA.Xm);
fprintf('  total (C - A, contre la Table 8)  : T %+.3f %% | I1 %+.3f %% | Pcu_r %+.3f %% | X_m %+.3f %%\n', ...
    100*(rC.Tem-rA.Tem)/rA.Tem,100*(rC.I1-rA.I1)/rA.I1,100*(rC.Pcu_r-rA.Pcu_r)/rA.Pcu_r,100*(rC.Xm-rA.Xm)/rA.Xm);
fprintf('  incertitude de pavage (Z6, 8 pavages cavites) : couple 0.49 %%, I1 1.02 %%\n');
%  les operateurs (AG) et solutions de champ (S) portent des centaines de Mo : on ne sauve que les scalaires
for k=1:2, RES.(sprintf('i%d',k))=rmfield(RES.(sprintf('i%d',k)),{'AG','S','Usurf'}); end
save('Z8_cavity_load.mat','rA','RES','slope','resid','sgn','src');
fprintf('\n  duree %.0f s\n=== Z8 termine ===\n',toc(t0));
diary off;

% ======================================================================
function [Fs,Fr]=mmf_from_phasors(I1c,I2c,M,W)
    p1=angle(I1c);
    i3=sqrt(2)*abs(I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
    Fs=W.slotMMF(i3);
    p2=angle(I2c); tb=2*pi*(0:M.Nr-1)/M.Nr;
    Fu=cumsum(cos(M.p*tb-p2).'); Fu=Fu-mean(Fu);
    c1=(2/M.Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
    Fr=-Fu*((3/2)*(4/pi)*(W.kw1*W.Nph/(2*M.p))*sqrt(2)*abs(I2c)/abs(c1));
end

function [res,info]=ec_loadfield(ctx,s,Xm_init,mode,sgn)
%  Copie de mec.equivalent_circuit dont la mise a jour de X_m et R_fe est
%  faite sur le champ EN CHARGE (FMM stator + onde de barres), avec ou sans
%  la source de cavite. Tout le reste est identique.
M=ctx.M; G=ctx.G; W=ctx.W; m=M.m; Uph=M.Uph; w=M.w; ws=w/M.p; Rs=ctx.Rs;
Xs_leak=ctx.Lk.Xs_leak;
Cg=mec.cage(M,G,W,ctx.Lk,s); Rr=Cg.Rr; Xr=Cg.Xr; sigd_r=Cg.sigd_r;
H=ctx.H; nH=numel(H.nu); snu=zeros(nH,1); Rrn=snu; Xrn=snu; frn=snu;
for k=1:nH
    if H.dir(k)>0, snu(k)=1-H.nu(k)*(1-s); else, snu(k)=1+H.nu(k)*(1-s); end
    Ck=mec.cage(M,G,W,ctx.Lk,max(abs(snu(k)),1e-3),H.nu(k),H.kw(k));
    Rrn(k)=Ck.Rr; Xrn(k)=Ck.Xr;
    srn=max(sin(pi*H.nu(k)*M.p/M.Nr)^2,1e-9); Req=Ck.Rring/(2*srn); frn(k)=Req/(Ck.Rbar+Req);
end
Xm=Xm_init; Rfe=50*ctx.Xm0; tolX=M.opt.sat_tol; itmax=120;   % point fixe plus long que la chaine (convergence lente en champ charge)
Ns=M.Ns; p=M.p; th=2*pi*(0:Ns-1)/Ns; tp=pi*G.Ds/(2*p);
I1=0; E1v=0; Im=0; I2=0; Zhk=zeros(nH,1); dX=inf; AG=ctx.AG; S=[]; E1f=NaN; Bg1=NaN;
for it=1:itmax
    Xs=Xs_leak; Xrp=Xr+sigd_r*Xm; Zrot=Rr/s+1i*Xrp;
    Yp=1/Zrot+1/(1i*Xm)+1/Rfe; Zp=1/Yp;
    Zh=0;
    for k=1:nH
        Xmn=Xm*H.sig(k);
        if Xmn<=0||abs(snu(k))<1e-6, Zhk(k)=0; continue; end
        Zhk(k)=1/(1/(1i*Xmn)+1/(Rrn(k)/snu(k)+1i*Xrn(k))); Zh=Zh+Zhk(k);
    end
    Ztot=Rs+1i*Xs+Zh+Zp; I1=Uph/Ztot; E1v=I1*Zp; Im=E1v/(1i*Xm); I2=E1v/Zrot;
    %  ---- mise a jour sur le champ EN CHARGE ----
    [Fs,Fr]=mmf_from_phasors(I1,I2,M,W);
    AG=ctx.AG;
    if strcmp(mode,'load+src')
        Iz=-sgn*(Fr([2:end 1])-Fr);            % Ampere : U_droite - U_gauche = -I_z
        AG=AG.set_source(Iz);
    end
    S=mec.solve_network(ctx.net,G,ctx.BH,AG,Fs,Fr,M.opt);
    Bg=S.Bgap_avg_i(:).'; Bg1=abs((2/Ns)*sum(Bg.*exp(-1j*p*th)));
    E1f=M.w*W.kw1*W.Nph*((2/pi)*Bg1*tp*M.L)/sqrt(2);
    Iem=abs(I1-I2);
    Xm_new=E1f/Iem;
    PfeS=mec.iron_losses(M,G,S,M.f); Pfe=PfeS.total;
    Rfe_new=m*abs(E1v)^2/max(Pfe,1);
    dX=abs(Xm_new-Xm)/max(Xm,eps); dR=abs(Rfe_new-Rfe)/max(Rfe,eps);
    Xm=Xm+0.6*(Xm_new-Xm); Rfe=Rfe+0.6*(Rfe_new-Rfe);
    if dX<tolX && dR<1e-2, break; end
end
I1m=abs(I1); I2m=abs(I2); Imm=abs(Im);
cosphi=cos(angle(I1)); Pin=m*Uph*real(I1);
Pag=m*I2m^2*Rr/s; Tem1=Pag/ws;
Pag_h=zeros(nH,1); T_h=zeros(nH,1);
for k=1:nH, Pag_h(k)=m*I1m^2*real(Zhk(k)); T_h(k)=Pag_h(k)*H.nu(k)*M.p/M.w*H.dir(k); end
Tem_h=sum(T_h); Pcu_r_h=sum(Pag_h.*snu); Pmech_h=sum(Pag_h.*(1-snu));
Tem=Tem1+Tem_h; Pcu_s=m*Rs*I1m^2; Pcu_r=m*Rr*I2m^2+Pcu_r_h; Pfe=m*abs(E1v)^2/Rfe;
Pmech_gross=Pag*(1-s)+Pmech_h; Pfw=mec.mech_losses(M,s); Padd=mec.stray_losses(M,I2m);
Pout=Pmech_gross-Pfw-Padd; eta=max(Pout,0)/Pin;
res.I1=I1m; res.I2=I2m; res.Im=Imm; res.I1c=I1; res.I2c=I2; res.E1=abs(E1v); res.cosphi=cosphi;
res.Tem=Tem; res.Tem1=Tem1; res.Tem_h=Tem_h; res.Pcu_s=Pcu_s; res.Pcu_r=Pcu_r; res.Pfe=Pfe; res.Pout=Pout; res.Pin=Pin; res.eta=eta;
res.Xm=Xm; res.Rfe=Rfe; res.iterFP=it;
info.E1f=E1f; info.Bg1=Bg1; info.Iem=abs(I1-I2); info.AG=AG; info.Usurf=S.Usurf; info.dX=dX; info.S=S;
end

function P=read_profile(f)
    fid=fopen(f); hdr=fgetl(fid); fclose(fid);
    tk=regexp(hdr,'"([^"]+)"','tokens'); names=cellfun(@(c)c{1},tk,'UniformOutput',false);
    D=readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    gc=@(pat)D(:,find(contains(names,pat),1));
    dist=gc('Distance'); A=gc('Flux_Lines'); Br=gc('Br'); Bt=gc('Bt');
    C=dist(end)*1e-3; n=numel(dist)-1;
    P.th=2*pi*dist(1:n)*1e-3/C; P.A=A(1:n); P.Br=Br(1:n); P.Bt=Bt(1:n);
end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
