%% RUN_FOURIER_GAP  -  P1 : entrefer HARMONIQUE (couronne de Laplace) — validation
%
%  Le couplage d'entrefer dent-a-dent (noyau gaussien sigma0) est remplace par
%  la solution EXACTE de Laplace dans la couronne d'air (mec.airgap_fourier,
%  M.opt.gap_fourier=1) : potentiels de surface projetes en serie de Fourier,
%  permeances harmoniques C_n=coth(n·ln(Rs/Rr)), D_n=1/sinh(n·ln(Rs/Rr)).
%  Carter, frange, composante tangentielle Bt et attenuation des harmoniques
%  de denture EMERGENT de la solution. C'est le correctif de fond du defaut
%  M6 identifie par l'etude critique (ETUDE_CRITIQUE_ANSYS_MEC.md, P1).
%
%  VALIDATION en 5 volets :
%    0. tests algebriques de la couronne (conservation, symetrie,
%       permeance homopolaire, invariance en r du couple MST) ;
%    1. a vide : Bg1 et Xm — gaussien vs FOURIER vs EF (~0,94 T / ~46 ohm) ;
%    2. profils Br/Bt au mi-entrefer vs ANSYS (cibles : nu=23/25 ~ 0,12/0,07
%       au lieu de 0,61/0,46 ; Bt_RMS ~ 0,13 T au lieu de 0) ;
%    3. en charge : couple MST harmonique vs EF (121,6 N.m) ;
%    4. ondulation de couple vs position (cible : << 88 % du gaussien).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
OUT=fileparts(mfilename('fullpath'));

fprintf('=== P1 : entrefer harmonique (couronne de Laplace) ===\n\n');
M=mec.machine_18_5kW(); ctx=mec.build_context(M);
G=ctx.G; W=ctx.W; BH=ctx.BH; p=M.p; Nr=M.Nr; mu0=4*pi*1e-7;
nc=6; nr=3;
thbar=2*pi*(0:Nr-1)'/Nr;
Rmid=0.5*(G.Rs+G.Rr);
thq=linspace(0,2*pi,2001); thq(end)=[];

%% 0. Tests algebriques de la couronne
Mf=M; Mf.opt.gap_fourier=1;
meT=mec.mesh_refined(Mf,G,nc,[0;0;0],0,nr,zeros(Nr,1));
AF=meT.gapF; MsF=AF.Ms; MrF=AF.Mr;
e1=max(abs(AF.Y*ones(MsF+MrF,1)))/max(abs(AF.Y(:)));         % conservation
e2=max(max(abs(AF.Y-AF.Y')))/max(abs(AF.Y(:)));              % symetrie
U0=[ones(MsF,1);zeros(MrF,1)];
Lam_num=sum(AF.Y(1:MsF,:)*U0);                               % permeance totale
Lam_ana=mu0*2*pi*M.L/log(G.Rs/G.Rr);
Us_t=cos(p*meT.th_s(:)); Ur_t=0.3*cos(p*(meT.th_r(:)-0.4));
[~,~,H1]=AF.field(Us_t,Ur_t,G.Rr+0.25*M.g,0);
[~,~,H2]=AF.field(Us_t,Ur_t,G.Rr+0.75*M.g,0);
T1=sum(H1.Brc.*H1.Btc+H1.Brs.*H1.Bts)*pi*(G.Rr+0.25*M.g)^2*M.L/mu0;
T2=sum(H2.Brc.*H2.Btc+H2.Brs.*H2.Bts)*pi*(G.Rr+0.75*M.g)^2*M.L/mu0;
fprintf('--- Tests algebriques ---\n');
fprintf('conservation (Y*1=0)        : %.1e\n',e1);
fprintf('symetrie                    : %.1e\n',e2);
fprintf('permeance homopolaire       : %.4e vs %.4e H (%+.2f %%)\n',...
    Lam_num,Lam_ana,100*(Lam_num-Lam_ana)/Lam_ana);
fprintf('invariance en r du couple   : %+.4e vs %+.4e (%+.1e rel)\n\n',...
    T1,T2,(T1-T2)/max(abs(T1),eps));

%% 1. A VIDE : Bg1 et Xm — gaussien vs Fourier vs EF
s0=1e-4; r0=mec.equivalent_circuit(ctx,s0,ctx.Xm0);
[i30,ib0]=inst_currents(M,W,r0);
% gaussien (reference interne)
tic;
meG=mec.mesh_refined(M,G,nc,i30,0,nr,ib0); SeG=mec.solve_mesh(meG,BH,M.opt);
tG=toc;
gi=meG.gapfirst:numel(meG.a);
phig=accumarray(meG.a(gi),SeG.Phi(gi),[meG.Ms 1]);
BcolG=phig./(meG.dth_s(:)*meG.Rg*M.L);
cG=(1/pi)*sum(BcolG.*exp(-1i*p*meG.th_s(:)).*meG.dth_s(:));
% Fourier
tic;
meF=mec.mesh_refined(Mf,G,nc,i30,0,nr,ib0); SeF=mec.solve_mesh(meF,BH,Mf.opt);
tF=toc;
[UsF,UrF]=surfU(meF,SeF);
[Br0,Bt0]=meF.gapF.field(UsF,UrF,G.Rs,thq);      % au rayon stator (pour Xm)
c0=fu(Br0,thq,p);
[BrM,BtM]=meF.gapF.field(UsF,UrF,Rmid,thq);      % au mi-entrefer (vs ANSYS)
cM=fu(BrM,thq,p);
Xm_of=@(Bg1) M.w*W.kw1*W.Nph*((2/pi)*Bg1*(pi*G.Ds/(2*p))*M.L)/sqrt(2)/r0.I1;
fprintf('--- 1. A VIDE (I0=%.2f A) --- [solve : gaussien %.1f s, Fourier %.1f s, %d it, res %.1e, conv %d]\n',...
    r0.I1,tG,tF,SeF.iter,SeF.res,SeF.converged);
fprintf('%22s | %10s %10s %8s\n','','gaussien','FOURIER','EF');
fprintf('%22s | %10.3f %10.3f %8.3f\n','Bg1 (mi-entrefer) [T]',abs(cG),abs(cM),0.942);
fprintf('%22s | %10.1f %10.1f %8.0f\n','Xm [ohm]',Xm_of(abs(cG)),Xm_of(abs(c0)),46);

%% 2. Profils vs ANSYS (a vide et en charge)
P0=read_profile(fullfile(ROOT,'transitoire','a vide','Calculator Expressions Plot 1.tab'));
P1p=read_profile(fullfile(ROOT,'transitoire','en charge','Calculator Expressions Plot 1.tab'));

% --- charge : point apparie au snapshot ANSYS ---
n1=1469.772776; s1=(M.ns*60-n1)/(M.ns*60);
r1=mec.equivalent_circuit(ctx,s1,ctx.Xm0);
[i31,ib1]=inst_currents(M,W,r1);
meF1=mec.mesh_refined(Mf,G,nc,i31,0,nr,ib1); SeF1=mec.solve_mesh(meF1,BH,Mf.opt);
[UsF1,UrF1]=surfU(meF1,SeF1);
[BrM1,BtM1]=meF1.gapF.field(UsF1,UrF1,Rmid,thq);

fprintf('\n--- 2. Contenu harmonique de Br au mi-entrefer ---\n');
fprintf('%16s | %8s %8s %8s | %8s %8s %8s\n','',...
    'EF vide','FOURIER','gauss.','EF chg','FOURIER','gauss.');
nus=[1 5 7 23 25];
Ag=harm_tab(BcolG,meG.th_s(:),meG.dth_s(:),nus*p);          % gaussien a vide
% gaussien en charge (reference interne)
meG1=mec.mesh_refined(M,G,nc,i31,0,nr,ib1); SeG1=mec.solve_mesh(meG1,BH,M.opt);
phig1=accumarray(meG1.a(meG1.gapfirst:end),SeG1.Phi(meG1.gapfirst:end),[meG1.Ms 1]);
BcolG1=phig1./(meG1.dth_s(:)*meG1.Rg*M.L);
Ag1=harm_tab(BcolG1,meG1.th_s(:),meG1.dth_s(:),nus*p);
for q=1:numel(nus)
    fprintf('%14s%2d | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f\n','nu=',nus(q),...
        abs(fu(P0.Br,P0.th,nus(q)*p)),abs(fu(BrM,thq,nus(q)*p)),Ag(q),...
        abs(fu(P1p.Br,P1p.th,nus(q)*p)),abs(fu(BrM1,thq,nus(q)*p)),Ag1(q));
end
fprintf('%16s | %8.3f %8.3f %8s | %8.3f %8.3f %8s\n','Bt RMS [T]',...
    rms(P0.Bt),rms(Bt0*0+BtM),'0',rms(P1p.Bt),rms(BtM1),'0');

%% 3. Couple en charge (MST harmonique)
Tmst=meF1.gapF.torque(UsF1,UrF1);
fprintf('\n--- 3. COUPLE en charge (s=%.4f) ---\n',s1);
fprintf('MST harmonique (Fourier) = %+8.1f N.m   (circuit %.1f, EF %.1f)\n',...
    Tmst,r1.Tem,121.6);

%% 4. Ondulation de couple vs position rotorique
Npos=19; th0=linspace(0,2*(2*pi/M.Ns),Npos);
psi1=angle(r1.I1c); psi2=angle(r1.I2c);
Ibar=r1.I2*(2*M.m*W.kw1*W.Nph)/Nr;
Tpos=zeros(Npos,1);
tic;
for k=1:Npos
    ph=p*th0(k);
    i3k=sqrt(2)*r1.I1*[cos(psi1+ph);cos(psi1+ph-2*pi/3);cos(psi1+ph+2*pi/3)];
    ibk=sqrt(2)*Ibar*cos(p*thbar-psi2);
    mek=mec.mesh_refined(Mf,G,nc,i3k,th0(k),nr,ibk);
    Sek=mec.solve_mesh(mek,BH,Mf.opt);
    [Usk,Urk]=surfU(mek,Sek);
    Tpos(k)=mek.gapF.torque(Usk,Urk);
end
tR=toc;
ond=100*(max(Tpos)-min(Tpos))/abs(mean(Tpos));
fprintf('\n--- 4. ONDULATION (%d positions, %.0f s) ---\n',Npos,tR);
fprintf('T moyen = %.1f N.m   ondulation = %.1f %%   (gaussien : ~88 %% ; physique ~1-5 %%)\n',...
    mean(Tpos),ond);

%% Figures
figure('Name','P1 profils','Position',[60 60 1150 640]);
subplot(2,2,1);
plot(P0.th*180/pi,P0.Br,'r-','LineWidth',0.9); hold on;
plot(al(thq,BrM,P0,p)*180/pi,BrM,'b.','MarkerSize',3);
grid on; xlim([0 90]); xlabel('\theta [deg]'); ylabel('B_r [T]');
title('B_r mi-entrefer, a vide (zoom 1 pole)'); legend('EF','MEC Fourier');
subplot(2,2,2);
plot(P0.th*180/pi,P0.Bt,'r-','LineWidth',0.9); hold on;
plot(al(thq,BtM,P0,p)*180/pi,BtM,'b.','MarkerSize',3);
grid on; xlim([0 90]); xlabel('\theta [deg]'); ylabel('B_t [T]');
title('B_t mi-entrefer, a vide — ABSENT du modele gaussien'); legend('EF','MEC Fourier');
subplot(2,2,3);
plot(P1p.th*180/pi,P1p.Br,'r-','LineWidth',0.9); hold on;
plot(al(thq,BrM1,P1p,p)*180/pi,BrM1,'b.','MarkerSize',3);
grid on; xlim([0 90]); xlabel('\theta [deg]'); ylabel('B_r [T]');
title('B_r mi-entrefer, en charge'); legend('EF','MEC Fourier');
subplot(2,2,4);
plot(th0*180/pi,Tpos,'b-o','LineWidth',1.3); hold on; yline(mean(Tpos),'k--');
grid on; xlabel('position rotor [deg mec]'); ylabel('T [N.m]');
title(sprintf('Couple MST harmonique (ondulation %.1f %%)',ond));
saveas(gcf,fullfile(OUT,'fourier_gap_validation.png'));
fprintf('\nFigure enregistree : fourier_gap_validation.png\n=== Termine ===\n');

%% ================= fonctions locales =================
function [i3,ib]=inst_currents(M,W,r)
    psi1=angle(r.I1c); psi2=angle(r.I2c);
    i3=sqrt(2)*r.I1*[cos(psi1);cos(psi1-2*pi/3);cos(psi1+2*pi/3)];
    Ibar=r.I2*(2*M.m*W.kw1*W.Nph)/M.Nr;
    ib=sqrt(2)*Ibar*cos(M.p*(2*pi*(0:M.Nr-1)'/M.Nr)-psi2);
end
function [Us,Ur]=surfU(me,Se)
    Us=Se.U(me.gapF.ids(1:me.Ms)); Ur=Se.U(me.gapF.ids(me.Ms+1:end));
end
function c=fu(y,th,k)
    n=numel(y); c=(2/n)*sum(y(:).*exp(-1i*k*th(:)));
end
function A=harm_tab(Bcol,th,dth,ks)
    A=arrayfun(@(k)abs((1/pi)*sum(Bcol.*exp(-1i*k*th).*dth)),ks);
end
function th2=al(thq,Bmec,Pans,p)
%  recalage de phase du MEC sur le fondamental ANSYS (affichage seulement)
    ca=(2/numel(Pans.Br))*sum(Pans.Br(:).*exp(-1i*p*Pans.th(:)));
    cm=(2/numel(Bmec))*sum(Bmec(:).*exp(-1i*p*thq(:)));
    th2=mod(thq(:)+(angle(cm)-angle(ca))/p,2*pi);
end
function P=read_profile(f)
    fid=fopen(f); hdr=fgetl(fid); fclose(fid);
    tk=regexp(hdr,'"([^"]+)"','tokens'); names=cellfun(@(c)c{1},tk,'UniformOutput',false);
    D=readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    gc=@(pat)D(:,find(contains(names,pat),1));
    dist=gc('Distance'); Br=gc('Br'); Bt=gc('Bt');
    C=dist(end)*1e-3; n=numel(dist)-1;
    P.th=2*pi*dist(1:n)*1e-3/C; P.Br=Br(1:n); P.Bt=Bt(1:n);
end
