%% RUN_Z5_SWEEP_CAV - chaine MAS complete avec la fermeture A CAVITES (33,16),
%  vrillage harmonique neutralise.  Copie de RUN_B10_B1_SKEWOFF (source de
%  B10_b1_skewoff.mat, caracteristique T(s), I(s) de la figure des
%  caracteristiques) ou seule la construction de l'operateur change :
%      A1 = mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1')
%  devient
%      cav = load('cavity_nO16.mat');
%      A1  = mec.airgap_dtn_tooth_cav(M,G,0,33,16,Nh,'p1',cav);
%  Sortie : Z5_sweep_cav.mat (memes variables sl, T, I, r, r0, rb que
%  B10_b1_skewoff.mat), transcript Z5_sweep_cav_out.txt.
clear; clc; t0=tic;
if isfile('Z5_sweep_cav_out.txt'), delete('Z5_sweep_cav_out.txt'); end
diary('Z5_sweep_cav_out.txt'); diary on;

M = mec.machine_18_5kW();
M.opt.skew_harm = 0;                       % configuration loyale (comme B10)
M.opt.verbose = false;
ctx = mec.build_context(M); G = ctx.G;
ctx.M.opt.skew_harm = 0;
nT=33; nO=16; Nh=8192; s_ch=0.0188;
cav=load(sprintf('cavity_nO%d.mat',nO));
A1 = mec.airgap_dtn_tooth_cav(M,G,0,nT,nO,Nh,'p1',cav);
ctx.AG = A1; ctx.Xm0 = mec.magnetizing(ctx,0.2).Xm;

fprintf('=== Z5 : chaine MAS complete, fermeture a CAVITES, vrillage harmonique NEUTRALISE ===\n');
fprintf('  pavage nT=%d nO=%d | N_h = %d | base P1 | skew_harm = 0 | cavites cavity_nO%d.mat\n',nT,nO,Nh,nO);
fprintf('  operateur %dx%d | X_m0 = %.3f ohm  (%.0f s)\n',size(A1.Y,1),size(A1.Y,2),ctx.Xm0,toc(t0));
fprintf('  references EF : famille TRANSITOIRE (121.63 N.m / 19.72 A ; 8.499 A a vide ; 104.31 N.m au calage)\n\n');

%% ---- 1. schema equivalent -------------------------------------------
r = mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
fprintf('  ---- 1. SCHEMA EQUIVALENT (point nominal s = %.4f) ----\n',s_ch);
fprintf('  %-24s %12s\n','parametre','cav (33,16)');
pr=@(l,a)fprintf('  %-24s %12.4f\n',l,a);
pr('Rs [ohm]',ctx.Rs);
pr('Rr'' [ohm]',r.Rr);
pr('Xsigma_s [ohm]',ctx.Lk.Xs_leak);
pr('Xsigma_r [ohm]',r.Xr);
pr('Xm sature [ohm]',r.Xm);
pr('Rfe [ohm]',r.Rfe);
fprintf('  X_m NON sature : %.3f ohm\n',ctx.Xm0);

%% ---- 2. point nominal ------------------------------------------------
fprintf('\n  ---- 2. POINT NOMINAL ----\n');
fprintf('  %-24s %12s %12s %10s\n','grandeur','cav (33,16)','EF','ecart');
q=@(l,a,f)fprintf('  %-24s %12.4f %12.4f %9.2f %%\n',l,a,f,100*(a-f)/f);
q('couple [N.m]',r.Tem,121.63);
q('courant I1 [A]',r.I1,19.72);
q('pertes fer [W]',r.Pfe,232.6);
fprintf('  %-24s %12.4f\n','cos(phi)',r.cosphi);
fprintf('  %-24s %12.4f\n','rendement',r.eta);
fprintf('  pertes Joule rotor [W]   %12.4f   (dont harmoniques %.4f) ; EF 488.5\n',r.Pcu_r,r.Pcu_r_h);

%% ---- 3. a vide et calage ---------------------------------------------
r0 = mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
rb = mec.equivalent_circuit(ctx,1.0 ,ctx.Xm0);
kbar = (2*M.m*ctx.W.kw1*ctx.W.Nph)/M.Nr;
fprintf('\n  ---- 3. A VIDE ----\n');
q('  I magnetisant [A]',r0.I1,8.499);
fprintf('  %-24s %12.4f   (EF InducedVoltage 2D : 381.5 V, contient la chute de fuite 2D)\n','  f.e.m. entrefer [V]',r0.E1);
q('  pertes fer [W]',r0.Pfe,249.3);
fprintf('\n  ---- 4. ROTOR BLOQUE (diagnostic) ----\n');
q('  couple [N.m]',rb.Tem,104.31);
q('  I1 [A]',rb.I1,104.4);
fprintf('  I barre [A]              %12.4f\n',rb.I2*kbar);
fprintf('  pertes Joule rotor [W]   %12.4f   (dont harmoniques %.4f)\n',rb.Pcu_r,rb.Pcu_r_h);

%% ---- 5. champs a mi-entrefer -----------------------------------------
fprintf('\n  ---- 5. CHAMPS A MI-ENTREFER (grille 2000 pts, phi = 0) ----\n');
Rm=0.5*(G.Rs+G.Rr); thq=linspace(0,2*pi,2001); thq(end)=[];
Rnl=mec.magnetizing(ctx,r0.Im);
[Br0,Bt0]=A1.field(Rnl.S.Usurf,Rm,thq); Br0=Br0(:).'; Bt0=Bt0(:).';
p1=angle(r.I1c);
i3=sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
Fs=ctx.W.slotMMF(i3);
p2=angle(r.I2c); tb=2*pi*(0:M.Nr-1)/M.Nr;
Fu=cumsum(cos(M.p*tb-p2).'); Fu=Fu-mean(Fu);
c1=(2/M.Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
Fr=-Fu*((3/2)*(4/pi)*(ctx.W.kw1*ctx.W.Nph/(2*M.p))*sqrt(2)*abs(r.I2c)/abs(c1));
SL=mec.solve_network(ctx.net,G,ctx.BH,A1,Fs,Fr,M.opt);
[Br1,Bt1]=A1.field(SL.Usurf,Rm,thq); Br1=Br1(:).'; Bt1=Bt1(:).';
ff=@(y)abs((2/numel(y))*sum(y.*exp(-1j*M.p*thq)));
q('Bg1 a vide [T]',ff(Br0),0.942);
q('Bt rms a vide [T]',sqrt(mean(Bt0.^2)),0.128);
q('Bg1 en charge [T]',ff(Br1),0.920);
q('Bt rms en charge [T]',sqrt(mean(Bt1.^2)),0.131);

%% ---- 6. caracteristique ----------------------------------------------
fprintf('\n  ---- 6. CARACTERISTIQUE (30 glissements) ----\n');
sl=[0.005:0.005:0.12,0.15,0.2,0.3,0.5,0.7,1.0];
T=nan(size(sl)); I=nan(size(sl)); Xp=ctx.Xm0; b6=0;
for k=1:numel(sl)
    rk=mec.equivalent_circuit(ctx,sl(k),Xp); Xp=rk.Xm;
    T(k)=rk.Tem; I(k)=rk.I1;
    c=mec.power_balance(rk,M); b6=max(b6,c.err_global);
end
[Tmx,im]=max(T);
fprintf('  bilan de puissance : %.2e sur %d glissements\n',b6,numel(sl));
fprintf('  decrochage : %.1f N.m a s = %.3f\n',Tmx,sl(im));
fprintf('  %6s %10s %10s\n','s','T [N.m]','I1 [A]');
for k=1:numel(sl)
    fprintf('  %6.3f %10.2f %10.2f\n',sl(k),T(k),I(k));
end
save('Z5_sweep_cav.mat','sl','T','I','r','r0','rb');
fprintf('\n  duree %.0f s\n=== Z5 termine ===\n',toc(t0));
diary off;
