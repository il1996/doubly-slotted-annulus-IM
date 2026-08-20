%% RUN_COMPARE_ANALYTIC  -  Comparaison A TROIS VOIES : ANALYTIQUE / MEC / EF
%
%  Confronte, sur la MEME machine (18,5 kW, 690 V, 4 poles, 48/44, g=0,243 mm,
%  M800-50A), les trois approches :
%    * ANALYTIQUE : methode de conception de Pyrhonen-Jokinen-Hrabovcova
%      (schema equivalent a parametres analytiques ; analytical_18_5kW,
%      version fonction de main_motor_IM.m — dossier conception) ;
%    * MEC        : reseau de reluctances non lineaire (mec.equivalent_circuit) ;
%    * EF (ANSYS) : reference Maxwell 2D (mec.ansys_ref + essais transitoires).
%
%  Sorties : (1) tableau du point nominal (parametres + performances + pertes) ;
%  (2) courbes couple-glissement et courant-glissement a trois voies ;
%  (3) ecarts chiffres vs EF.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
addpath('C:\Users\hp\Desktop\Matlab program\conception\New Folder');   % analytique
OUT=fileparts(mfilename('fullpath'));

fprintf('==================================================================\n');
fprintf('   COMPARAISON A TROIS VOIES : ANALYTIQUE / MEC / EF (18,5 kW)\n');
fprintf('==================================================================\n\n');

% --- les trois modeles ---
A   = analytical_18_5kW();                 % methode analytique (Pyrhonen)
M   = mec.machine_18_5kW(); ctx=mec.build_context(M); ref=mec.ansys_ref();

% point de comparaison = glissement de l'essai EF en charge
s_cmp = 0.0188;
rM = mec.equivalent_circuit(ctx,s_cmp,ctx.Xm0);
Cg = mec.cage(M,ctx.G,ctx.W,ctx.Lk,s_cmp);

% grandeurs EF de reference (essai en charge, moyennes de regime etabli)
EF.I=19.73; EF.T=121.63; EF.Xm=46; EF.Rs=0.445; EF.Rr=0.44;
EF.Pfe=232.6; EF.Pcus=520.2; EF.Pcur=488.5; EF.eta=0.916; EF.cosphi=0.864;

%% ======== 1. PARAMETRES DU SCHEMA EQUIVALENT ========
fprintf('----- 1. PARAMETRES DU SCHEMA EQUIVALENT -----\n');
fprintf('%-22s %12s %12s %12s\n','grandeur','ANALYTIQUE','MEC','EF (ANSYS)');
row('Rs [ohm]',       A.Rs,      ctx.Rs,   EF.Rs);
row('Rr'' [ohm]',     A.Rr,      rM.Rr,    EF.Rr);
row('Xm [ohm]',       A.Xm,      rM.Xm,    EF.Xm);
row('X_sigma_s [ohm]',A.Xssigma, ctx.Lk.Xs_leak, NaN);
row('X_sigma_r [ohm]',A.Xrsigma, rM.Xr,    NaN);
row('Rfe [ohm]',      A.Rfe,     rM.Rfe,   1740);

%% ======== 2. POINT NOMINAL (au glissement EF s=0,0188) ========
fprintf('\n----- 2. POINT NOMINAL (s = %.4f, glissement EF en charge) -----\n',s_cmp);
TA=interp1(A.ss,A.T,s_cmp);  IA=interp1(A.ss,A.i1,s_cmp);
etaA=interp1(A.ss,A.eff,s_cmp); cosA=interp1(A.ss,A.cosfi,s_cmp);
fprintf('%-22s %12s %12s %12s | %10s %10s\n','grandeur','ANALYTIQUE','MEC','EF','e_ana%','e_mec%');
row3('Couple T [N.m]', TA, rM.Tem, EF.T);
row3('Courant I1 [A]', IA, rM.I1,  EF.I);
row3('Rendement eta',  etaA, rM.eta, EF.eta);
row3('cos(phi)',       cosA, rM.cosphi, EF.cosphi);
fprintf('  -- pertes [W] --\n');
row3('Pertes fer',     A.Pfe, rM.Pfe, EF.Pfe);
row3('Pertes Cu stator', A.m*IA^2*A.Rs, rM.Pcu_s, EF.Pcus);
row3('Pertes Cu rotor',  A.m*interp1(A.ss,A.i2,s_cmp)^2*A.Rr, rM.Pcu_r, EF.Pcur);

%% ======== 3. POINTS NOMINAUX PROPRES (chaque methode a sa charge) ========
fprintf('\n----- 3. POINT NOMINAL PROPRE (chaque methode a P_utile = 18,5 kW) -----\n');
fprintf('%-22s %12s %12s %12s\n','grandeur','ANALYTIQUE','MEC','EF');
sM=ref.anchor.s_rated; rMr=mec.equivalent_circuit(ctx,sM,ctx.Xm0);
sbrk=0.05:0.01:0.2; Tbrk=arrayfun(@(s)mec.equivalent_circuit(ctx,s,ctx.Xm0).Tem,sbrk);
[TbM,ibM]=max(Tbrk);
fprintf('%-22s %12.4f %12.4f %12.4f\n','glissement nominal',A.s_nom,sM,0.0188);
fprintf('%-22s %12.1f %12.1f %12.1f\n','couple nominal [N.m]',A.Tn,rMr.Tem,EF.T);
fprintf('%-22s %12.1f %12.1f %12.1f\n','couple decrochage [N.m]',A.Tb,TbM,ref.anchor.T_break);
fprintf('%-22s %12.4f %12.4f %12.4f\n','glissement decrochage',A.sb,sbrk(ibM),ref.anchor.s_break);

%% ======== 4. COURBES A TROIS VOIES ========
sV=[0.005:0.005:0.05,0.06:0.02:0.2,0.25:0.05:1.0];
TM=zeros(size(sV)); IM=TM;
for k=1:numel(sV)
    r=mec.equivalent_circuit(ctx,sV(k),ctx.Xm0); TM(k)=r.Tem; IM(k)=r.I1;
end

figure('Name','Comparaison 3 voies','Position',[60 60 1150 760]);
subplot(2,2,1);
plot(A.ss,A.T,'g-','LineWidth',1.6); hold on;
plot(sV,TM,'b-o','LineWidth',1.3,'MarkerSize',4);
plot(ref.s,ref.T,'r--s','LineWidth',1.2); grid on;
xlabel('glissement s'); ylabel('Couple [N.m]');
legend('Analytique','MEC','EF (ANSYS)','Location','best');
title('(a) Couple–glissement');

subplot(2,2,2);
plot(A.ss,A.i1,'g-','LineWidth',1.6); hold on;
plot(sV,IM,'b-o','LineWidth',1.3,'MarkerSize',4);
plot(ref.s2,ref.I,'r--s','LineWidth',1.2); grid on;
xlabel('glissement s'); ylabel('Courant I_1 [A]');
legend('Analytique','MEC','EF (ANSYS)','Location','best');
title('(b) Courant–glissement');

subplot(2,2,3);
plot(A.ss,A.eff,'g-','LineWidth',1.6); hold on;
etaM=zeros(size(sV));
for k=1:numel(sV), r=mec.equivalent_circuit(ctx,sV(k),ctx.Xm0); etaM(k)=r.eta; end
plot(sV,etaM,'b-o','LineWidth',1.3,'MarkerSize',4);
plot(s_cmp,EF.eta,'rs','MarkerSize',10,'LineWidth',2); grid on; ylim([0 1]);
xlabel('glissement s'); ylabel('rendement'); xlim([0 0.2]);
legend('Analytique','MEC','EF (nominal)','Location','best');
title('(c) Rendement (zone utile)');

subplot(2,2,4);
plot(A.ss,A.cosfi,'g-','LineWidth',1.6); hold on;
cosM=zeros(size(sV));
for k=1:numel(sV), r=mec.equivalent_circuit(ctx,sV(k),ctx.Xm0); cosM(k)=r.cosphi; end
plot(sV,cosM,'b-o','LineWidth',1.3,'MarkerSize',4);
plot(s_cmp,EF.cosphi,'rs','MarkerSize',10,'LineWidth',2); grid on; ylim([0 1]);
xlabel('glissement s'); ylabel('cos(\phi)'); xlim([0 0.2]);
legend('Analytique','MEC','EF (nominal)','Location','best');
title('(d) Facteur de puissance (zone utile)');
saveas(gcf,fullfile(OUT,'compare_3voies.png'));

fprintf('\nFigure enregistree : compare_3voies.png\n');
fprintf('\n=== Termine ===\n');

%% ================= fonctions locales =================
function row(lab,va,vm,ve)
    if isnan(ve)
        fprintf('%-22s %12.4f %12.4f %12s\n',lab,va,vm,'-');
    else
        fprintf('%-22s %12.4f %12.4f %12.4f\n',lab,va,vm,ve);
    end
end
function row3(lab,va,vm,ve)
    fprintf('%-22s %12.3f %12.3f %12.3f | %+9.1f %+9.1f\n',...
        lab,va,vm,ve,100*(va-ve)/ve,100*(vm-ve)/ve);
end
