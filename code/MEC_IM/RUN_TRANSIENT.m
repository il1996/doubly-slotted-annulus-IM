%% RUN_TRANSIENT  -  Démarrage direct transitoire du moteur asynchrone (MEC->dq)
%
%  Simule le démarrage direct sur le réseau (direct-on-line) par le modèle
%  dynamique dq dont les paramètres proviennent du réseau MEC (Lm saturée,
%  Rr'/fuites de la cage, Rs), couplé à l'équation mécanique
%       J dOmega/dt = Tem - TL - B*Omega.
%
%  Validation contre le transitoire EF (ANSYS Maxwell 2D, essai à vide) :
%  vitesse et couple électromagnétique en fonction du temps.
%
%  Exécuter après RUN_MEC_IM (mêmes modules).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));

fprintf('=== Démarrage transitoire (MEC -> dq) - moteur 18,5 kW à vide ===\n\n');

M   = mec.machine_18_5kW();
ctx = mec.build_context(M);

% Simulation du démarrage à vide
tend = 0.4;
out  = mec.dq_startup(ctx, struct('tend',tend,'TL',0));

fprintf('Parametres dq (issus du MEC) :\n');
fprintf('  Rs=%.3f  Rr''=%.3f ohm   Lm=%.1f mH  Lls=%.2f mH  Llr=%.2f mH\n', ...
    out.params.Rs, out.params.Rr, out.params.Lm*1e3, out.params.Lls*1e3, out.params.Llr*1e3);
fprintf('  J=%.3f kg.m2   B=%.3f N.m.s/rad\n\n', out.params.J, out.params.B);

% ---- Référence EF (ANSYS, essai à vide) ----
base = 'C:\Users\hp\Desktop\ANSYS résultat 18.5KW\transitoire\a vide';
spd = load_tab(fullfile(base,'la vitesse en fonction du temps.tab'));
trq = load_tab(fullfile(base,'Torque Plot 1.tab'));   % [t, LoadTorque, Torque]

% temps de montée à 95 % de la vitesse synchrone
n95 = 0.95*1500;
t95_mec = out.t(find(out.n_rpm>=n95,1));
t95_fem = spd(find(spd(:,2)>=n95,1),1);
fprintf('Temps de montee a 95%% (1425 rpm) : MEC=%.3f s   EF=%.3f s\n', t95_mec, t95_fem);
fprintf('Vitesse finale : MEC=%.1f rpm   EF=%.1f rpm\n', out.n_rpm(end), spd(end,2));
fprintf('Couple crete transitoire : MEC=%.0f N.m\n', max(out.Tem));
fprintf('Courant crete (inrush) phase : MEC=%.0f A\n\n', max(abs(out.ia)));

% ---- Figures ----
figure('Name','Demarrage transitoire MEC vs EF','Position',[80 80 1100 700]);

subplot(2,2,1);
plot(out.t, out.n_rpm,'b-','LineWidth',1.4); hold on;
plot(spd(:,1), spd(:,2),'r--','LineWidth',1.1);
grid on; xlim([0 tend]); xlabel('temps [s]'); ylabel('vitesse [tr/min]');
legend('MEC->dq','EF (ANSYS)','Location','southeast'); title('Vitesse - démarrage à vide');

subplot(2,2,2);
plot(out.t, out.Tem,'b-','LineWidth',1.0); hold on;
plot(trq(:,1), trq(:,3),'r--','LineWidth',1.0);
grid on; xlim([0 tend]); xlabel('temps [s]'); ylabel('couple [N.m]');
legend('MEC->dq','EF (ANSYS)','Location','northeast'); title('Couple électromagnétique');

subplot(2,2,3);
plot(out.t, out.ia,'b-','LineWidth',0.8);
grid on; xlim([0 tend]); xlabel('temps [s]'); ylabel('courant phase A [A]');
title('Courant statorique (inrush)');

subplot(2,2,4);
plot(out.n_rpm, out.Tem,'b-','LineWidth',0.8); grid on;
xlabel('vitesse [tr/min]'); ylabel('couple [N.m]');
title('Trajectoire couple-vitesse (démarrage)');

saveas(gcf, fullfile(fileparts(mfilename('fullpath')),'transitoire_MEC.png'));
fprintf('Figure enregistree : transitoire_MEC.png\n');

%% ---- Démarrage EN CHARGE (profil de charge TL(t) relevé d'ANSYS) ----
Dc = 'C:\Users\hp\Desktop\ANSYS résultat 18.5KW\transitoire\en charge';
trqc = load_tab(fullfile(Dc,'Plot 1.tab'));        % [t, -LoadTorque, Torque]
spdc = load_tab(fullfile(Dc,'Speed Plot 1.tab'));
curc = load_tab(fullfile(Dc,'Winding Plot 4.tab'));
oc = mec.dq_startup(ctx, struct('tend',2.0,'TL',[trqc(:,1),trqc(:,2)]));
mfin = @(x,t) mean(x(t>1.5));
fprintf('\n--- EN CHARGE (TL(t) ANSYS, palier %.1f N.m) ---\n', mfin(trqc(:,2),trqc(:,1)));
fprintf('t95 : MEC=%.3f s   EF=%.3f s\n', ...
    oc.t(find(oc.n_rpm>=n95,1)), spdc(find(spdc(:,2)>=n95,1),1));
fprintf('regime etabli : n=%.1f (EF %.1f) tr/min ; T=%.1f (EF %.1f) N.m ; I1=%.2f (EF %.2f) A rms\n', ...
    mfin(oc.n_rpm,oc.t), mfin(spdc(:,2),spdc(:,1)), ...
    mfin(oc.Tem,oc.t),  mfin(trqc(:,3),trqc(:,1)), ...
    sqrt(mean(oc.ia(oc.t>1.5).^2)), sqrt(mean(curc(curc(:,1)>1.5,2).^2)));

figure('Name','Demarrage EN CHARGE MEC vs EF','Position',[100 100 1100 380]);
subplot(1,3,1);
plot(oc.t,oc.n_rpm,'b-','LineWidth',1.3); hold on; plot(spdc(:,1),spdc(:,2),'r--');
grid on; xlim([0 2]); xlabel('temps [s]'); ylabel('vitesse [tr/min]');
legend('MEC->dq','EF','Location','southeast'); title('Vitesse — en charge');
subplot(1,3,2);
plot(oc.t,oc.Tem,'b-'); hold on; plot(trqc(:,1),trqc(:,3),'r--');
grid on; xlim([0 2]); xlabel('temps [s]'); ylabel('couple [N.m]');
title('Couple — en charge (TL(t) ANSYS)');
subplot(1,3,3);
plot(oc.t,oc.ia,'b-'); hold on; plot(curc(:,1),curc(:,2),'r--');
grid on; xlim([0 2]); xlabel('temps [s]'); ylabel('i_A [A]');
title('Courant — en charge'); legend('MEC','EF');
saveas(gcf, fullfile(fileparts(mfilename('fullpath')),'transitoire_charge_MEC.png'));
fprintf('Figure enregistree : transitoire_charge_MEC.png\n');
fprintf('\n=== Termine ===\n');

% ================== utilitaires ==================
function A = load_tab(path)
%LOAD_TAB  Charge un fichier .tab ANSYS (en-tête entre guillemets, tab).
    if ~isfile(path)
        warning('Fichier EF absent : %s', path); A = [NaN NaN]; return;
    end
    A = readmatrix(path,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
end
