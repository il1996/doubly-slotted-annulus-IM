%% RUN_MEC_IM  -  Programme principal : etude MEC d'un moteur asynchrone
%
%  Etude electromagnetique complete d'un moteur asynchrone triphase a cage
%  par la methode des circuits magnetiques equivalents (MEC), avec :
%    * reseau de reluctances dent-par-dent, non lineaire, Newton exact ;
%    * saturation (courbe B-H M800-50A) ; entrefer physique (Carter+frange) ;
%    * cage barre + anneaux, effet de peau 1D EXACT (barre trapezoidale) ;
%    * fuites (encoche, bec, tetes de bobines) ;
%    * branches HARMONIQUES explicites (couples parasites asynchrones) ;
%    * pertes fer (Bertotti par region), cuivre, mecaniques ;
%    * couplage MEC <-> circuit equivalent (Xm et Rfe satures) ;
%    * bilan de puissance systematique (B6) ;
%    * validation chiffree contre EF (ANSYS Maxwell 2D).
%
%  GEOMETRIE : cotes reelles fournies (profils tc6 stator / tcr rotor).
%
%  Autres scripts : RUN_TRANSIENT (demarrage), RUN_MESH (convergence 2D),
%  RUN_RIPPLE (ondulation), RUN_FEM_IDENT (sigma0), RUN_C1 (Foucault),
%  RUN_HARMONICS / RUN_BARSKIN / RUN_LEAKAGE / RUN_SLOTLEAK (diagnostics).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));

fprintf('==============================================================\n');
fprintf('   MEC - Moteur asynchrone triphase a cage (18,5 kW, 690 V)\n');
fprintf('==============================================================\n\n');

%% 1. Machine, contexte, reference
M   = mec.machine_18_5kW();
ctx = mec.build_context(M);
ref = mec.ansys_ref();
G   = ctx.G;

fprintf('Geometrie (cotes reelles) :\n');
fprintf('  entrefer g = %.3f mm   kC = %.3f   g_eff = %.3f mm\n', ...
        M.g*1e3, ctx.AG.kC, ctx.AG.g_eff*1e3);
fprintf('  hs = %.3f mm   hr = %.3f mm\n', G.hs*1e3, G.hr*1e3);
fprintf('  bts = %.3f mm (dent parallele, err %.1e)   btr = %.3f mm\n', ...
        G.bts*1e3, G.bts_parallel_err, G.btr*1e3);
fprintf('  hys = %.3f mm   hyr = %.3f mm   Abar = %.2f mm2\n', ...
        G.hys*1e3, G.hyr*1e3, G.Abar*1e6);
fprintf('Bobinage : kw1 = %.4f   Nph = %.0f   q = %d\n', ctx.W.kw1, ctx.W.Nph, ctx.W.q);
fprintf('Rs = %.4f ohm    Xm0 (non sature) = %.2f ohm\n\n', ctx.Rs, ctx.Xm0);

%% 2. Balayage en glissement
s_list = [0.005:0.005:0.12, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0];
N = numel(s_list); Res = cell(N,1); Xm_prev = ctx.Xm0;
tic;
for k = 1:N
    Res{k} = mec.equivalent_circuit(ctx, s_list(k), Xm_prev);
    Xm_prev = Res{k}.Xm;
end
fprintf('Balayage de %d glissements en %.1f s.\n', N, toc);

gv = @(f) cellfun(@(r) r.(f), Res);
T = gv('Tem'); T1 = gv('Tem1'); Th = gv('Tem_h');
I1 = gv('I1'); PF = gv('cosphi'); eta = gv('eta');
Pfe = gv('Pfe'); Pcus = gv('Pcu_s'); Pcur = gv('Pcu_r'); Pfw = gv('Pfw'); Padd = gv('Padd');
Xm = gv('Xm'); nrpm = gv('n_rpm');

%% 3. Bilan de puissance (B6)
maxerr = 0;
for k = 1:N
    c = mec.power_balance(Res{k}, M);
    maxerr = max(maxerr, c.err_global);
end
fprintf('Bilan de puissance (B6) : erreur max = %.2e  -> %s\n\n', ...
        maxerr, ternary(maxerr<0.02,'OK','ECHEC'));

%% 4. Point nominal
sr = ref.anchor.s_rated;
rr = mec.equivalent_circuit(ctx, sr, ctx.Xm0);
fprintf('----- POINT NOMINAL (s = %.3f) -----\n', sr);
fprintf('  Vitesse        = %.0f tr/min\n', rr.n_rpm);
fprintf('  Courant I1     = %.2f A     (EF %.2f)\n', rr.I1, interp1(ref.s2,ref.I,sr));
fprintf('  Couple Tem     = %.1f N.m   (EF %.1f)\n', rr.Tem, interp1(ref.s,ref.T,sr));
fprintf('    dont fondamental %.1f + parasites %.1f N.m\n', rr.Tem1, rr.Tem_h);
fprintf('  Facteur puiss. = %.3f     Rendement = %.3f\n', rr.cosphi, rr.eta);
fprintf('  Xm (sature)    = %.2f ohm  (EF ~%.0f)\n', rr.Xm, ref.anchor.Xm);
fprintf('  Bg1            = %.3f T    (EF ~%.2f)\n', rr.Bg1, ref.anchor.Bg1);
fprintf('  P_in=%.0f W  P_out=%.0f W\n', rr.Pin, rr.Pout);
fprintf('  Pertes: Cu_s=%.0f  Fe=%.0f  Cu_r=%.0f  Mec=%.0f  Suppl=%.0f W\n\n', ...
        rr.Pcu_s, rr.Pfe, rr.Pcu_r, rr.Pfw, rr.Padd);

%% 5. Validation chiffree vs EF
fprintf('----- VALIDATION vs EF (ANSYS) -----\n');
fprintf('%6s | %8s %8s %7s | %8s %8s %7s\n','s','T_MEC','T_EF','err%','I_MEC','I_EF','err%');
for s = [0.005 0.01 0.02 0.05 0.10 0.20 0.50 1.0]
    r = mec.equivalent_circuit(ctx, s, ctx.Xm0);
    Tf = interp1(ref.s, ref.T, s); If = interp1(ref.s2, ref.I, s);
    fprintf('%6.3f | %8.1f %8.1f %6.1f | %8.1f %8.1f %6.1f\n', ...
        s, r.Tem, Tf, 100*(r.Tem-Tf)/Tf, r.I1, If, 100*(r.I1-If)/If);
end
[Tmax, imax] = max(T);
fprintf('\nCouple de decrochage MEC = %.1f N.m a s=%.3f (EF %.1f a s=%.3f)\n', ...
        Tmax, s_list(imax), ref.anchor.T_break, ref.anchor.s_break);

%% 6. Figures
figure('Name','MEC vs EF','Position',[80 80 1100 700]);
subplot(2,3,1); plot(s_list,T,'b-o','LineWidth',1.4); hold on;
plot(ref.s,ref.T,'r--s','LineWidth',1.1); grid on;
xlabel('glissement s'); ylabel('Couple [N.m]'); legend('MEC','EF','Location','best');
title('Couple-glissement');

subplot(2,3,2); plot(s_list,I1,'b-o','LineWidth',1.4); hold on;
plot(ref.s2,ref.I,'r--s','LineWidth',1.1); grid on;
xlabel('glissement s'); ylabel('I_1 [A]'); legend('MEC','EF','Location','best');
title('Courant-glissement');

subplot(2,3,3); plot(s_list,T1,'b-','LineWidth',1.4); hold on;
plot(s_list,Th,'m-','LineWidth',1.4); plot(s_list,T,'k--','LineWidth',1.2);
grid on; xlabel('glissement s'); ylabel('Couple [N.m]');
legend('fondamental','parasites','net','Location','best'); title('Couples parasites');

subplot(2,3,4); plot(nrpm,T,'b-o','LineWidth',1.4); grid on;
xlabel('vitesse [tr/min]'); ylabel('Couple [N.m]'); title('Caracteristique mecanique');

subplot(2,3,5); plot(s_list,eta,'b-o','LineWidth',1.4); hold on;
plot(s_list,PF,'m-^','LineWidth',1.2); grid on; ylim([0 1]);
xlabel('glissement s'); legend('rendement','cos\phi','Location','best');
title('Rendement et cos\phi');

subplot(2,3,6); area(s_list,[Pcus, Pfe, Pcur, Pfw, Padd]); grid on;
xlabel('glissement s'); ylabel('Pertes [W]');
legend('Cu stator','Fer','Cu rotor','Mecaniques','Suppl.','Location','best');
title('Repartition des pertes');
saveas(gcf, fullfile(fileparts(mfilename('fullpath')),'resultats_MEC.png'));
fprintf('\nFigure enregistree : resultats_MEC.png\n');
fprintf('\n=== Termine ===\n');

function o = ternary(c,a,b), if c, o=a; else, o=b; end, end
