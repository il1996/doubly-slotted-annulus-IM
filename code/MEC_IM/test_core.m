% Test d'intégration du cœur MEC : magnétisation à vide.
clear; clc;
addpath(fileparts(mfilename('fullpath')));   % rend +mec visible

M  = mec.machine_18_5kW();
G  = mec.geometry(M);
W  = mec.winding(M);
mat= mec.mat_M800_50A();
BH = mec.bh(mat);
net= mec.build_network(M,G);

fprintf('--- Geometrie ---\n');
fprintf('taus=%.3f mm  taur=%.3f mm\n',G.taus*1e3,G.taur*1e3);
fprintf('bts =%.3f mm  btr =%.3f mm\n',G.bts*1e3,G.btr*1e3);
fprintf('hys =%.3f mm  hyr =%.3f mm\n',G.hys*1e3,G.hyr*1e3);
fprintf('A_ts=%.2f mm2 A_ys=%.2f mm2\n',G.A_ts*1e6,G.A_ys*1e6);
fprintf('kw1=%.4f  Nph=%.0f  q=%d\n',W.kw1,W.Nph,W.q);

% Courant magnetisant d'essai (instant phase A maximale)
I0 = 8.32;                     % A rms
i3 = sqrt(2)*I0*[1; -0.5; -0.5];
Fs = W.slotMMF(i3);

AG = mec.airgap_permeance(M,G,0);
fprintf('--- Entrefer ---\n');
fprintf('kC=%.4f  g_eff=%.3f mm  P0=%.3e H\n',AG.kC,AG.g_eff*1e3,AG.P0);
fprintf('Psum/P0 (min..max)=%.4f..%.4f\n',min(AG.Psum_stator)/AG.P0,max(AG.Psum_stator)/AG.P0);

S = mec.solve_network(net,G,BH,AG,Fs,[],M.opt);
fprintf('--- Solveur ---\n');
fprintf('Newton: %d iterations, residu=%.2e, converged=%d\n',S.iter,S.res,S.converged);

% Fondamental spatial de l'induction MOYENNE d'entrefer (p paires de poles)
Ns = M.Ns; p = M.p;
th = 2*pi*(0:Ns-1)/Ns;
Bg = S.Bgap_avg_i(:).';
c  = (2/Ns)*sum(Bg.*exp(-1j*p*th));
Bg1 = abs(c);
tp  = pi*G.Ds/(2*p);
Phi1 = (2/pi)*Bg1*tp*M.L;
lam1 = W.kw1*W.Nph*Phi1;
E1   = M.w*lam1/sqrt(2);
Xm   = E1/I0;
fprintf('--- Grandeurs fondamentales ---\n');
fprintf('Bg1 = %.4f T   (cible ~1.01 T)\n',Bg1);
fprintf('Phi1= %.4f mWb\n',Phi1*1e3);
fprintf('E1  = %.2f V\n',E1);
fprintf('Xm  = %.2f ohm (cible ~46 ohm)   Lm=%.1f mH\n',Xm,Xm/M.w*1e3);

% Cartes de B (max par region)
kinds = {'ts','ys','tr','yr'};
for kk=1:numel(kinds)
    m = strcmp(S.kind,kinds{kk});
    fprintf('B_%s : moy=%.3f  max=%.3f T\n',kinds{kk},mean(abs(S.Biron(m))),max(abs(S.Biron(m))));
end

% --- Sonde LINEAIRE (fer quasi infini) pour valider le modele d'entrefer ---
fprintf('\n--- Sonde lineaire (I0=0.5 A) ---\n');
I0l = 0.5; i3l = sqrt(2)*I0l*[1;-0.5;-0.5];
Sl = mec.solve_network(net,G,BH,AG,W.slotMMF(i3l),[],M.opt);
Bgl = Sl.Bgap_i(:).';
Bg1l = abs((2/Ns)*sum(Bgl.*exp(-1j*p*th)));
% analytique : F1 crete fondamental, Bg1_lin = mu0*F1/g_eff
mu0=4*pi*1e-7;
F1 = (3/2)*(4/pi)*(W.kw1*W.Nph/(2*p))*sqrt(2)*I0l;
Bg1_analytic = mu0*F1/AG.g_eff;
fprintf('MEC   Bg1(lin)/I0 = %.4f T/A (bec de dent)\n',Bg1l/I0l);
fprintf('Analyt Bg1(lin)/I0 = %.4f T/A  (F1=%.1f A)\n',Bg1_analytic/I0l,F1);
fprintf('rapport MEC/analytique = %.4f\n',Bg1l/Bg1_analytic);

% Fondamental de la FMM de dents et B moyen d'entrefer (flux / (taus*L))
Fs_l = W.slotMMF(i3l);
Fs_fund = abs((2/Ns)*sum(Fs_l(:).'.*exp(-1j*p*th)));
fprintf('Fs_fund(MEC)=%.2f A   F1(analyt)=%.2f A   ratio=%.3f\n',Fs_fund,F1,Fs_fund/F1);
Bavg = Sl.Phi_gap_stator(:).'/(G.taus*M.L);
Bavg1 = abs((2/Ns)*sum(Bavg.*exp(-1j*p*th)));
fprintf('B_moy_entrefer fond/I0 = %.4f T/A  (vs analyt %.4f)\n',Bavg1/I0l,Bg1_analytic/I0l);
Ak = W.C*i3l;  Ak_fund = abs((2/Ns)*sum(Ak(:).'.*exp(-1j*p*th)));
fprintf('A_k fondamental = %.3f A\n',Ak_fund);
