%% RUN_STEPPING  -  Pertes supplementaires par un modele PAS-A-PAS (physique)
%
%  OBJET. Reproduire EXPLICITEMENT, a partir de la physique, les pertes
%  supplementaires (~380 W) actuellement introduites comme une CONSTANTE dans
%  ANSYS. Mecanisme dominant vise : les PERTES PAR PULSATION DE DENTURE — les
%  encoches rotoriques modulent le flux des dents statoriques a Nr*n (~1078 Hz
%  au nominal) et les encoches statoriques modulent celui des dents rotoriques
%  a Ns*n. Ces harmoniques sont invisibles a un modele fondamental.
%
%  METHODE (pas-a-pas en POSITION = en temps en regime etabli) :
%    1. le champ 2D charge est resolu a N positions rotoriques successives,
%       avec l'avance SYNCHRONE des courants statoriques et l'avance a la
%       frequence de glissement des courants de barres ;
%    2. on enregistre la forme d'onde B(t) de CHAQUE element de fer ;
%    3. les pertes sont calculees dans le DOMAINE TEMPOREL :
%           P_Foucault = (sigma*d^2/12) * <(dB/dt)^2> * Volume
%       forme EXACTE pour une tole mince, valable pour une forme d'onde
%       QUELCONQUE. On n'utilise PAS de FFT : le glissement rend la denture
%       non commensurable avec la periode electrique (44*n/50 = 21,6 : non
%       entier), ce qui provoquerait une fuite spectrale.
%    4. on compare a la meme perte calculee sur le FONDAMENTAL seul :
%       la difference EST la perte supplementaire de pulsation de denture.
%
%  sigma et d sont les grandeurs PHYSIQUES de la tole (cf. C1) : aucun
%  coefficient ajuste n'intervient dans ce calcul.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== Pertes supplementaires par modele pas-a-pas ===\n\n');

M = mec.machine_18_5kW(); ctx = mec.build_context(M);
G = ctx.G; W = ctx.W; BH = ctx.BH; p = M.p; Nr = M.Nr;
sig = M.fe.sigma; d = M.fe.d; rho = M.fe.rho;

% ---- point nominal ----
s = 0.0188;                                  % glissement ANSYS mesure
r = mec.equivalent_circuit(ctx, s, ctx.Xm0);
wm = (M.w/p)*(1-s);                          % vitesse mecanique [rad/s]
fprintf('Point nominal : s=%.4f, I1=%.2f A, I2=%.2f A, n=%.0f tr/min\n', ...
        s, r.I1, r.I2, wm*60/(2*pi));

psi1 = angle(r.I1c); psi2 = angle(r.I2c);
Ibar = r.I2*(2*M.m*W.kw1*W.Nph)/Nr;
thbar = 2*pi*(0:Nr-1)'/Nr;

% ---- balayage : une periode ELECTRIQUE = pi/p*... en position mecanique ----
% t = theta/wm ; une periode electrique T=1/f correspond a dtheta = wm*T
nc = 4; nr = 3;
Npos = 360;                                   % ~16 points par passage d'encoche
Tel  = 1/M.f;
th   = linspace(0, wm*Tel, Npos+1); th(end) = [];   % positions mecaniques
tt   = th/wm;                                        % temps correspondant
fprintf('Balayage : %d positions sur %.1f deg mec (= 1 periode electrique)\n', ...
        Npos, (wm*Tel)*180/pi);
fprintf('Passages d''encoche rotorique sur la fenetre : %.1f\n\n', Nr*(wm*Tel)/(2*pi));

% ---- pas-a-pas ----
mesh0 = mec.mesh_refined(M,G,nc,sqrt(2)*r.I1*[1;-0.5;-0.5],0,nr,zeros(Nr,1));
nB = numel(mesh0.a); isFe = mesh0.iron(:);
Ball = zeros(Npos, nnz(isFe));
tic;
for k = 1:Npos
    ph  = M.w*tt(k);                                   % avance elec. synchrone
    phs = s*M.w*tt(k);                                 % avance a la freq. de glissement
    i3k = sqrt(2)*r.I1*[cos(psi1+ph); cos(psi1+ph-2*pi/3); cos(psi1+ph+2*pi/3)];
    ibk = sqrt(2)*Ibar*cos(p*thbar - psi2 - phs);
    me  = mec.mesh_refined(M,G,nc,i3k,th(k),nr,ibk);
    Se  = mec.solve_mesh(me,BH,M.opt);
    Ball(k,:) = Se.B(isFe).';
end
fprintf('Pas-a-pas : %d resolutions en %.1f s\n\n', Npos, toc);

% ---- geometrie des elements de fer ----
lFe = mesh0.l(isFe); AFe = mesh0.A(isFe);
VolFe = lFe.*AFe;                                   % volume de fer [m3]
mFe   = rho*VolFe;
kind  = mesh0.iron;                                  %#ok<NASGU>
% masque stator / rotor : les noeuds stator sont numerotes en premier
isStator = mesh0.a(isFe) <= mesh0.Ms*mesh0.Ls;

% ---- dB/dt (domaine temporel, sans FFT) ----
dt   = tt(2)-tt(1);
dBdt = zeros(size(Ball));
for j = 1:size(Ball,2)
    b = Ball(:,j);
    dBdt(:,j) = (circshift(b,-1) - circshift(b,1))/(2*dt);   % periodique
end

% ---- pertes Foucault : forme EXACTE en temps ----
kEddy = sig*d^2/12;
P_ed_tot = kEddy * mean(dBdt.^2,1).' .* VolFe;      % [W] par element

% ---- reference FONDAMENTALE : meme formule sur le 1er harmonique seul ----
% extraction du fondamental de chaque B(t) (projection sur exp(j*w*t))
c1 = (2/Npos)*sum(Ball.*exp(-1i*M.w*tt.'),1);        % amplitude complexe
B1 = abs(c1).';                                      % amplitude du fondamental
P_ed_fond = kEddy * (M.w^2*B1.^2/2) .* VolFe;        % <(dB/dt)^2>=w^2*B1^2/2

% ---- pertes par EXCES : meme traitement temporel <|dB/dt|^1.5> ----
%  Bertotti : p_exc = ke*(f*B)^1.5 [W/kg] pour une sinusoide. En temporel,
%  p_exc = ke_v*<|dB/dt|^1.5> avec ke_v identifie sur le cas sinusoidal :
%     ke*(f*B)^1.5 = ke_v*(w*B)^1.5*<|cos|^1.5>  =>  ke_v = ke/((2*pi)^1.5*0.8244)
ke_v = M.iron.ke/((2*pi)^1.5*0.8244);
P_ex_tot  = ke_v * mean(abs(dBdt).^1.5,1).' .* mFe;
P_ex_fond = ke_v * ((M.w*B1).^1.5*0.8244) .* mFe;

P_tot   = sum(P_ed_tot);
P_fond  = sum(P_ed_fond);
P_extra = P_tot - P_fond;
Pex_tot = sum(P_ex_tot); Pex_fond = sum(P_ex_fond);

fprintf('--- Pertes par courants de Foucault (sigma, d PHYSIQUES) ---\n');
fprintf('%-34s %10s %10s %10s\n','','TOTAL','fondam.','SUPPL.');
fprintf('%-34s %10.1f %10.1f %10.1f\n','Stator [W]', ...
        sum(P_ed_tot(isStator)), sum(P_ed_fond(isStator)), ...
        sum(P_ed_tot(isStator))-sum(P_ed_fond(isStator)));
fprintf('%-34s %10.1f %10.1f %10.1f\n','Rotor [W]', ...
        sum(P_ed_tot(~isStator)), sum(P_ed_fond(~isStator)), ...
        sum(P_ed_tot(~isStator))-sum(P_ed_fond(~isStator)));
fprintf('%-34s %10.1f %10.1f %10.1f\n','TOTAL [W]', P_tot, P_fond, P_extra);
fprintf('  facteur de majoration harmonique (Foucault) = %.2f\n', P_tot/P_fond);

fprintf('\n--- Pertes par EXCES (coeff. ke de Bertotti, traitement temporel) ---\n');
fprintf('%-34s %10.1f %10.1f %10.1f\n','TOTAL [W]', Pex_tot, Pex_fond, Pex_tot-Pex_fond);
fprintf('  facteur de majoration harmonique (exces)    = %.2f\n', Pex_tot/Pex_fond);

fprintf('\n=================== BILAN ===================\n');
fprintf('Pertes SUPPLEMENTAIRES produites par la physique (pas-a-pas) :\n');
fprintf('   Foucault (pulsation de denture) = %6.1f W\n', P_extra);
fprintf('   Exces    (pulsation de denture) = %6.1f W\n', Pex_tot-Pex_fond);
fprintf('   ---------------------------------------------\n');
fprintf('   TOTAL produit par le modele     = %6.1f W\n', P_extra+(Pex_tot-Pex_fond));
fprintf('   Constante imposee dans ANSYS    = %6.1f W\n', 381.5);
fprintf('   => le modele n''en explique que    %5.1f %%\n', ...
        100*(P_extra+(Pex_tot-Pex_fond))/381.5);

% ---- figure : formes d'onde ----
[~,jw] = max(std(Ball,0,1));                        % element le plus module
figure('Name','Pas-a-pas : formes d''onde','Position',[100 100 950 380]);
subplot(1,2,1);
plot(tt*1e3, Ball(:,jw),'b-','LineWidth',1.2); grid on;
xlabel('temps [ms]'); ylabel('B [T]');
title('B(t) de l''element le plus module (pulsation de denture)');
subplot(1,2,2);
plot(tt*1e3, dBdt(:,jw),'r-','LineWidth',1.0); grid on;
xlabel('temps [ms]'); ylabel('dB/dt [T/s]');
title('dB/dt  (c''est son carre moyen qui fait la perte)');
saveas(gcf, fullfile(fileparts(mfilename('fullpath')),'pas_a_pas.png'));
fprintf('\nFigure enregistree : pas_a_pas.png\n');
