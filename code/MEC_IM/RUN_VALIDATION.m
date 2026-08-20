%% RUN_VALIDATION  -  Comparaison rigoureuse MEC vs ANSYS Maxwell (complet)
%
%  Confronte le modele MEC a l'ENSEMBLE des performances extraites d'ANSYS
%  Maxwell 2D : couple, vitesse, courant, puissances, rendement, facteur de
%  puissance, ET la decomposition des pertes (fer, cuivre stator, barres
%  rotor, anneaux, mecaniques) — a vide, en charge et rotor bloque.
%
%  METHODE. Les grandeurs ANSYS sont des TRANSITOIRES instantanes fortement
%  ondules : on les moyenne sur une fenetre de REGIME ETABLI (derniere
%  seconde), et l'on recalcule rendement et cos(phi) a partir des puissances
%  MOYENNES (les valeurs instantanees d'ANSYS oscillent entre 0,57 et 0,93 et
%  ne sont pas comparables telles quelles). Le glissement du point de
%  comparaison est deduit de la vitesse moyenne ANSYS : le MEC est evalue
%  EXACTEMENT au meme glissement.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
ROOT = 'C:\Users\hp\Desktop\ANSYS résultat 18.5KW';

fprintf('==================================================================\n');
fprintf('     VALIDATION MEC vs ANSYS Maxwell 2D  -  18,5 kW / 690 V\n');
fprintf('==================================================================\n');

M = mec.machine_18_5kW(); ctx = mec.build_context(M);

%% ---------- helpers ----------
rd  = @(f) readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
% moyenne sur la fenetre de regime etabli [t0, tend]
    function v = avgw(A, col, t0)
        m = A(:,1) >= t0;  v = mean(A(m,col),'omitnan');
    end
    function v = rmsw(A, col, t0)
        m = A(:,1) >= t0;  v = sqrt(mean(A(m,col).^2,'omitnan'));
    end

%% ================= EN CHARGE (point nominal) =================
D = fullfile(ROOT,'transitoire','en charge');
t0 = 1.0;                                  % fenetre de regime etabli [s]

spd = rd(fullfile(D,'Speed Plot 1.tab'));
trq = rd(fullfile(D,'Plot 1.tab'));
pw  = rd(fullfile(D,'Output Variables Plot 3.tab'));   % Pem, Pout, Pin [kW]
cl  = rd(fullfile(D,'Loss Plot 1.tab'));               % CoreLoss [W]
sl  = rd(fullfile(D,'Loss Plot 2.tab'));               % SolidLoss, StrandedLossR [kW]
pro = rd(fullfile(D,'Output Variables Plot 4.tab'));   % Pro [W] (pertes meca)
cur = rd(fullfile(D,'Winding Plot 4.tab'));            % courants instantanes
rng = rd(fullfile(D,'End Connection Plot 3.tab'));     % RingSolidLoss [kW]
ec  = rd(fullfile(D,'End Connection Plot 1.tab'));     % I_barre, I_anneau [kA]

n_ans   = avgw(spd,2,t0);
s_ans   = (M.ns*60 - n_ans)/(M.ns*60);
T_ans   = avgw(trq,3,t0);
TL_ans  = avgw(trq,2,t0);
Pem_ans = avgw(pw,2,t0)*1e3;  Pout_ans = avgw(pw,3,t0)*1e3;  Pin_ans = avgw(pw,4,t0)*1e3;
Pfe_ans = avgw(cl,2,t0);
Pbar_ans= avgw(sl,2,t0)*1e3;  Pcus_ans = avgw(sl,3,t0)*1e3;
Pring_ans = avgw(rng,2,t0)*1e3;
Pmec_ans= avgw(pro,2,t0);
I_ans   = rmsw(cur,2,t0);
eta_ans = Pout_ans/Pin_ans;
PF_ans  = Pin_ans/(M.m*M.Uph*I_ans);

% --- MEC evalue AU MEME glissement ---
r = mec.equivalent_circuit(ctx, s_ans, ctx.Xm0);
Pcur_mec = r.Pcu_r;                        % barres + anneaux (reporte)

fprintf('\n########## EN CHARGE (regime etabli, t > %.1f s) ##########\n',t0);
fprintf('Vitesse ANSYS moyenne = %.1f tr/min  ->  glissement s = %.4f\n',n_ans,s_ans);
fprintf('(couple de charge applique = %.1f N.m)\n\n',TL_ans);
fprintf('%-26s %12s %12s %10s\n','Grandeur','MEC','ANSYS','ecart %');
prt('Couple electromag. [N.m]', r.Tem, T_ans);
prt('Courant I1 [A rms]',       r.I1,  I_ans);
prt('Courant de barre [A rms]', r.Ibar,  1e3*rmsw(ec,2,t0));
prt('Courant d''anneau [A rms]',r.Iring, 1e3*rmsw(ec,3,t0));
prt('Puissance entree [W]',     r.Pin, Pin_ans);
prt('Puissance sortie [W]',     r.Pout,Pout_ans);
prt('Rendement [-]',            r.eta, eta_ans);
prt('Facteur de puissance [-]', r.cosphi, PF_ans);
fprintf('  -- decomposition des pertes --\n');
prt('Pertes fer [W]',           r.Pfe,  Pfe_ans);
prt('Pertes Joule stator [W]',  r.Pcu_s,Pcus_ans);
prt('Pertes rotor (barres+ann.) [W]', Pcur_mec, Pbar_ans+Pring_ans);
prt('   dont barres [W]',       r.Pbars, Pbar_ans);
prt('   dont anneaux [W]',      r.Pring, Pring_ans);
prt('Pertes mecaniques [W]',    r.Pfw,  Pmec_ans);
% pertes supplementaires en charge : cote EF = residu non identifie du bilan ;
% cote MEC = terme de pertes supplementaires (mec.stray_losses, loi I2^2)
stray_ans = (Pin_ans-Pout_ans) - (Pfe_ans+Pcus_ans+Pbar_ans+Pring_ans+Pmec_ans);
prt('Pertes suppl. en charge [W]', r.Padd, stray_ans);

% --- bilan de pertes : (Pin - Pout) vs somme des pertes (stray incluses) ---
sumL_ans = Pfe_ans + Pcus_ans + Pbar_ans + Pring_ans + Pmec_ans + stray_ans;
sumL_mec = r.Pfe + r.Pcu_s + r.Pcu_r + r.Pfw + r.Padd;
fprintf('  -- bilan de pertes --\n');
fprintf('%-26s %12.1f %12.1f\n','Pin - Pout [W]', r.Pin-r.Pout, Pin_ans-Pout_ans);
fprintf('%-26s %12.1f %12.1f\n','Somme pertes (stray incl.)', sumL_mec, sumL_ans);
fprintf('%-26s %12.1f %12.1f  <== residu (doit etre ~0)\n','Ecart [W]', ...
        (r.Pin-r.Pout)-sumL_mec, (Pin_ans-Pout_ans)-sumL_ans);
fprintf('%-26s %12.1f %12.1f  (%.1f %% / %.1f %% Pn, IEC 60034-2-1)\n','  dont suppl. en charge', ...
        r.Padd, stray_ans, 100*r.Padd/M.Pn, 100*stray_ans/M.Pn);

%% ================= A VIDE =================
D0 = fullfile(ROOT,'transitoire','a vide');
spd0 = rd(fullfile(D0,'la vitesse en fonction du temps.tab'));
cur0 = rd(fullfile(D0,'Winding Plot 4.tab'));
cl0  = rd(fullfile(D0,'Loss Plot 1.tab'));
sl0  = rd(fullfile(D0,'Loss Plot 2.tab'));
vi0  = rd(fullfile(D0,'Winding Plot 2.tab'));   % tension induite
ec0  = rd(fullfile(D0,'End Connection Plot 1.tab'));
t00  = 1.0;
n0   = avgw(spd0,2,t00);  s0 = max((M.ns*60-n0)/(M.ns*60), 1e-4);
I0_ans = rmsw(cur0,2,t00);
Pfe0_ans = avgw(cl0,2,t00);
E0_ans = rmsw(vi0,2,t00);
r0 = mec.equivalent_circuit(ctx, s0, ctx.Xm0);

fprintf('\n########## A VIDE (regime etabli) ##########\n');
fprintf('Vitesse ANSYS = %.1f tr/min  -> s = %.5f\n\n',n0,s0);
fprintf('%-26s %12s %12s %10s\n','Grandeur','MEC','ANSYS','ecart %');
prt('Courant a vide I0 [A rms]', r0.I1, I0_ans);
prt('Pertes fer [W]',            r0.Pfe, Pfe0_ans);
prt('f.e.m. induite [V rms]',    r0.E1,  E0_ans);
prt('Courant de barre [A rms]',  r0.Ibar,  1e3*rmsw(ec0,2,t00));
prt('Courant d''anneau [A rms]', r0.Iring, 1e3*rmsw(ec0,3,t00));

%% ================= ROTOR BLOQUE =================
Db = fullfile(ROOT,'transitoire','rotor bloqué');
curb = rd(fullfile(Db,'Winding Plot 3.tab'));
ecb  = rd(fullfile(Db,'End Connection Plot 1.tab'));
trqb = rd(fullfile(Db,'Torque Plot 1.tab'));   % [kN.m]
pwb  = rd(fullfile(Db,'Output Variables Plot 2.tab'));
clb  = rd(fullfile(Db,'Loss Plot 1.tab'));
slb  = rd(fullfile(Db,'Loss Plot 2.tab'));
t0b  = 0.5*max(curb(:,1));
Ib_ans = rmsw(curb,2,t0b);
Tb_ans = avgw(trqb,2,t0b)*1e3;
Pinb_ans = avgw(pwb,4,t0b)*1e3;
Pfeb_ans = avgw(clb,2,t0b);
Pbarb_ans= avgw(slb,2,t0b)*1e3;  Pcusb_ans = avgw(slb,3,t0b)*1e3;
rb = mec.equivalent_circuit(ctx, 1.0, ctx.Xm0);

fprintf('\n########## ROTOR BLOQUE (s = 1) ##########\n\n');
fprintf('%-26s %12s %12s %10s\n','Grandeur','MEC','ANSYS','ecart %');
prt('Couple de demarrage [N.m]', rb.Tem, Tb_ans);
prt('Courant I [A rms]',         rb.I1,  Ib_ans);
prt('Courant de barre [A rms]',  rb.Ibar,  1e3*rmsw(ecb,2,t0b));
prt('Courant d''anneau [A rms]', rb.Iring, 1e3*rmsw(ecb,3,t0b));
prt('Puissance entree [W]',      rb.Pin, Pinb_ans);
prt('Pertes fer [W]',            rb.Pfe, Pfeb_ans);
prt('Pertes Joule stator [W]',   rb.Pcu_s, Pcusb_ans);
prt('Pertes rotor [W]',          rb.Pcu_r, Pbarb_ans);

fprintf('\n=== Termine ===\n');

%% ---------- affichage ----------
function prt(lab, mec_v, ans_v)
    if isnan(mec_v)
        fprintf('%-26s %12s %12.2f %10s\n', lab, '-', ans_v, '-');
    else
        e = 100*(mec_v-ans_v)/ans_v;
        fprintf('%-26s %12.3f %12.3f %+9.1f\n', lab, mec_v, ans_v, e);
    end
end
