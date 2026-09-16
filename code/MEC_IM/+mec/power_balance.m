function chk = power_balance(res, M)
%POWER_BALANCE  Vérification de cohérence énergétique (test B6).
%
%   chk = mec.power_balance(res,M) vérifie, à chaque point de
%   fonctionnement, le bilan de puissance :
%       Pin = Pcu_stator + Pfe + Pag           (entrée -> entrefer)
%       Pag = Pcu_rotor + Pmec_brute           (partage au glissement)
%       Pin = Pout + (Pcu_s + Pfe + Pcu_r + Pfw)
%   et l'identité couple/puissance  Pag = Tem * Omega_synchrone.
%
%   C'est le test le plus rentable du programme (défaut P5/B6 de Silva,
%   qui ne présente JAMAIS de bilan de puissance) : il détecte toute
%   incohérence de bookkeeping ou de référence d'impédance.
%
%   Sorties : chk.err_global, chk.err_torque (erreurs relatives).

ws = M.w/M.p;
Padd = 0; if isfield(res,'Padd'), Padd = res.Padd; end   % pertes suppl. en charge
sumLoss = res.Pcu_s + res.Pfe + res.Pcu_r + res.Pfw + Padd;
lhs = res.Pin;
rhs = res.Pout + sumLoss;
chk.err_global = abs(lhs - rhs)/max(abs(lhs),1);

% Identité couple/puissance sur le FONDAMENTAL (les harmoniques ont leur
% propre vitesse synchrone Omega_s/nu : Pag_nu = T_nu * Omega_s,nu).
if isfield(res,'Tem1'), T1 = res.Tem1; else, T1 = res.Tem; end
chk.err_torque = abs(res.Pag - T1*ws)/max(abs(res.Pag),1);

chk.Pin = res.Pin; chk.Pout = res.Pout; chk.sumLoss = sumLoss;
chk.ok = (chk.err_global < 0.02) && (chk.err_torque < 1e-6);
end
