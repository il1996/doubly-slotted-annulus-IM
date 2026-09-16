function Rs = stator_resistance(M, G, W, Lk)
%STATOR_RESISTANCE  Résistance de phase du bobinage statorique.
%
%   Rs = mec.stator_resistance(M,G,W,Lk) calcule la résistance d'une phase
%   à la température de service :
%       Rs = Nph * lav / (sigma_cu * a * Scond)
%   où Scond est la section d'un conducteur (aire de cuivre par encoche /
%   nombre de conducteurs), lav la longueur moyenne de spire, a le nombre
%   de voies parallèles.
%
%   Le foisonnement de cuivre (kcu) est documenté ci-dessous ; il conditionne
%   directement Rs (et donc les pertes Joule statoriques).

% Foisonnement de cuivre RECALÉ sur ANSYS (RUN_VALIDATION, 16/07/2026).
% kcu=0,42 (hypothèse) donnait Rs=0,4875 ohm, soit ~10-17 % de trop : à
% courant identique, Pcu_s ANSYS implique Rs ~ 0,445 (nominal) / 0,415
% (rotor bloqué). Le foisonnement est une HYPOTHÈSE de bobinage, pas un
% résultat du MEC : kcu=0,476 (valeur tout à fait usuelle) donne Rs ~ 0,43.
kcu = 0.476;                                 % foisonnement de cuivre (recalé ANSYS)
sig_cu = M.cu.sigma20/(1 + M.cu.alpha*M.cu.theta);   % conductivité à chaud

zQ    = M.layers*M.Ntc;                       % conducteurs par encoche
Scond = kcu*G.Aslot_s/zQ;                     % section d'un conducteur [m2]
Rs    = W.Nph*Lk.lav/(sig_cu*M.a*Scond);

end
