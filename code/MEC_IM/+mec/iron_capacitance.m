function Cm = iron_capacitance(M, l, A, mu)
%IRON_CAPACITANCE  Capacité magnétique de Perho pour une région de fer (C1).
%
%   Cm = mec.iron_capacitance(M, l, A, mu) renvoie la capacité magnétique
%   C_m d'une branche/région de tube de flux de fer feuilleté, qui embarque
%   les courants de FOUCAULT dans le réseau magnétique (Perho 2002, ch. 5,
%   éq. 79/82/88) :
%
%       Phi = DeltaV_m / R_m  -  C_m * d(DeltaV_m)/dt          (régime réel)
%       Y_m = 1/R_m - j*omega*C_m                              (harmonique)
%
%   Le signe MOINS est imposé par la physique (les courants induits
%   S'OPPOSENT au flux : le flux est réduit et déphasé). C'est aussi ce que
%   donne le développement de la solution exacte 1D d'une tôle (ci-dessous).
%
%   FORMULATION. Perho écrit C_m = 1/(R_z*R_m^2) où R_z est la résistance de
%   la branche axiale court-circuitée (éq. 79), et pour le fer la remplace par
%   une résistance équivalente calée sur le coefficient de pertes du matériau
%   (éq. 87-88). Ici on la dérive directement des grandeurs PHYSIQUES de la
%   tôle (conductivité sigma, épaisseur d), ce qui rend la validation
%   INDÉPENDANTE (cf. RUN_C1) :
%
%       C_m = mu^2 * sigma * d^2 * A / (12 * l)
%
%   Cette forme reproduit exactement le développement basse fréquence de la
%   perméabilité complexe d'une tôle :
%       mu_eff = mu * tanh(gamma*d/2)/(gamma*d/2) ,  gamma = sqrt(j*omega*mu*sigma)
%       mu_eff ~ mu * (1 - j*omega*mu*sigma*d^2/12)          (d << delta_tole)
%   et équivaut à la relation classique  kc*rho = pi^2*sigma*d^2/6.
%
%   ATTENTION (choix assumé) : on utilise sigma et d PHYSIQUES, et NON le
%   coefficient de Bertotti kc de mec.iron_losses. kc est AJUSTÉ sur la perte
%   TOTALE (hystérésis + Foucault + excès) : il surestime les Foucault
%   classiques d'un facteur ~8 sur cette nuance, et l'hystérésis/l'excès ne
%   produisent PAS la même rétroaction sur le flux. Perho le dit : « les pertes
%   par hystérésis restent hors du modèle, d'origine différente ».
%   => C1 modélise la rétroaction des Foucault CLASSIQUES ; les pertes totales
%   restent données par mec.iron_losses (Bertotti).
%
%   Entrées : l (longueur de la branche [m]), A (section perpendiculaire au
%   flux [m2], foisonnement inclus), mu (perméabilité locale [H/m], issue du
%   point de fonctionnement saturé). l, A, mu peuvent être des vecteurs.
%
%   Voir aussi : mec.solve_mesh_complex, RUN_C1, C1_capacite_magnetique_Perho.md

sig = M.fe.sigma;  d = M.fe.d;
Cm = (mu.^2) * sig * d^2 .* A ./ (12*l);
end
