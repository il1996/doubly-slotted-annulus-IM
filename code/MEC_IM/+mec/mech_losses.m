function Pfw = mech_losses(M, s)
%MECH_LOSSES  Pertes mécaniques (frottement + ventilation).
%
%   Pfw = mec.mech_losses(M,s) estime les pertes mécaniques à un
%   glissement s, avec une loi en carré de la vitesse :
%       Pfw = kfw * Pn * ((1-s))^2
%   kfw (fraction de la puissance nominale à vitesse synchrone) est
%   documenté ci-dessous. Modèle simple mais suffisant : les pertes
%   mécaniques ne dépendent que de la vitesse, pas de la charge.

% kfw RECALÉ sur ANSYS (RUN_VALIDATION, 16/07/2026) : la valeur initiale 0,008
% (hypothèse) donnait 142 W contre 98 W mesurés par ANSYS (+45 %). ANSYS
% fournit Pro = 98,2 W a n = 1471,8 tr/min, d'ou
%     kfw = 98,2 / (Pn * (1-s)^2) = 0,0055
% C'est une DONNÉE D'ENTRÉE (frottement + ventilation), pas un résultat du
% MEC : elle doit venir de l'essai ou de l'EF, comme ici.
kfw = 0.0055;                % ~0,55 % de Pn à vitesse synchrone (recalé ANSYS)
Pfw = kfw * M.Pn * (1-s).^2;
end
