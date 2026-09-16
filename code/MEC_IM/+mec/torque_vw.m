function Tvw = torque_vw(ctx, res)
%TORQUE_VW  Couple électromagnétique par travaux virtuels (champ chargé).
%
%   Tvw = mec.torque_vw(ctx, res) calcule le couple à partir du CHAMP MEC
%   chargé (stator + rotor excités) par la méthode des travaux virtuels :
%       T = dW_co/dtheta   à courants constants,
%   évaluée par différence finie de la co-énergie magnétique du réseau
%   entre deux positions rotoriques voisines. C'est la voie recommandée
%   (dérivée de perméance d'entrefer) plutôt que le tenseur de Maxwell
%   monocouche (défaut M6 de Silva) — indépendante du schéma équivalent,
%   elle sert de CONTRE-VALIDATION du couple de circuit res.Tem.
%
%   DOMAINE DE VALIDITÉ : validé au voisinage du point de fonctionnement
%   (faible glissement, écart < 6 % vs circuit et vs EF). À fort glissement
%   la reconstruction de la phase spatiale rotorique et l'entrefer
%   monocouche (1 élément/dent — défaut M6 de Silva) dégradent l'estimation
%   (le couple de circuit reste la référence). Amélioration : entrefer
%   C1 + travaux virtuels analytiques dLambda/dtheta (A5), maillage
%   tangentiel raffiné (A2).
%
%   Reconstruction du champ chargé à l'instant t=0 :
%     * FMM statorique : courants instantanés des 3 phases (phaseur I1) ;
%     * FMM rotorique  : onde de barres reconstruite pour reproduire la
%       FMM fondamentale rotorique (phaseur I2), avec la phase spatiale
%       relative issue des phaseurs.
%
%   Voir aussi : mec.equivalent_circuit, mec.solve_network.

M = ctx.M; G = ctx.G; W = ctx.W; BH = ctx.BH; net = ctx.net;
p = M.p; Nr = M.Nr;

% ---- FMM statorique (instant t=0) ----
psi1 = angle(res.I1c);
i3 = sqrt(2)*abs(res.I1c)*[cos(psi1); cos(psi1-2*pi/3); cos(psi1+2*pi/3)];
Fs = W.slotMMF(i3);

% ---- FMM rotorique : onde de barres calée sur le phaseur I2 ----
psi2 = angle(res.I2c);
thr  = 2*pi*(0:Nr-1)/Nr;                      % angles mécaniques des barres
ib   = cos(p*thr - psi2);                     % courants de barre (unité)
Fr_u = cumsum(ib(:)); Fr_u = Fr_u - mean(Fr_u);
% fondamental de la FMM rotorique unitaire
c    = (2/Nr)*sum(Fr_u(:).'.*exp(-1j*p*thr));
% cible : FMM fondamentale rotorique (formule stator-référée)
Fr1  = (3/2)*(4/pi)*(W.kw1*W.Nph/(2*p))*sqrt(2)*abs(res.I2c);
Fr   = -Fr_u * (Fr1/abs(c));                  % signe : rotor s'oppose (charge)

% ---- Couple moyen : <dW_co/dtheta> sur un pas d'encoche (anti-ondulation)
% Le couple instantané par travaux virtuels est bruité par la denture
% (entrefer monocouche, 1 élément/dent — défaut M6). On moyenne la
% dérivée de co-énergie sur un pas d'encoche statorique.
dth = 0.1*pi/180;                             % pas de différence finie [rad]
Np  = 12;                                     % positions sur un pas d'encoche
th0 = linspace(0, 2*pi/M.Ns, Np+1); th0(end) = [];
Tk = zeros(Np,1);
for kk = 1:Np
    AGp = mec.airgap_permeance(M, G, th0(kk)+dth);
    AGm = mec.airgap_permeance(M, G, th0(kk)-dth);
    Sp = mec.solve_network(net, G, BH, AGp, Fs, Fr, M.opt);
    Sm = mec.solve_network(net, G, BH, AGm, Fs, Fr, M.opt);
    Tk(kk) = (Sp.Wco - Sm.Wco)/(2*dth);
end
Tvw = mean(Tk);                               % periodicite = 1 (machine complete)
end
