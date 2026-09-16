function P = iron_losses(M, G, S, f)
%IRON_LOSSES  Pertes fer par région à partir du champ MEC (Bertotti).
%
%   P = mec.iron_losses(M,G,S,f) calcule les pertes fer statoriques à
%   partir de la carte d'induction S.Biron (une valeur de B crête par
%   branche de fer) issue de la résolution du réseau, avec le modèle de
%   Bertotti à trois termes :
%       p = kh*f*B^2 + kc*(f*B)^2 + ke*(f*B)^1.5     [W/kg]
%   sommé par branche, pondéré par la masse de fer de la branche et par
%   les coefficients de majoration de région (dents kd, culasses ky).
%
%   Contrairement au post-traitement global de Silva (B d'élément lissé),
%   les pertes sont ici évaluées par région sur le B LOCAL de chaque
%   branche (dents et culasses séparées), ce qui respecte la distribution
%   spatiale de saturation.
%
%   Sorties :
%       P.total, P.teeth_s, P.yoke_s : pertes fer [W]
%       P.byregion : détail
%
%   Note : les pertes fer rotoriques (fréquence de glissement s*f) sont
%   négligeables au voisinage du synchronisme et non comptées ici.

rho = M.fe.rho;
kh = M.iron.kh; kc = M.iron.kc; ke = M.iron.ke;
kd = M.iron.kd; ky = M.iron.ky;

B    = abs(S.Biron(:));
kind = S.kind(:);

% masse de fer par branche : rho * section(kFe) * longueur
Aall = zeros(size(B)); lall = zeros(size(B)); maj = ones(size(B));
% on reconstruit section/longueur par région
for k = 1:numel(B)
    switch kind{k}
        case 'ts', Aall(k)=G.A_ts; lall(k)=G.l_ts; maj(k)=kd;
        case 'ys', Aall(k)=G.A_ys; lall(k)=G.l_ys; maj(k)=ky;
        case 'tr', Aall(k)=G.A_tr; lall(k)=G.l_tr; maj(k)=1;  % rotor ~ negl.
        case 'yr', Aall(k)=G.A_yr; lall(k)=G.l_yr; maj(k)=1;
    end
end
mass = rho .* Aall .* lall;

% densité de pertes de Bertotti [W/kg]
pk = kh*f.*B.^2 + kc*(f.*B).^2 + ke*(f.*B).^1.5;
pbr = maj .* pk .* mass;

isST = strcmp(kind,'ts'); isSY = strcmp(kind,'ys');
P.teeth_s = sum(pbr(isST));
P.yoke_s  = sum(pbr(isSY));
P.total   = P.teeth_s + P.yoke_s;      % stator seul (rotor négligé)
P.mass_teeth_s = sum(mass(isST));
P.mass_yoke_s  = sum(mass(isSY));
end
