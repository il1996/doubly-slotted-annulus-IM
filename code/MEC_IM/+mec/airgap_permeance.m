function AG = airgap_permeance(M, G, theta_rot)
%AIRGAP_PERMEANCE  Perméances d'entrefer dent-à-dent [MODÈLE DE COMPARAISON].
%
%   ------------------------------------------------------------------
%   AVERTISSEMENT — CETTE FONCTION N'EST PLUS SUR LE CHEMIN ACTIF.
%   Le modèle ferme l'entrefer par un COEFFICIENT DE CARTER (mec.carter)
%   et normalise ensuite les perméances de sorte que leur somme par dent
%   vaille exactement P0 = mu0*L*tau_s/g_eff. La réactance magnétisante
%   qui en sort est donc juste PAR CONSTRUCTION : sa comparaison avec la
%   référence EF ne valide rien. C'est précisément ce que l'opérateur de
%   Dirichlet-to-Neumann (mec.airgap_dtn_tooth) remplace, et l'écart
%   entre les deux est le résultat rapporté en §7 de RUN_ARTICLE.
%   Ne la rebranchez sur ctx.AG que pour PRODUIRE cette comparaison.
%   ------------------------------------------------------------------
%
%   AG = mec.airgap_permeance(M,G,theta_rot) construit les perméances de
%   couplage entre chaque dent statorique i et les dents rotoriques j
%   qu'elle voit, pour une position mécanique rotorique theta_rot [rad].
%
%   PRINCIPE (correction du défaut M1 de Silva 2023) :
%   la fonction de couplage est paramétrée en grandeurs PHYSIQUES —
%   entrefer g, ouvertures d'encoche bs0/br0, largeurs de becs de dents —
%   et NON en arcs d'éléments de maillage. Elle est donc invariante au
%   raffinement tangentiel (test : Ns,Nr x2 -> Xm et couple stables).
%
%   Modèle :
%     * entrefer effectif  g_eff = kC * g   (coefficient de Carter,
%       ouvertures stator ET rotor) ;
%     * perméance totale attribuée à une dent stator :
%           P0 = mu0 * L * tau_s / g_eff        (conservation de flux) ;
%     * répartition sur les dents rotoriques en regard par un noyau de
%       recouvrement trapézoïdal (becs de dents, largeur physique
%       tau - b0) ÉLARGI d'une frange de portée lambda_fr ~ g (flux de
%       frange gouverné par l'entrefer, non par le maillage) ;
%     * normalisation par dent stator => somme(P_ij) = P0 exactement
%       (Xm correct par construction), noyau C^1 (couple dérivable).
%
%   Sorties :
%     AG.i, AG.j, AG.P : listes (branches d'entrefer) stator i <-> rotor j
%     AG.kC            : coefficient de Carter total
%     AG.g_eff         : entrefer effectif [m]
%     AG.P0            : perméance moyenne par dent stator [H]
%     AG.Psum_stator   : contrôle (doit valoir ~P0 pour chaque i)
%
%   Voir aussi : mec.build_network, mec.carter.

mu0 = 4*pi*1e-7;
Ns = M.Ns;  Nr = M.Nr;  L = M.L;  g = M.g;

% ------------------------------------------------------------------
%  1. Coefficient de Carter (ouvertures stator et rotor)
% ------------------------------------------------------------------
kCs = mec.carter(G.taus, M.bs0, g);
kCr = mec.carter(G.taur, M.br0, g);
kC  = kCs * kCr;
g_eff = kC * g;
AG.kC = kC;  AG.g_eff = g_eff;

% ------------------------------------------------------------------
%  2. Positions périphériques (au rayon moyen d'entrefer)
% ------------------------------------------------------------------
Rg = (G.Rs + G.Rr)/2;                 % rayon moyen d'entrefer
Cper = 2*pi*Rg;                       % circonférence moyenne
xs = ((0:Ns-1)+0.5) * (Cper/Ns);      % centres de dents stator [m]
xr = ((0:Nr-1)+0.5) * (Cper/Nr) + theta_rot*Rg;   % centres dents rotor [m]

% Largeurs physiques des becs de dents (fer face à l'entrefer)
w_s = G.taus - M.bs0;                 % bec de dent stator [m]
w_r = G.taur - M.br0;                 % bec de dent rotor  [m]
lambda_fr = pi*g;                     % portée physique de frange ~ g [m]

% Perméance totale par dent stator (conservation)
P0 = mu0 * L * G.taus / g_eff;
AG.P0 = P0;

% ------------------------------------------------------------------
%  3. Noyau de recouvrement + frange, par couple (i,j)
% ------------------------------------------------------------------
% On ne connecte que les dents rotoriques proches (fenêtre +/- 1.5 pas).
win = 1.5 * max(G.taus, G.taur);      % demi-fenêtre de recherche [m]

ii = []; jj = []; PP = [];
Psum = zeros(Ns,1);
for i = 1:Ns
    % distance périphérique signée (repli circulaire) vers chaque dent rotor
    d = xr - xs(i);
    d = mod(d + Cper/2, Cper) - Cper/2;      % dans [-C/2, C/2)
    near = find(abs(d) <= win);
    if isempty(near), continue; end
    k = overlap_kernel(d(near), w_s, w_r, lambda_fr);
    ks = sum(k);
    if ks <= 0, continue; end
    w = k / ks;                              % poids normalisés (somme=1)
    P = P0 * w;                              % perméances (somme=P0)
    ii = [ii; i*ones(numel(near),1)];        %#ok<AGROW>
    jj = [jj; near(:)];                      %#ok<AGROW>
    PP = [PP; P(:)];                         %#ok<AGROW>
    Psum(i) = sum(P);
end

AG.i = ii;  AG.j = jj;  AG.P = PP;
AG.Psum_stator = Psum;
AG.Rg = Rg;

end

% ======================================================================
function k = overlap_kernel(d, w_s, w_r, lambda)
%OVERLAP_KERNEL  Recouvrement trapézoïdal de deux becs + frange gaussienne.
%   d      : distances périphériques (vecteur) [m]
%   w_s,w_r: largeurs des becs stator/rotor [m]
%   lambda : portée de frange [m]
d = abs(d(:));
% recouvrement de deux rectangles de largeurs w_s et w_r :
ov = min((w_s + w_r)/2 - d, min(w_s, w_r));
ov = max(ov, 0);
% frange au-delà du contact (transition C^1) :
gap = max(d - (w_s + w_r)/2, 0);
fr  = min(w_s, w_r) * exp(-(gap./lambda).^2);
k = ov + fr;
end
