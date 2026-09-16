function Lk = leakage(M, G, W)
%LEAKAGE  Inductances de fuite statoriques et coefficients différentiels.
%
%   Lk = mec.leakage(M,G,W) calcule les composantes de fuite magnétique
%   NON prises en compte par le réseau d'entrefer principal :
%     * fuite d'encoche statorique (perméance de fente)
%     * fuite de bec de dent / zig-zag (tooth-tip)
%     * fuite de têtes de bobines (end-winding)
%     * fuite différentielle (harmoniques d'espace) via un coefficient
%       sigma_delta appliqué à Xm dans le circuit équivalent.
%   et de même les perméances de fuite GÉOMÉTRIQUES rotoriques (hors
%   effet de peau, traité dans mec.cage).
%
%   Formulation classique (Pyrhönen, Machines électriques) :
%       L_sigma = (4*m/Q) * mu0 * L * Nph^2 * lambda
%   pour chaque composante lambda (perméance spécifique adimensionnelle).
%
%   Sorties :
%       Lk.Xs_slot, Lk.Xs_tip, Lk.Xs_ew : réactances de fuite stator [ohm]
%       Lk.Xs_leak  : fuite stator hors différentielle [ohm]
%       Lk.sigd_s   : coefficient de fuite différentielle stator [-]
%       Lk.lam_r_slot, Lk.lam_r_tip : perméances de fuite rotor (géom.)
%       Lk.sigd_r   : coefficient de fuite différentielle rotor [-]
%
%   Voir aussi : mec.cage, mec.equivalent_circuit.

mu0 = 4*pi*1e-7;
m = M.m; Ns = M.Ns; Nr = M.Nr; p = M.p; L = M.L; g = M.g;
Nph = W.Nph; kw1 = W.kw1; q = W.q;

% Facteur de raccourcissement pour la fuite d'encoche (double couche)
beta = M.yq/(Ns/(2*p));            % pas relatif (5/6)
k1   = (1 + 3*beta)/4;             % correction double couche (partie basse)

% ---------------- Fuite d'encoche statorique --------------------------
%  P4 : permeance IDENTIFIEE PAR FEM (mec.fem_slot_leakage / RUN_SLOTLEAK,
%  banc Neumann valide a 0,0-0,2 % vs h/(3b)) : lam_FEM = 2.722 a courant
%  UNIFORME (k1=1), contre 2.143 pour les formules analytiques (qui
%  sous-estiment le biseau reel). La correction double couche k1 ne
%  s'applique qu'a la part de CORPS (ou loge le cuivre).
bs_avg = (M.bs1 + M.bs2)/2;
lam_body_u = M.hs2/(3*bs_avg);                 % corps, courant uniforme
lam_s_FEM  = 2.722;                            % IDENTIFIE (RUN_SLOTLEAK)
lam_s_slot = k1*lam_body_u + (lam_s_FEM - lam_body_u);
Ls_slot = (4*m/Ns)*mu0*L*Nph^2*lam_s_slot;

% ---------------- Fuite de bec de dent / zig-zag (stator) -------------
lam_s_tip = (5*(g/M.bs0))/(5 + 4*(g/M.bs0));
Ls_tip = (4*m/Ns)*mu0*L*Nph^2*lam_s_tip;

% ---------------- Fuite de têtes de bobines ---------------------------
tp   = pi*G.Ds/(2*p);
Wew  = tp*(M.yq/(Ns/(2*p)));       % largeur d'une bobine (portée)
lw   = 1.2*Wew + 0.05;             % longueur d'une tête (empirique 5.2)
lew  = lw - Wew;
lam_lew = 0.50; lam_w = 0.20;      % perméances de têtes (Table 4.1)
lam_w_av = (2*lew*lam_lew + Wew*lam_w)/lw;
Ls_ew = (4*m/Ns)*mu0*q*Nph^2*lw*lam_w_av;

% ---------------- Fuite d'extremite 3D identifiee (P4) ----------------
%  Les bancs d'encoche 2D ne voient ni la geometrie reelle des tetes de
%  bobines ni l'environnement 3D des anneaux : le RESIDU est identifie
%  GLOBALEMENT sur les Xsigma impliques par l'EF a s = 0.2/0.5/1.0
%  (RUN_P4_XSIGMA — seul le TOTAL est identifiable depuis les grandeurs
%  terminales ; il est porte cote stator par convention).
Lk.Xs_end3D = 0.426;                % [ohm] IDENTIFIE (RUN_P4_XSIGMA :
                                    %  delta 0.478/0.471/0.301 a s=0.2/0.5/1 ;
                                    %  geometrie d'encoches figure ANSYS 19/07)
if isfield(M,'opt') && isfield(M.opt,'Xs_end3D') && ~isempty(M.opt.Xs_end3D)
    Lk.Xs_end3D = M.opt.Xs_end3D;   % override (etudes de sensibilite)
end

Lk.Xs_slot = M.w*Ls_slot;
Lk.Xs_tip  = M.w*Ls_tip;
Lk.Xs_ew   = M.w*Ls_ew;
Lk.Xs_leak = M.w*(Ls_slot + Ls_tip + Ls_ew) + Lk.Xs_end3D;
Lk.lw = lw;  Lk.lav = 2*(L + lw);   % longueur moyenne de spire (résistance)

% ---------------- Fuite différentielle (harmonique) stator ------------
Lk.sigd_s = diff_leakage_coeff(q, beta, p, Ns, kw1);

% ---------------- Perméances de fuite rotor (géométriques) ------------
%  P4 : encoche rotor IDENTIFIEE PAR FEM (RUN_SLOTLEAK) : total barre +
%  isthme = 1.543. L'isthme vaut EXACTEMENT hr0/br0 = 0.500 (champ uniforme
%  au-dessus de la barre) => CORPS reel = 1.043 — la formule a largeur
%  moyenne (1.427) surestime le trapeze a haut large. Le terme empirique
%  0.4*(br1/br0) = 1.161 N'EXISTE PAS (le FEM implique -0.384) : c'etait
%  une rustine qui masquait la fuite d'extremite 3D (voir Xs_end3D).
lam_r_FEM = 1.543;                             % IDENTIFIE (RUN_SLOTLEAK)
lam_r_ist = M.hr0/M.br0;                       % isthme (exact)
Lk.lam_r_slot = lam_r_FEM - lam_r_ist;         % corps de barre (skin -> kX)
Lk.lam_r_tip  = lam_r_ist ...
              + (5*(g/M.br0))/(5 + 4*(g/M.br0));% zig-zag rotor
% Fuite différentielle rotor (harmoniques de barres)
Lk.sigd_r = 0.5*(pi*p/Nr)^2 / 3;               % approx. cage (Pyrhönen)

end

% ======================================================================
function sigd = diff_leakage_coeff(q, beta, p, Ns, kw1)
%DIFF_LEAKAGE_COEFF  Coefficient de fuite différentielle d'un bobinage.
%   Somme des harmoniques d'espace ν = 6k±1 pondérés par (kw_nu/nu)^2.
sigd = 0;
for nu = [5 7 11 13 17 19 23 25]
    % facteur de distribution/raccourcissement de l'harmonique nu
    alse = 2*pi*p/Ns;
    kdn = sin(nu*q*alse/2)/(q*sin(nu*alse/2));
    kpn = sin(nu*beta*pi/2);
    kwn = kdn*kpn;
    sigd = sigd + (kwn/(nu*kw1))^2;
end
end
