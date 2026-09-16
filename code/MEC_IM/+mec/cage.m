function C = cage(M, G, W, Lk, s, nu, kw_nu)
%   [Variante harmonique] mec.cage(M,G,W,Lk,s_nu,nu,kw_nu) renvoie l'impédance
%   de cage vue par l'harmonique d'espace d'ordre nu (nu*p paires de pôles,
%   facteur de bobinage kw_nu) : le report utilise kw_nu au lieu de kw1, et la
%   contribution d'anneau 1/(2*sin^2(pi*nu*p/Nr)) au lieu de 1/(2*sin^2(pi*p/Nr)).
%   Indispensable aux couples parasites (mec.harmonics / RUN_HARMONICS) : le
%   report FONDAMENTAL y serait faux.
%CAGE  Impédance de la cage rotorique reportée au stator (effet de peau).
%
%   C = mec.cage(M,G,W,Lk,s) calcule, pour un glissement s :
%       C.Rr     : résistance rotorique reportée au stator [ohm]
%       C.Xr     : réactance de fuite rotorique reportée [ohm]
%       C.kR,C.kX: facteurs d'effet de peau (barre profonde)
%       C.Rbar, C.Rring : résistances physiques barre / segment d'anneau
%       C.xi     : hauteur de barre réduite h/delta
%
%   MODÈLE — cage barre + anneaux, avec effet de peau analytique :
%     * profondeur de peau a la fréquence rotorique  delta = sqrt(2/(w_r*mu0*sigma))
%     * facteurs de barre profonde (Field/Emde) :
%         kR = xi (sinh2xi + sin2xi)/(cosh2xi - cos2xi)
%         kX = (3/2xi)(sinh2xi - sin2xi)/(cosh2xi - cos2xi)
%       kR majore la résistance de barre (refoulement du courant en tête),
%       kX minore la perméance de fuite d'encoche rotorique.
%     * report au stator : facteur 4*m*(kw1*Nph)^2 / Nr  (cage = Nr phases
%       à 1 conducteur), anneaux ramenés par 1/(2 sin^2(pi*p/Nr)).
%
%   Cette formulation flux-de-dent / courant-de-barre (Ostović 5.2.1)
%   évite l'ambiguïté d'une « FMM de barre » et respecte la dualité
%   FMM<->FEM (cohérence énergétique). Règle de discrétisation associée
%   (si passage à une cage multicouche dans le réseau) : >= 3 couches par
%   profondeur de peau a la fréquence maximale.
%
%   Voir aussi : mec.leakage, mec.equivalent_circuit.

mu0 = 4*pi*1e-7;
m = M.m; Nr = M.Nr; p = M.p; L = M.L;
Nph = W.Nph; kw1 = W.kw1;
% --- variante harmonique (defaut = fondamental) ---
if nargin < 6 || isempty(nu),    nu = 1;      end
if nargin < 7 || isempty(kw_nu), kw_nu = kw1; end
kw_use = abs(kw_nu);          % facteur de bobinage de l'harmonique
pn     = nu*p;                % paires de poles de l'harmonique

% Conductivités a chaud
sig_al = M.al.sigma20/(1 + M.al.alpha*M.al.theta);

% Facteur de réduction d'inclinaison (skew) rotor — C5
%  Vrillage réel = 1 pas d'encoche STATOR (ssp=1, soit 7,5° mécaniques).
%  Pour l'harmonique d'espace nu, l'angle de vrillage ÉLECTRIQUE est
%  multiplié par nu :  ksq_nu = sin(nu*asq)/(nu*asq)  ->  ~0,997 au
%  fondamental, ~0,04 pour nu=23/25 : les harmoniques de denture ne se
%  couplent quasiment plus à la cage — c'est précisément le rôle du
%  vrillage (leur branche rotorique devient quasi ouverte : la branche
%  harmonique redevient une réactance quasi pure, sans perte ni couple).
%  Désactivable par M.opt.skew_harm = 0 (comparaison avec une référence
%  EF NON vrillée).
taus = G.taus; tp = pi*G.Ds/(2*p);
asq1 = 1*(taus/tp)*pi/2;            % demi-angle électrique fondamental
skewh = ~(isfield(M,'opt') && isfield(M.opt,'skew_harm') ...
          && ~isempty(M.opt.skew_harm) && M.opt.skew_harm==0);
if skewh, asq = nu*asq1; else, asq = asq1; end
ksq = sin(asq)/asq;
C.ksq = ksq;

% ---------------- Effet de peau (fréquence rotorique) -----------------
%  Résolution 1D EXACTE pour la forme RÉELLE (trapézoïdale) de l'encoche
%  rotorique (mec.bar_skin), au lieu des formules de Field-Emde qui supposent
%  une barre RECTANGULAIRE. Le haut de barre étant large (br1=5,80) et le fond
%  étroit (br2=2,04), le refoulement du courant coûte moins cher :
%  kR(trapèze) < kR(rectangle) — jusqu'à -7 % a s=1 (cf. RUN_BARSKIN).
%  mec.bar_skin est validé a 0,00 % contre Field-Emde sur un profil rectangulaire.
wr = 2*pi*M.f*abs(s);
if wr < 1e-9
    xi = 0; kR = 1; kX = 1;
else
    Sk = mec.bar_skin(M, G, abs(s), struct('shape','trapz','N',600));
    kR = Sk.kR; kX = Sk.kX; xi = Sk.xi;
end
C.xi = xi; C.kR = kR; C.kX = kX;

% ---------------- Résistances physiques -------------------------------
Lbar = L + 2*M.ring.ler;                       % longueur de barre
Rbar_dc = Lbar/(sig_al*G.Abar);                % résistance DC de barre
Rbar = Rbar_dc*kR;                             % avec effet de peau
% Segment d'anneau entre deux barres
Rring_seg = (pi*G.Dring/Nr)/(sig_al*G.Aring);
% ---- P2 : effet de peau de l'ANNEAU a la frequence rotorique ----------
%  L'anneau est un conducteur massif EN AIR (pas en encoche) : le champ de
%  son propre courant penetre par les DEUX faces -> modele 1D bilateral de
%  Field-Emde applique a la DEMI-hauteur, xi_ring = (h/2)/delta, avec
%  h = hauteur radiale de la section (~ hauteur de barre, l'anneau faisant
%  face aux extremites de barres ; la largeur axiale Aring/h est plus
%  grande). Negligeable en zone utile (delta(s=0,02) ~ 90 mm >> h) ; actif a
%  s=1 (delta = 12,7 mm) et surtout aux frequences des branches HARMONIQUES
%  (|s_nu|*f ~ 1,2 kHz -> delta ~ 2,6 mm). ANSYS le resout en conducteur
%  massif (RingSolidLoss) ; kX d'anneau non corrige (part interne faible
%  du forfait lam_ring).
if wr < 1e-9
    kRing = 1;
else
    delta_r = sqrt(2/(wr*mu0*sig_al));
    xr2 = 0.5*G.hr/delta_r;                 % (h/2)/delta, h ~ hauteur de barre
    kRing = xr2*(sinh(2*xr2)+sin(2*xr2))/max(cosh(2*xr2)-cos(2*xr2),1e-12);
    kRing = max(kRing,1);
end
C.kRing = kRing;
C.Rring_dc = Rring_seg;
Rring_seg = Rring_seg*kRing;
% Contribution d'anneau ramenée a une barre (pn = nu*p paires de pôles)
sr = sin(pi*pn/Nr)^2;
if sr < 1e-9, sr = 1e-9; end            % garde (harmonique en phase avec la cage)
Rring_eq = Rring_seg/(2*sr);
C.Rbar = Rbar; C.Rring = Rring_seg;
C.Rbar_dc = Rbar_dc;

% Résistance rotorique par « phase de cage » (barre + anneaux)
Rrot_bar = Rbar + Rring_eq;

% Report au stator (kw_use = kw1 au fondamental, kw_nu pour un harmonique)
kref = (4*m/Nr)*(kw_use*Nph)^2;
C.Rr = kref * Rrot_bar / ksq^2;

% ---------------- Fuite rotorique reportée ----------------------------
% perméance de fuite : encoche (skin -> kX) + bec/zig-zag + anneau
lam_ring = 0.66;                               % surplomb/anneau (Table 4.1)
lam_r = Lk.lam_r_slot*kX + Lk.lam_r_tip + lam_ring;
% Report cohérent de la fuite (même facteur que la résistance) :
%   Lsigma_r' = (4*m/Nr)*mu0*L*(kw_use*Nph)^2 * lambda_r
Lr_prime = (4*m/Nr)*mu0*L*(kw_use*Nph)^2*lam_r/ksq^2;
C.Xr = M.w*Lr_prime;
C.lam_r = lam_r;

% Ajout de la fuite différentielle rotor (proportionnelle a Xm, ajoutée
% dans le circuit équivalent) : on renvoie juste le coefficient.
C.sigd_r = Lk.sigd_r;

end
