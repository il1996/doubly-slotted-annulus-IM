function M = machine_18_5kW()
%MACHINE_18_5KW  Jeu de paramètres complet d'un moteur asynchrone à cage 18,5 kW.
%
%   M = mec.machine_18_5kW() renvoie une structure regroupant TOUS les
%   paramètres primaires de la machine (géométriques, électriques,
%   magnétiques, mécaniques). Toute grandeur dérivée est calculée dans
%   mec.geometry ; ce fichier ne contient que des données d'entrée.
%
%   Machine de référence : moteur asynchrone triphasé à cage,
%   4 pôles, 48 encoches statoriques / 44 barres rotoriques, tôles
%   M800-50A, alimenté sous 690 V (couplage étoile), 50 Hz.
%   Géométrie issue du modèle EF validé (emdlab / ANSYS Maxwell 2D).
%
%   Ancrages FEM du projet (à retrouver) :
%     Bg1 ~ 1,01 T ; I0 ~ 8,32 A ; Xm ~ 46 ohm ; Rs ~ 0,101 ohm.
%
%   Convention d'unités : tout est en SI (m, A, V, T, H, ohm, kg) sauf
%   mention contraire ; les cotes primaires sont saisies en mm puis
%   converties ci-dessous.
%
%   Voir aussi : mec.geometry, mec.mat_M800_50A, RUN_MEC_IM.

mm = 1e-3;

% ----------------------------------------------------------------------
%  1. Plaque signalétique et alimentation
% ----------------------------------------------------------------------
M.name   = 'IM 18.5 kW - 690 V - 4 poles - 48/44';
M.Pn     = 18.5e3;      % puissance nominale utile [W]
M.Ull    = 690;         % tension composée (ligne-ligne) [V]
M.connection = 'Y';     % couplage étoile
M.Uph    = M.Ull/sqrt(3);   % tension simple (phase-neutre) [V]
M.f      = 50;          % fréquence d'alimentation [Hz]
M.m      = 3;           % nombre de phases
M.p      = 2;           % nombre de paires de pôles (=> 4 pôles)
M.w      = 2*pi*M.f;    % pulsation électrique [rad/s]
M.ns     = M.f/M.p;     % vitesse synchrone [tr/s]  (=25 tr/s = 1500 tr/min)

% ----------------------------------------------------------------------
%  2. Géométrie primaire (paquet, entrefer, diamètres)
% ----------------------------------------------------------------------
M.L      = 164.782448*mm;   % longueur active du paquet [m]
M.g      = 0.243*mm;        % entrefer mécanique [m]  (cote réelle fournie)
M.Ds     = 163.7775*mm;     % diamètre d'alésage statorique (bore) [m]
M.Dso    = 266.5506*mm;     % diamètre extérieur statorique [m]
M.Dsh    = 74.2665*mm;      % diamètre d'arbre [m]
M.kFe    = 0.97;            % coefficient de foisonnement des tôles [-]

% ----------------------------------------------------------------------
%  3. Encoches statoriques  (cotes réelles — profil « tc6 »)
% ----------------------------------------------------------------------
%   Profil (de l'entrefer vers la culasse) :
%     - ouverture (isthme)  : largeur Bs0 sur la hauteur Hs0
%     - biseau (wedge)      : Bs0 -> Bs1 sur la hauteur Hs1
%     - corps trapézoïdal   : Bs1 (haut) -> Bs2 (fond) sur la hauteur Hs2
%                             (l'encoche S'ÉLARGIT vers le fond)
%     - fond                : coins arrondis de rayon Rs (ajoute Rs à la hauteur)
%   NB : l'élargissement compense exactement la croissance du pas d'encoche
%        => la DENT STATORIQUE EST À FLANCS PARALLÈLES (cf. mec.geometry).
M.Ns     = 48;              % nombre d'encoches statoriques
M.bs0    = 2.000*mm;        % largeur de l'ouverture d'encoche (isthme) [m]
M.bs1    = 5.2360*mm;       % largeur d'encoche côté entrefer [m]
M.bs2    = 8.4724*mm;       % largeur d'encoche côté fond [m]
M.Rs     = M.bs2/3;         % rayon d'arrondi de fond d'encoche [m]
M.hs0    = 0.50*mm;         % hauteur de l'isthme (ouverture) [m]
M.hs1    = 2.50*mm;         % hauteur du biseau (wedge) [m]
M.hs2    = 24.7242*mm;      % hauteur droite de l'encoche [m]

% ----------------------------------------------------------------------
%  4. Encoches / barres rotoriques (cotes réelles — profil « tcr », goutte)
% ----------------------------------------------------------------------
%   Profil (de l'entrefer vers l'arbre), CONFORME A LA FIGURE ANSYS (19/07) :
%     - ouverture (isthme)  : largeur Br0 sur la hauteur Hr0 (+ Hr01)
%     - ÉPAULEMENT          : le cercle supérieur (diamètre Br1) est TANGENT
%       au bas de l'ouverture ; son centre est à Br1/2 sous celle-ci
%     - partie effilée      : tangente aux deux cercles, entraxe Hr1
%     - cercle inférieur    : diamètre Br2  (l'encoche SE RÉTRÉCIT vers le fond)
M.Nr     = 44;              % nombre de barres rotoriques
M.br0    = 2.000*mm;        % largeur de l'ouverture d'encoche rotorique [m]
M.br1    = 5.8043*mm;       % diamètre du cercle supérieur de barre [m]
M.br2    = 2.0382*mm;       % diamètre du cercle inférieur de barre [m]
M.hr0    = 1.00*mm;         % hauteur de l'isthme rotorique [m]
M.hr01   = 0.00*mm;         % hauteur additionnelle d'isthme [m]
M.hr1    = 16.7898*mm;      % entraxe des deux cercles de barre [m]

% ----------------------------------------------------------------------
%  5. Bobinage statorique (double couche, pas raccourci 5/6)
% ----------------------------------------------------------------------
M.Ntc    = 18;      % nombre de conducteurs (spires) par côté de bobine
M.layers = 2;       % nombre de couches
M.a      = 2;       % nombre de voies parallèles par phase
M.yq     = 10;      % pas de bobinage [pas d'encoche]  (plein = 12 => 5/6)
% Motif d'encoches (machine complète) reproduisant le modèle EF :
%   couche du bas (down) et couche du haut (up), phase A ; B et C décalés.
pA1 = [1 2 3 4 13 14 15 16];  pA1 = [pA1, pA1+24];      % down, phase A
sgn = [1 1 1 1 -1 -1 -1 -1];  sgn = [sgn, sgn];          % signes (aller/retour)
w48 = @(x) mod(x-1,48)+1;
pA2 = [47 48 1 2 11 12 13 14]; pA2 = [pA2, w48(pA2+24)]; % up, phase A
M.wind.pDown = pA1;   M.wind.sDown = sgn;
M.wind.pUp   = pA2;   M.wind.sUp   = sgn;
M.wind.shift = [0 8 4];   % décalage d'encoches phases [A B C] ; C inversée
M.wind.Cinv  = true;      % la phase C est connectée en inverse (modèle EF)

% ----------------------------------------------------------------------
%  6. Cage rotorique
% ----------------------------------------------------------------------
M.ring.ler   = 0.0*mm;      % jeu axial paquet <-> anneau (un côté) [m]
M.ring.frac  = 1.0;         % fraction de la section de barre reportée à l'anneau
% (section et diamètre moyen de l'anneau calculés dans mec.geometry)

% ----------------------------------------------------------------------
%  7. Propriétés des matériaux conducteurs
% ----------------------------------------------------------------------
M.cu.sigma20 = 57e6;    % conductivité cuivre à 20 degC [S/m]
M.cu.alpha   = 3.81e-3; % coeff. de température du cuivre [1/K]
M.cu.theta   = 80;      % échauffement de service du bobinage [K]
M.cu.rho     = 8900;    % masse volumique du cuivre [kg/m3]

M.al.sigma20 = 37e6;    % conductivité aluminium (barres) à 20 degC [S/m]
M.al.alpha   = 3.7e-3;  % coeff. de température de l'aluminium [1/K]
M.al.theta   = 80;      % échauffement de service de la cage [K]
M.al.rho     = 2700;    % masse volumique de l'aluminium [kg/m3]

M.fe.rho     = 7700;    % masse volumique des tôles [kg/m3]
M.fe.sigma   = 1.03e6;  % conductivité électrique des tôles [S/m] (M800-50A)
M.fe.d       = 0.50e-3; % épaisseur d'une tôle [m] (nuance ...-50A = 0,50 mm)
% NB : sigma et d gouvernent les courants de Foucault CLASSIQUES (et donc la
% rétroaction sur le flux, cf. C1 / capacité magnétique de Perho). Ils sont
% distincts des coefficients de Bertotti ci-dessous, qui sont AJUSTÉS sur la
% perte totale (hystérésis + Foucault + excès) et ne sont pas utilisables pour
% la rétroaction. Voir C1_capacite_magnetique_Perho.md.

% ----------------------------------------------------------------------
%  8. Coefficients de pertes fer (Bertotti, M800-50A)  et pertes méca
% ----------------------------------------------------------------------
%   pFe = kh*f*B^2 + kc*(f*B)^2 + ke*(f*B)^1.5     [W/kg]
%   Calés pour ~8 W/kg à 1,5 T / 50 Hz (nuance M800-50A) :
%   répartition hystérésis 60 % / Foucault 30 % / excès 10 %.
M.iron.kh = 0.043;      % coeff. hystérésis  [W/(kg) / (T^2 Hz)]
M.iron.kc = 4.3e-4;     % coeff. courants de Foucault [W/(kg)/(T Hz)^2]
M.iron.ke = 1.2e-3;     % coeff. pertes en excès [W/(kg)/(T Hz)^1.5]
% Majorations de région RECALÉES sur ANSYS (RUN_VALIDATION, 16/07/2026).
% Valeurs initiales kd=1,8 / ky=1,6 (reprises de Perho 2002, AUTRE lot de
% tôles — faiblesse signalée par le défaut P3/C2) : elles surestimaient les
% pertes fer de +34 % (312 vs 233 W en charge) et +30 % (325 vs 249 à vide).
% Recalage global par le facteur 0,755 = moyenne des deux rapports ANSYS/MEC,
% en conservant le rapport kd/ky. NB : avec deux points de fonctionnement à
% flux quasi identique, on ne peut PAS séparer kh/kc/ke ni distinguer un
% excès de coefficients de Bertotti d'un excès de majoration : c'est un
% recalage GLOBAL assumé, pas une identification des trois termes.
M.iron.kd = 1.36;       % majoration pertes dents (recalé ANSYS) [-]
M.iron.ky = 1.21;       % majoration pertes culasses (recalé ANSYS) [-]

M.mech.kfric = 15;      % facteur de pertes mécaniques (Table 9.2 Pyrhönen)
% pMech ~ kfric * (D/... ) : formule appliquée dans mec.losses.

% Dynamique mécanique (démarrage transitoire) :
M.mech.J   = 0.170;     % inertie rotor + arbre + ventilateur [kg.m2]
% (calé sur le run-up EF a vide ~0.13 s ; fer rotor seul ~0.09, le reste
%  arbre/anneaux/ventilateur/accouplement)
M.mech.B   = 0.010;     % coefficient de frottement visqueux [N.m.s/rad]
M.mech.TL0 = 0.0;       % couple de charge à vide [N.m]

% ----------------------------------------------------------------------
%  8bis. Pertes supplémentaires en charge (stray load losses)
% ----------------------------------------------------------------------
%  Terme de pertes EN CHARGE non capturé par le réseau fondamental (pulsation
%  de denture, surface, courants inter-barres, extrémités 3D) : sans lui le
%  MEC donne 94,1 % au nominal (irréaliste). Le modèle de champ pas-à-pas
%  (RUN_STEPPING) n'explique que ~17,5 W (4,6 %) de ces pertes par la physique
%  de denture ; les mécanismes dominants (surface/inter-barres/3D) sont hors
%  d'un modèle 2D fondamental. Introduites comme CONSTANTE de pertes suppl. =
%  ALLOCATION IEC 60034-2-1 : Pn_ref = 381,5 W (~2 % Pin pour cette taille),
%  valeur COMMUNE aux trois méthodes de l'étude comparative (analytique
%  PLL = 381,5 W ; résidu non identifié de l'EF ANSYS = 381,5 W) -> rendements
%  comparés à hypothèse de pertes suppl. IDENTIQUE (choix utilisateur 20/07).
%  NB : à 381,5 W le rendement MEC vaut ~92,1 %, AU-DESSUS de la plage
%  expérimentale 90,5-91,3 %, car les pertes IDENTIFIÉES du MEC sont ~215 W
%  sous l'EF ; un calage à Pn_ref ~ 606 W placerait le rendement au milieu de
%  la plage (90,9 %). Loi quadratique en courant de charge (IEEE 112 / IEC
%  60034-2-1), plafonnée à 150 % de charge (cf. mec.stray_losses).
M.stray.Pn_ref   = 381.5;    % pertes suppl. de référence au nominal [W] (IEC 60034-2-1)
M.stray.I2_ref   = 16.463;   % courant rotor rapporté de référence [A rms]
M.stray.load_cap = 1.5;      % plafond de charge (150 %) de la loi quadratique [-]

% ----------------------------------------------------------------------
%  9. Options de résolution (documentées ; valeurs par défaut)
% ----------------------------------------------------------------------
M.opt.newton_tol   = 1e-9;   % tolérance sur le résidu de bilan de flux [Wb]
M.opt.newton_itmax = 60;     % itérations Newton maxi
M.opt.newton_relax = 1.0;    % amortissement (line-search auto si <1 diverge)
M.opt.sat_tol      = 1e-4;   % tolérance point fixe MEC<->circuit sur Xm [-]
M.opt.sat_itmax    = 40;     % itérations maxi du point fixe magnétisant
M.opt.nbar_layers  = 6;      % couches radiales de la barre (effet de peau)
M.opt.verbose      = true;

end
