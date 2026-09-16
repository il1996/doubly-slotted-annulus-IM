function G = geometry(M)
%GEOMETRY  Grandeurs géométriques dérivées de la machine.
%
%   G = mec.geometry(M) calcule, à partir des cotes primaires de M
%   (mec.machine_*), toutes les grandeurs géométriques nécessaires au
%   réseau de réluctances et aux modèles de fuites/pertes :
%   diamètres, pas d'encoche, largeurs et hauteurs de dents, épaisseurs
%   de culasses, sections et longueurs de tubes de flux (dents, culasses),
%   section de barre et d'anneau, longueurs de fer.
%
%   Toutes les sections de fer intègrent le foisonnement kFe.
%   Les largeurs de dents sont évaluées à mi-hauteur d'encoche
%   (valeur représentative d'un tube de flux radial d'élément) — grandeur
%   PHYSIQUE, indépendante de toute discrétisation (principe d'invariance).
%
%   Voir aussi : mec.build_network, mec.leakage, mec.cage.

G = struct();
L    = M.L;
kFe  = M.kFe;
G.L  = L;      % longueur active du paquet [m]

% ---------------- Diamètres et rayons de référence --------------------
G.Ds   = M.Ds;                  % alésage stator
G.Rs   = M.Ds/2;                % rayon d'alésage
G.Dr   = M.Ds - 2*M.g;          % diamètre extérieur rotor
G.Rr   = G.Dr/2;
G.Dso  = M.Dso;
G.Dsh  = M.Dsh;

% ---------------- Hauteurs totales d'encoches (cotes réelles) ---------
% Stator : isthme + biseau + corps + arrondi de fond (rayon Rs)
G.hs   = M.hs0 + M.hs1 + M.hs2 + M.Rs;
% Rotor (goutte, CONSTRUCTION CONFORME A LA FIGURE ANSYS du 19/07/2026) :
% l'ouverture Br0 debouche sur un EPAULEMENT : le cercle superieur (diam.
% Br1) est TANGENT au bas de l'ouverture, son centre est donc a Br1/2 sous
% celle-ci (et non a sqrt((Br1/2)^2-(Br0/2)^2), variante « cercle recoupe »
% utilisee avant le 19/07 : delta +0,18 mm sur hr).
G.hc_r = M.br1/2;                             % ouverture -> centre cercle sup.
G.hc_r_cut = sqrt((M.br1/2)^2 - (M.br0/2)^2); % (ancienne variante, trace)
G.hr   = M.hr0 + M.hr01 + G.hc_r + M.hr1 + M.br2/2;

% ---------------- Pas d'encoche (au niveau de l'entrefer) -------------
G.taus = pi*G.Ds/M.Ns;          % pas d'encoche statorique [m]
G.taur = pi*G.Dr/M.Nr;          % pas d'encoche rotorique  [m]

% ---------------- Largeurs de dents (géométrie réelle) ----------------
% STATOR : l'élargissement Bs1->Bs2 compense la croissance du pas d'encoche
% => dent à FLANCS PARALLÈLES. On le vérifie aux deux extrémités du corps.
r_top   = G.Ds/2 + M.hs0 + M.hs1;                 % haut du corps d'encoche
r_bot   = G.Ds/2 + M.hs0 + M.hs1 + M.hs2;         % fond du corps d'encoche
bt_top  = 2*pi*r_top/M.Ns - M.bs1;
bt_bot  = 2*pi*r_bot/M.Ns - M.bs2;
G.bts   = 0.5*(bt_top + bt_bot);                  % largeur de dent stator [m]
G.bts_parallel_err = abs(bt_top-bt_bot)/G.bts;    % contrôle (~0 si parallèle)

% ROTOR : la dent s'élargit vers l'arbre (l'encoche en goutte se rétrécit).
% Largeur EFFECTIVE au sens de la réluctance : moyenne harmonique entre le
% centre du cercle supérieur et celui du cercle inférieur (les deux extrêmes).
r_up    = G.Dr/2 - (M.hr0 + M.hr01 + G.hc_r);     % centre cercle supérieur
r_lo    = r_up - M.hr1;                            % centre cercle inférieur
bt_up   = 2*pi*r_up/M.Nr - M.br1;                  % dent la plus ÉTROITE
bt_lo   = 2*pi*r_lo/M.Nr - M.br2;                  % dent la plus large
G.btr   = 2/(1/bt_up + 1/bt_lo);                   % moyenne harmonique [m]
G.btr_min = bt_up;  G.btr_max = bt_lo;

% ---------------- Épaisseurs de culasses ------------------------------
G.hys  = M.Dso/2 - (G.Ds/2 + G.hs);           % épaisseur culasse stator [m]
G.hyr  = (G.Dr/2 - G.hr) - M.Dsh/2;           % épaisseur culasse rotor  [m]

% ---------------- Diamètres moyens de culasses ------------------------
G.Dys  = M.Dso - G.hys;                        % diamètre moyen culasse stator
G.Dyr  = M.Dsh + G.hyr;                        % diamètre moyen culasse rotor

% ---------------- Longueurs de fer ------------------------------------
G.lfe  = kFe * L;               % longueur de fer effective [m]

% ---------------- Sections des tubes de flux (fer) --------------------
% Dents : section = largeur de dent x longueur de fer.
G.A_ts = G.bts * G.lfe;         % section d'une dent stator [m2]
G.A_tr = G.btr * G.lfe;         % section d'une dent rotor  [m2]
% Culasses : section = épaisseur de culasse x longueur de fer.
G.A_ys = G.hys * G.lfe;         % section de culasse stator [m2]
G.A_yr = G.hyr * G.lfe;         % section de culasse rotor  [m2]

% ---------------- Longueurs des tubes de flux -------------------------
% Dents : longueur radiale = hauteur d'encoche.
G.l_ts = G.hs;                  % longueur d'une dent stator [m]
G.l_tr = G.hr;                  % longueur d'une dent rotor  [m]
% Culasses : longueur d'un segment = arc entre deux dents adjacentes.
G.l_ys = pi*G.Dys/M.Ns;         % segment de culasse stator [m]
G.l_yr = pi*G.Dyr/M.Nr;         % segment de culasse rotor  [m]

% ---------------- Sections d'encoche / de barre (géométrie réelle) ----
% ENCOCHE STATOR (aire disponible pour le cuivre) :
%   biseau Bs0->Bs1 sur Hs1  +  corps trapézoïdal Bs1->Bs2 sur Hs2
%   +  fond : bande Bs2 x Rs dont les DEUX coins sont arrondis au rayon Rs
%   (chaque coin retire Rs^2*(1-pi/4)).
G.Aslot_s = 0.5*(M.bs0 + M.bs1)*M.hs1 ...
          + 0.5*(M.bs1 + M.bs2)*M.hs2 ...
          + M.bs2*M.Rs - 2*M.Rs^2*(1 - pi/4);              % [m2]

% BARRE ROTOR (goutte, figure ANSYS) : le cercle superieur etant TANGENT au
% bas de l'ouverture (epaulement), il n'est PLUS ampute d'un segment :
% demi-cercle superieur PLEIN + trapeze d'entraxe Hr1 + demi-cercle inferieur.
r1 = M.br1/2;  r2 = M.br2/2;
G.Abar = 0.5*pi*r1^2 + 0.5*pi*r2^2 + (r1 + r2)*M.hr1;       % [m2]
G.Aseg_r = 0;                                  % plus de segment (epaulement)

% ---------------- Anneau de court-circuit -----------------------------
% Diamètre moyen de l'anneau ~ diamètre au fond des barres.
G.Dring = G.Dr - 2*G.hr + G.hr;      % ~ diamètre moyen barre
G.Dring = G.Dr - G.hr;               % (diamètre moyen barre, plus robuste)
% Section d'anneau : dimensionnée pour porter le courant d'anneau ;
% règle usuelle Aring ~ (Nr/(pi*p)) * Abar (Pyrhönen) — ajustable.
G.Aring = (M.Nr/(pi*M.p)) * G.Abar * M.ring.frac;  % [m2]
% Longueur d'un segment d'anneau entre deux barres.
G.lring_seg = pi*G.Dring/M.Nr;       % [m]

% ---------------- Ouvertures d'encoche (pour Carter et fuites) --------
G.bs0 = M.bs0;    G.br0 = M.br0;
G.g   = M.g;

% ---------------- Aires de pôle (contrôle) ----------------------------
G.taup = pi*G.Ds/(2*M.p);            % pas polaire [m]
G.Apole = G.taup * L;                % aire d'un pôle à l'entrefer [m2]

end
