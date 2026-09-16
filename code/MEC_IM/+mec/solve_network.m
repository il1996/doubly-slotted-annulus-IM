function S = solve_network(net, G, BH, AG, Fs, Fr, opt)
%SOLVE_NETWORK  Résolution non linéaire du réseau de réluctances (Newton).
%
%   S = mec.solve_network(net,G,BH,AG,Fs,Fr,opt) résout le réseau
%   magnétique pour :
%       net : topologie fer (mec.build_network)
%       BH  : interpolants matériau (mec.bh)
%       AG  : perméances d'entrefer (mec.airgap_permeance)
%       Fs  : FMM de dents statoriques (Ns x 1) [A]  (bobinage)
%       Fr  : FMM de dents rotoriques  (Nr x 1) [A]  (barres) ou []
%       opt : options (tolérance, itérations)
%
%   MÉTHODE — Newton-Raphson EXACT sur les potentiels nodaux :
%   pour chaque branche de fer, la loi d'Ampère locale s'écrit
%       Delta U_eff = H(B) * l ,   B = Phi/A ,
%   donc B = Bof(DeltaU_eff / l) est explicite (courbe B(H) monotone).
%   La réluctance DIFFÉRENTIELLE dPhi/dDeltaU = A/(l*dH/dB) fournit le
%   jacobien analytique J = Ainc' * diag(Pdiff) * Ainc (Perho 2002).
%   Ce solveur remplace le point fixe sous-relaxé de Silva (kr=0,05-0,3,
%   100-250 it.) par un Newton amorti qui converge en ~5-10 itérations,
%   avec critère d'arrêt sur le RÉSIDU de bilan de flux (défaut B2), pas
%   sur la simple stagnation de B.
%
%   Sorties :
%       S.U        : potentiels nodaux [A]
%       S.Biron    : induction dans chaque branche de fer [T]
%       S.Phiron   : flux dans chaque branche de fer [Wb]
%       S.kind     : région de chaque branche de fer
%       S.Pgap     : flux dans chaque branche d'entrefer [Wb]
%       S.Bgap_i   : induction radiale moyenne par dent stator [T]
%       S.iter, S.res : diagnostic de convergence
%
%   Voir aussi : mec.build_network, mec.bh, mec.airgap_permeance.

if nargin < 7 || isempty(opt), opt = struct(); end
tol   = getfielddef(opt,'newton_tol',1e-9);
itmax = getfielddef(opt,'newton_itmax',60);

Nn = net.Nnodes;
if isempty(Fr), Fr = zeros(net.Nr,1); end

% ------------------------------------------------------------------
%  Assemblage des listes de branches (fer + entrefer)
% ------------------------------------------------------------------
% Fer
ia = net.iron.a; ib = net.iron.b;
il = net.iron.l; iA = net.iron.A;
nFe = net.iron.n;
% Source FMM par branche de fer
E = zeros(nFe,1);
mS = net.iron.srcS>0; E(mS) = Fs(net.iron.srcS(mS));
mR = net.iron.srcR>0; E(mR) = Fr(net.iron.srcR(mR));

% ------------------------------------------------------------------
%  Entrefer : DEUX formes admises
% ------------------------------------------------------------------
%  (a) FORME DENSE (defaut) — AG.Y est l'operateur de Dirichlet-to-Neumann
%      de la couronne d'air (mec.airgap_fourier). Ce n'est pas un jeu de
%      branches a deux bornes : c'est un bloc de couplage qui relie TOUS
%      les noeuds de surface entre eux. Il ne peut donc pas passer par la
%      matrice d'incidence et s'ajoute DIRECTEMENT sur le sous-bloc
%      [surface x surface] du systeme nodal.
%
%      SIGNE. L'en-tete d'airgap_fourier pose "flux SORTANT des noeuds =
%      AF.Y*[Us;Ur]" ; la loi des noeuds s'ecrit ici "somme des flux
%      sortants nulle" et les branches de fer contribuent par
%      Phi = G*(Ua-Ub+F) avec incidence +1 en a. Le bloc s'ajoute donc
%      AVEC LE MEME SIGNE que les conductances de branche. Une erreur ici
%      converge vers une solution fausse sans rien signaler : elle est
%      verrouillee par le test T4 de RUN_VERIFY_OPERATOR (limites
%      analytiques de la couronne).
%
%      ORDRE DES NOEUDS. AF.Y suppose [stator ; rotor]. Toute permutation
%      silencieuse produit un couplage croise faux.
%
%  (b) FORME LISTES — AG.i / AG.j / AG.P, permeances dent-a-dent de
%      mec.airgap_permeance. Conservee UNIQUEMENT comme modele de
%      comparaison (elle contient le coefficient de Carter et normalise
%      les permeances de sorte que Xm soit juste par construction).
denseGap = isfield(AG,'Y') && ~isempty(AG.Y);

if denseGap
    %  net.idx.TS / TR sont des FONCTIONS d'indexation, pas des tableaux :
    %  on les evalue sur la liste complete des dents.
    iTS = net.idx.TS((1:net.Ns).');
    iTR = net.idx.TR((1:net.Nr).');
    idxSurf = [iTS(:); iTR(:)];
    MsSurf  = numel(iTS);
    if size(AG.Y,1) ~= numel(idxSurf)
        error('mec:solve_network:gapSize', ...
            ['AG.Y est %dx%d mais le reseau porte %d noeuds de surface ' ...
             '(%d stator + %d rotor).'],size(AG.Y,1),size(AG.Y,2), ...
            numel(idxSurf),MsSurf,numel(iTR));
    end
    nAg  = 0;
    allA = ia; allB = ib;
else
    ga = net.idx.TS(AG.i);      % bec dent stator
    gb = net.idx.TR(AG.j);      % bec dent rotor
    gP = AG.P;
    nAg = numel(gP);
    allA = [ia; ga(:)];
    allB = [ib; gb(:)];
end

% Incidence globale (Nb x Nn), branche k : +1 en a, -1 en b
Nb   = nFe + nAg;
rows = (1:Nb).';
Ainc = sparse([rows;rows],[allA;allB],[ones(Nb,1);-ones(Nb,1)],Nb,Nn);

% Masque fer / air
isFe = [true(nFe,1); false(nAg,1)];
Eall = [E; zeros(nAg,1)];

% Réduction : on retire le nœud de référence
keep = true(Nn,1); keep(net.ref) = false;
kk = find(keep);
Ared = Ainc(:,kk);

% Bloc dense reduit : le noeud de reference doit etre retire de la LIGNE
% ET de la COLONNE du bloc, comme il l'est du reste du systeme.
if denseGap
    Ybig = sparse(Nn,Nn);
    Ybig(idxSurf,idxSurf) = AG.Y;               %#ok<SPRIX>
    Yred = Ybig(kk,kk);
else
    Yred = sparse(numel(kk),numel(kk));
end
%  Flux SOURCE du bloc d'entrefer (courant de barre dans les cavites
%  d'encoche, mec.airgap_dtn_tooth_cav.set_source) : flux sortant des noeuds
%  de surface independant des potentiels, ajoute au bilan de flux avec le
%  meme signe que AG.Y*U. Absent (vide) dans toutes les chaines anterieures.
fsrc = zeros(numel(kk),1);
if denseGap && isfield(AG,'f') && ~isempty(AG.f)
    fbig = zeros(Nn,1); fbig(idxSurf) = AG.f(:); fsrc = fbig(kk);
end

% ------------------------------------------------------------------
%  Boucle de Newton amortie
% ------------------------------------------------------------------
U = zeros(Nn,1);              % init : potentiels nuls
res = inf; it = 0;
mu0 = 4*pi*1e-7;
while it < itmax
    it = it + 1;
    dU_eff = Ainc*U + Eall;                 % Delta U_eff par branche

    Phi   = zeros(Nb,1);
    Pdiff = zeros(Nb,1);

    % --- branches de fer (non linéaires) ---
    Hfe = dU_eff(isFe)./il;                  % champ H = DeltaU/l
    Bfe = BH.Bof(Hfe);                       % induction B(H)
    Phi(isFe)   = Bfe .* iA;
    dHdB        = BH.dHdB(Bfe);              % réluctance diff. (pente)
    Pdiff(isFe) = iA ./ (il .* max(dHdB,1e-12));

    % --- branches d'entrefer (linéaires) ---
    if ~denseGap
        Phi(~isFe)   = gP .* dU_eff(~isFe);
        Pdiff(~isFe) = gP;
    end

    % Résidu de bilan de flux (KCL magnétique) sur nœuds libres
    r = Ared.' * Phi + Yred*U(kk) + fsrc;
    res = max(abs(r));
    if res < tol, break; end

    % Jacobien et pas de Newton. Le bloc d'entrefer est LINEAIRE : il entre
    % dans le jacobien tel quel, sans reassemblage.
    J = Ared.' * spdiags(Pdiff,0,Nb,Nb) * Ared + Yred;
    dx = -(J \ r);

    % Amortissement (line-search simple sur la norme du résidu)
    lam = 1.0;
    Ufull = U;
    rnorm0 = norm(r);
    for ls = 1:20
        Utry = U;  Utry(kk) = U(kk) + lam*dx;
        dUe = Ainc*Utry + Eall;
        Ph = zeros(Nb,1);
        Ph(isFe) = BH.Bof(dUe(isFe)./il).*iA;
        if ~denseGap, Ph(~isFe) = gP.*dUe(~isFe); end
        rt = Ared.'*Ph + Yred*Utry(kk) + fsrc;
        if norm(rt) < rnorm0 || lam < 1e-4
            Ufull = Utry; break;
        end
        lam = lam/2;
    end
    U = Ufull;
end

% ------------------------------------------------------------------
%  Post-traitement
% ------------------------------------------------------------------
dU_eff = Ainc*U + Eall;
Bfe = BH.Bof(dU_eff(isFe)./il);
S.U       = U;
S.Biron   = Bfe;
S.Phiron  = Bfe .* iA;
S.kind    = net.iron.kind;
S.iter    = it;
S.res     = res;
S.converged = res < max(tol*1e3, 1e-6);
S.denseGap = denseGap;

Vol_fe = il .* iA;                              % volume de chaque branche fer
Wco_fe = sum( BH.wco(Bfe) .* Vol_fe );

Ns = net.Ns;
if denseGap
    %  Flux d'entrefer : il n'y a plus de branche a deux bornes. Le flux
    %  SORTANT de chaque noeud de surface est directement Y*U_surf ; celui
    %  qui traverse l'entrefer depuis une dent statorique est donc la
    %  composante correspondante du produit.
    Usurf = U(idxSurf);
    Fsurf = AG.Y*Usurf;
    if isfield(AG,'f') && ~isempty(AG.f), Fsurf = Fsurf + AG.f(:); end   % flux source (courant de barre)
    S.Usurf = Usurf;
    S.Us    = Usurf(1:MsSurf);
    S.Ur    = Usurf(MsSurf+1:end);
    Phi_gap_stator = Fsurf(1:MsSurf);
    S.Pgap  = Phi_gap_stator;
    %  Co-energie de la couronne : forme quadratique du bloc, 1/2 U'YU.
    Wco_ag = 0.5*(Usurf.'*Fsurf);
else
    S.Pgap = gP .* dU_eff(~isFe);       % flux d'entrefer par branche
    Wco_ag = 0.5 * sum( S.Pgap .* dU_eff(~isFe) );
    Phi_gap_stator = accumarray(AG.i(:), S.Pgap, [Ns,1]);
end
S.Wco = Wco_fe + Wco_ag;
S.Phi_gap_stator = Phi_gap_stator;
A_tip_stator = (G.taus - G.bs0)*G.L;     % section de bec stator (air, sans kFe)
S.Bgap_i     = Phi_gap_stator ./ A_tip_stator;   % B au bec de dent (concentré)
S.Bgap_avg_i = Phi_gap_stator ./ (G.taus*G.L);   % B moyen d'entrefer (pas d'encoche)

end

% ---- utilitaires locaux ----
function v = getfielddef(s,f,d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
