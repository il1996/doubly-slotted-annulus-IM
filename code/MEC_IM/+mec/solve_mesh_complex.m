function S = solve_mesh_complex(mesh, omega, opt)
%SOLVE_MESH_COMPLEX  Réseau de réluctances COMPLEXE avec capacité magnétique (C1).
%
%   S = mec.solve_mesh_complex(mesh, omega, opt) résout, en régime
%   harmonique à la pulsation omega, un réseau de réluctances dont les
%   branches de fer portent une CAPACITÉ MAGNÉTIQUE (Perho ch. 5) qui
%   embarque les courants de Foucault DANS le réseau magnétique — sans
%   circuit électrique additionnel. C'est l'implémentation de C1, qui lève
%   le défaut P6 (Silva déclare la chose impossible ; sa propre réf. [16] la fait).
%
%   Admittance magnétique par branche :
%       fer  :  Y = 1/R_m - j*omega*C_m        (Foucault : flux réduit/déphasé)
%       air  :  Y = P                          (réel)
%   Le système nodal est LINÉAIRE (perméabilité figée au point de
%   fonctionnement) : une seule résolution complexe, pas de Newton.
%
%   mesh (champs requis, en plus de ceux de mec.solve_mesh) :
%       .Rm  : réluctance de chaque branche de fer [A/Wb]   (= nu*l/A)
%       .Cm  : capacité magnétique de chaque branche de fer [Wb/A/s]
%              (0 pour les branches sans Foucault)
%
%   Sorties : S.U (potentiels complexes), S.Phi (flux complexes),
%   S.Peddy (perte Foucault par branche [W]), S.Peddy_tot,
%   S.Y (admittances). La perte d'une branche vaut
%       P = 1/2 * C_m * omega^2 * |DeltaV_m|^2      (DeltaV_m en amplitude)
%   (dérivée de Perho éq. 77/79 ; cf. C1_capacite_magnetique_Perho.md §4ter).
%
%   NB : les pertes par hystérésis et par excès ne sont PAS dans ce modèle
%   (origine différente — Perho §5.3) ; elles restent données par
%   mec.iron_losses. C1 apporte la RÉTROACTION des Foucault sur le flux.
%
%   Voir aussi : mec.iron_capacitance, mec.solve_mesh, RUN_C1.

if nargin<3 || isempty(opt), opt=struct(); end
Nn = mesh.Nnodes;  Nb = numel(mesh.a);
a = mesh.a(:); b = mesh.b(:);
isFe = logical(mesh.iron(:));
E = mesh.E(:);

% --- admittances de branche ---
Y = zeros(Nb,1);
Y(~isFe) = mesh.P(~isFe);                       % air : perméance réelle
Rm = mesh.Rm(:); Cm = mesh.Cm(:);
Y(isFe) = 1./Rm(isFe) - 1i*omega*Cm(isFe);      % fer : Foucault (signe -)

% --- incidence et assemblage nodal ---
rows = (1:Nb).';
Ainc = sparse([rows;rows],[a;b],[ones(Nb,1);-ones(Nb,1)],Nb,Nn);
keep = true(Nn,1); keep(mesh.ref) = false; kk = find(keep);
Ared = Ainc(:,kk);
Gy = spdiags(Y,0,Nb,Nb);

% Phi = Y.*(Ainc*U + E) ; bilan nodal : Ared'*Phi = 0
K = Ared.' * Gy * Ared;
rhs = -(Ared.' * (Y.*E));
U = zeros(Nn,1);
U(kk) = K \ rhs;

% --- post-traitement ---
dUe = Ainc*U + E;                                % DeltaV_m effectif (complexe)
Phi = Y .* dUe;
Peddy = zeros(Nb,1);
Peddy(isFe) = 0.5 * Cm(isFe) * omega^2 .* abs(dUe(isFe)).^2;

S.U = U;  S.Phi = Phi;  S.dU = dUe;  S.Y = Y;
S.Peddy = Peddy;  S.Peddy_tot = sum(Peddy);
S.omega = omega;
end
