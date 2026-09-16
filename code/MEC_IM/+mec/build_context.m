function ctx = build_context(M)
%BUILD_CONTEXT  Assemble tous les objets pré-calculés de la machine.
%
%   ctx = mec.build_context(M) construit une fois pour toutes les
%   structures réutilisées par le solveur de performance :
%       ctx.M,G,W,BH,net,AG,Lk
%       ctx.Rs    : résistance statorique de phase [ohm]
%       ctx.Xm0   : réactance magnétisante NON saturée (référence) [ohm]
%
%   Voir aussi : mec.equivalent_circuit, mec.magnetizing.

ctx.M   = M;
ctx.G   = mec.geometry(M);
ctx.W   = mec.winding(M);
ctx.mat = mec.mat_M800_50A();
ctx.BH  = mec.bh(ctx.mat);
ctx.net = mec.build_network(M, ctx.G);

% ------------------------------------------------------------------
%  ENTREFER : operateur de Dirichlet-to-Neumann de la couronne d'air.
% ------------------------------------------------------------------
%  Chemin par DEFAUT du modele. Il remplace mec.airgap_permeance, qui
%  fermait l'entrefer par un coefficient de Carter ET normalisait les
%  permeances de sorte que leur somme par dent vaille exactement
%  P0 = mu0*L*taus/g_eff : la reactance magnetisante y etait juste PAR
%  CONSTRUCTION, ce qui retirait toute valeur de validation a sa
%  comparaison avec la reference. Avec l'operateur, Xm devient une
%  PREDICTION.
%
%  Les arcs sont ceux des DENTS (Ms = Ns, Mr = Nr), pas des colonnes du
%  maillage fin : le reseau de performance porte un noeud de bec par dent.
%
%  POSITION ROTORIQUE. Contrairement au cas a une seule surface encochee,
%  le bloc Y_sr depend ici de la position relative des deux dentures. Le
%  schema equivalent represente un regime MOYEN : on moyenne donc les
%  grandeurs de sortie sur un pas dentaire (mec.magnetizing), et non la
%  matrice elle-meme. ctx.AG est la position de reference theta = 0 ;
%  ctx.AGpos(k) donne l'operateur aux positions du balayage.
%  ARCS DE SURFACE. Ce sont les FACES DE DENT vues par l'entrefer, donc le
%  pas moins l'ouverture d'encoche -- et NON la largeur du CORPS de dent
%  bts, qui sert au maillage interne. La distinction est decisive : bts ne
%  couvre que 41 % du pas alors que la face en couvre 86 %. Avec les arcs
%  du corps, les noeuds ne pavent pas l'alesage, le mode homopolaire est
%  amoindri (somme des poids 0,41 au lieu de 1) et la surface est vue comme
%  bien plus encochee qu'elle ne l'est. C'est le meme remappage que
%  mec.mesh_refined applique a sa grille de surface.
%  GRILLE DE SURFACE FINE + CONDENSATION. L'operateur ne peut pas etre
%  construit directement sur un noeud par dent : sans degre de liberte au
%  droit de l'ouverture d'encoche, la projection y impose phi = 0 et
%  l'effet d'encoche disparait (coefficient de Carter implicite mesure a
%  1.02 au lieu de 1.27). On pave donc la surface -- nT colonnes de face,
%  nO colonnes d'ouverture -- et on condense sur les noeuds dentaires par
%  complement de Schur, avec Phi = 0 sur les colonnes d'ouverture.
nT = getdef(M,'gap_nT',6);  nO = getdef(M,'gap_nO',2);
ctx.AG  = mec.airgap_dtn_tooth(M, ctx.G, 0, nT, nO, getdef(M,'gap_Nh',[]));
ctx.gap = struct('nT',nT,'nO',nO,'Nh',ctx.AG.Nh,'taur',2*pi/M.Nr);
%  Modele de COMPARAISON, jamais utilise par defaut (voir RUN_ARTICLE §7).
ctx.AGcarter = mec.airgap_permeance(M, ctx.G, 0);

ctx.Lk  = mec.leakage(M, ctx.G, ctx.W);
ctx.Rs  = mec.stator_resistance(M, ctx.G, ctx.W, ctx.Lk);
ctx.H   = mec.harmonics(M, ctx.W, 49);   % spectre d'espace (couples parasites)

% Réactance magnétisante non saturée (petit courant) — référence
R0 = mec.magnetizing(ctx, 0.2);
ctx.Xm0 = R0.Xm;

end

% ---- utilitaire local ----
function v = getdef(M,f,d)
if isfield(M,'opt') && isfield(M.opt,f) && ~isempty(M.opt.(f))
    v = M.opt.(f);
else
    v = d;
end
end
