function id = fem_slot_leakage(M, G, opt)
%FEM_SLOT_LEAKAGE  Identification FEM de la perméance de fuite d'encoche.
%
%   id = mec.fem_slot_leakage(M,G,opt) identifie par éléments finis
%   magnétostatiques 2D (PDE Toolbox) la perméance de fuite lambda d'une
%   encoche (rotor ou stator), à comparer aux formules analytiques
%   (Pyrhönen) utilisées par mec.leakage. Même démarche que
%   mec.fem_airgap_ident pour sigma0 : on mesure au lieu de supposer.
%
%   MOTIVATION : le diagnostic RUN_LEAKAGE montre que la réactance de fuite
%   du modèle est 8-13 % TROP FAIBLE au fort glissement (Xsig_MEC=3,26 vs
%   Xsig_EF=3,74 a s=0,5). Le maillon faible identifié est le terme
%   empirique 0.4*(br1/br0) ajouté dans lam_r_tip. On le remplace ici par
%   une identification.
%
%   BANC : le domaine est l'ENCOCHE SEULE (barre + isthme). Le fer qui
%   l'entoure est supposé infiniment perméable, ce qui impose sur ses parois
%   B TANGENTIEL NUL (le champ entre perpendiculairement au fer) — soit une
%   condition de NEUMANN homogène sur A (naturelle). Le flux traverse donc
%   l'encoche horizontalement et se referme par le fer EXTÉRIEUR au domaine,
%   qui ne coûte rien : c'est exactement l'hypothèse de la perméance d'encoche
%   classique. Une seule condition de Dirichlet (A=0 en haut, coté entrefer)
%   fixe la jauge et exclut le flux d'entrefer.
%   La perméance est extraite par l'ENERGIE (indépendante du chemin) :
%       L' = 2*W/I^2   [H/m]        lambda = L'/mu0     [-]
%
%   VALIDATION : sur une encoche RECTANGULAIRE, lambda doit valoir h/(3b)
%   (cf. test_slotleak). NB : un banc qui enferme le flux (A=0 sur TOUT le
%   contour) est FAUX — le flux y boucle en traversant l'encoche deux fois et
%   lambda est surestimé de +9 a +60 % (erreur croissante avec h/b).
%
%   opt.side  : 'rotor' (défaut) | 'stator'
%   opt.mur   : perméabilité relative du fer (défaut 5000)
%   opt.Hmax  : taille de maille (défaut ouverture/4)
%
%   Sortie id : .lambda (FEM), .L (H/m), .W (J/m), .nelem, .lambda_ana
%   (valeur analytique du modèle, pour comparaison).
%
%   Voir aussi : mec.leakage, RUN_LEAKAGE, mec.fem_airgap_ident.

if nargin<3, opt=struct(); end
side = getdef(opt,'side','rotor');
mur  = getdef(opt,'mur',5000);
mu0  = 4*pi*1e-7;

% ---- cotes de l'encoche selon le coté ----
switch lower(side)
    case 'rotor'
        tau=G.taur; b0=M.br0; h0=M.hr0; b1=M.br1; b2=M.br2; hb=M.hr1;
        % analytique du modele (corps + isthme + pont), hors zig-zag
        lam_ana = hb/(3*0.5*(b1+b2)) + h0/b0 + 0.4*(b1/b0);
    case 'stator'
        tau=G.taus; b0=M.bs0; h0=M.hs0+M.hs1; b1=M.bs2; b2=M.bs1; hb=M.hs2;
        lam_ana = M.hs2/(3*0.5*(M.bs1+M.bs2)) + M.hs1/((M.bs0+M.bs1)/2) + M.hs0/M.bs0;
    otherwise, error('side inconnu');
end
Hmax  = getdef(opt,'Hmax', b0/4);
yb0 = 0; yb1 = hb; Hd = hb + h0;            % encoche seule : barre puis isthme

% ---- géométrie (decsg) : barre (trapèze) + isthme ; PAS de fer maillé ----
model = createpde('electromagnetic','magnetostatic');
model.VacuumPermeability = mu0;
BAR = [2;4; -b2/2; b2/2; b1/2; -b1/2;  yb0; yb0; yb1; yb1];   % trapèze
IST = [3;4; -b0/2; b0/2; b0/2; -b0/2;  yb1; yb1; Hd; Hd];     % isthme
n = max(numel(BAR),numel(IST));
BAR(end+1:n)=0; IST(end+1:n)=0;
gd = [BAR, IST];
ns = char('BAR','IST')';
dl = decsg(gd,'BAR+IST',ns);
geometryFromEdges(model,dl);
generateMesh(model,'Hmax',Hmax,'GeometricOrder','linear');

% ---- air partout (le fer est HORS domaine) + source dans la barre ----
me = model.Mesh; Nodes = me.Nodes;
I = 1;                                    % courant unitaire (linéaire)
Abar = 0.5*(b1+b2)*hb;
Jbar = I/Abar;
electromagneticProperties(model,'RelativePermeability',1);
for f = 1:model.Geometry.NumFaces
    el = findElements(me,'region','Face',f);
    if isempty(el), continue; end
    nn = me.Elements(:,el(1));
    xc = mean(Nodes(1,nn)); yc = mean(Nodes(2,nn));
    if inbar(xc,yc)
        electromagneticSource(model,'CurrentDensity',Jbar,'Face',f);
    end
end

% ---- CL : A=0 sur la seule arête du HAUT (jauge + exclusion du flux
%      d'entrefer). Partout ailleurs : Neumann naturel = B tangentiel nul,
%      soit la paroi de fer infiniment perméable.
eb = [];
for e = 1:model.Geometry.NumEdges
    nid = findNodes(me,'region','Edge',e);
    if all(abs(Nodes(2,nid)-Hd)<1e-9), eb(end+1)=e; end %#ok<AGROW>
end
electromagneticBC(model,'MagneticPotential',0,'Edge',eb);

Rz = solve(model);
Bx = Rz.MagneticFluxDensity.Bx(:); By = Rz.MagneticFluxDensity.By(:);

% ---- énergie : W = integrale B^2/(2*mu) ----
P = Rz.Mesh.Nodes; T = Rz.Mesh.Elements(1:3,:);
xe = reshape(P(1,T),3,[]); ye = reshape(P(2,T),3,[]);
Ar = 0.5*abs((xe(2,:)-xe(1,:)).*(ye(3,:)-ye(1,:))-(xe(3,:)-xe(1,:)).*(ye(2,:)-ye(1,:)));
B2 = mean(Bx(T).^2 + By(T).^2, 1);
W = sum( B2./(2*mu0) .* Ar );         % J/m (encoche = air uniquement)

id.W = W;  id.L = 2*W/I^2;  id.lambda = id.L/mu0;
id.lambda_ana = lam_ana;  id.nelem = size(T,2);
id.side = side;

    function b = inbar(x,y)
        if y<yb0 || y>yb1, b=false; return; end
        t = (y-yb0)/max(yb1-yb0,eps);
        bw = b2 + (b1-b2)*t;
        b = abs(x) <= bw/2;
    end
end

function v=getdef(s,f,d), if isfield(s,f)&&~isempty(s.(f)),v=s.(f);else,v=d;end, end
