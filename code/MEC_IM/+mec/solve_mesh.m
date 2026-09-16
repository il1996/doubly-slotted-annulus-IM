function S = solve_mesh(mesh, BH, opt)
%SOLVE_MESH  Solveur de Newton générique sur un maillage de réluctances.
%
%   S = mec.solve_mesh(mesh, BH, opt) résout un réseau de réluctances décrit
%   par des LISTES DE BRANCHES génériques (fer non linéaire + air linéaire),
%   par le même Newton-Raphson exact (réluctance différentielle dH/dB) que
%   mec.solve_network, mais sans hypothèse sur la topologie. Il sert le
%   maillage polaire raffinable (mec.mesh_refined) — implémentation du
%   mailleur non-paramétrique absent de la toolbox de Silva (A3).
%
%   mesh (champs requis) :
%       Nnodes         : nombre de nœuds
%       ref            : nœud de référence (potentiel nul)
%       a, b           : nœuds extrémités de chaque branche (Nb x 1)
%       iron           : logique, true = branche de fer non linéaire
%       l, A           : longueur/section des branches de fer
%       P              : perméance des branches d'air
%       E              : source de FMM par branche [A] (0 si aucune)
%
%   Sorties : S.U (potentiels), S.B (induction des branches fer),
%   S.Phi (flux par branche), S.Wco (co-énergie), S.iter, S.res.

if nargin < 3 || isempty(opt), opt = struct(); end
tol   = getdef(opt,'newton_tol',1e-9);
itmax = getdef(opt,'newton_itmax',80);

Nn = mesh.Nnodes;  Nb = numel(mesh.a);
a = mesh.a(:); b = mesh.b(:);
isFe = logical(mesh.iron(:));
l = mesh.l(:); A = mesh.A(:); P = mesh.P(:); E = mesh.E(:);

% Incidence (Nb x Nn) : +1 en a, -1 en b
rows = (1:Nb).';
Ainc = sparse([rows;rows],[a;b],[ones(Nb,1);-ones(Nb,1)],Nb,Nn);

keep = true(Nn,1); keep(mesh.ref) = false;
kk = find(keep);
Ared = Ainc(:,kk);

% admittance lineaire additionnelle (entrefer harmonique P1, mesh.Yx)
hasYx = isfield(mesh,'Yx') && ~isempty(mesh.Yx);
if hasYx, Yx = mesh.Yx; Yxr = Yx(kk,:); Yxkk = Yx(kk,kk); end
% source nodale optionnelle (aimant du DtN etendu, mesh.Isrc) — inactif si absent
hasIs = isfield(mesh,'Isrc') && ~isempty(mesh.Isrc);
if hasIs, Isr = mesh.Isrc(kk); end

U = zeros(Nn,1); res = inf; it = 0;
lFe = l(isFe); AFe = A(isFe); PAir = P(~isFe);
while it < itmax
    it = it + 1;
    dUe = Ainc*U + E;
    Phi = zeros(Nb,1); Pd = zeros(Nb,1);
    Hfe = dUe(isFe)./lFe;
    Bfe = BH.Bof(Hfe);
    Phi(isFe) = Bfe.*AFe;
    Pd(isFe)  = AFe./(lFe.*max(BH.dHdB(Bfe),1e-12));
    Phi(~isFe) = PAir.*dUe(~isFe);
    Pd(~isFe)  = PAir;
    r = Ared.'*Phi;
    if hasYx, r = r + Yxr*U; end
    if hasIs, r = r - Isr; end
    res = max(abs(r));
    if res < tol, break; end
    J = Ared.'*spdiags(Pd,0,Nb,Nb)*Ared;
    if hasYx, J = J + Yxkk; end
    dx = -(J\r);
    % line-search
    lam = 1; rn0 = norm(r);
    for ls = 1:25
        Ut = U; Ut(kk) = U(kk) + lam*dx;
        dq = Ainc*Ut + E;
        Ph = zeros(Nb,1);
        Ph(isFe) = BH.Bof(dq(isFe)./lFe).*AFe;
        Ph(~isFe) = PAir.*dq(~isFe);
        rt = Ared.'*Ph;
        if hasYx, rt = rt + Yxr*Ut; end
        if hasIs, rt = rt - Isr; end
        if norm(rt) < rn0 || lam < 1e-5, U = Ut; break; end
        lam = lam/2;
    end
end

dUe = Ainc*U + E;
Bfe = BH.Bof(dUe(isFe)./lFe);
S.U = U;
S.B = zeros(Nb,1);  S.B(isFe) = Bfe;
S.Phi = zeros(Nb,1);
S.Phi(isFe) = Bfe.*AFe;
S.Phi(~isFe) = PAir.*dUe(~isFe);
S.Wco = sum(BH.wco(Bfe).*(lFe.*AFe)) + 0.5*sum(S.Phi(~isFe).*dUe(~isFe));
if hasYx, S.Wco = S.Wco + 0.5*(U.'*(Yx*U)); end
S.iter = it; S.res = res;
S.converged = res < max(tol*1e3,1e-6);
end

function v = getdef(s,f,d)
if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end
end
