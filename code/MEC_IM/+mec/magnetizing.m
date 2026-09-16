function R = magnetizing(ctx, Im)
%MAGNETIZING  Branche magnétisante saturée par le réseau de réluctances.
%
%   R = mec.magnetizing(ctx, Im) résout le réseau MEC non linéaire pour un
%   courant magnétisant efficace Im [A] (courant à vide, réparti en
%   triphasé équilibré à l'instant de phase A maximale) et renvoie :
%       R.Xm    : réactance magnétisante SATURÉE (sécante) [ohm]
%       R.Lm    : inductance magnétisante saturée [H]
%       R.E1    : f.e.m. induite fondamentale (efficace) [V]
%       R.Bg1   : induction d'entrefer fondamentale (crête) [T]
%       R.S     : solution de champ complète (pour cartes B et pertes fer)
%
%   ctx regroupe les objets pré-calculés : ctx.M,G,W,BH,net,AG.
%   La saturation des dents et culasses est prise en compte exactement
%   (Newton), ce qui rend Xm dépendant de Im (courbe de magnétisation).
%   C'est la contribution centrale du MEC au modèle de la machine :
%   Xm(Im) sature, là où un schéma équivalent classique le suppose constant.
%
%   Voir aussi : mec.solve_network, mec.equivalent_circuit.

M = ctx.M; G = ctx.G; W = ctx.W; BH = ctx.BH; net = ctx.net; AG = ctx.AG;

% Courants triphasés instantanés (phase A maximale)
i3 = sqrt(2)*Im*[1; -0.5; -0.5];
Fs = W.slotMMF(i3);

% Résolution du champ (rotor passif : pas de FMM de barre à vide)
S = mec.solve_network(net, G, BH, AG, Fs, [], M.opt);

% Fondamental spatial de l'induction MOYENNE d'entrefer (p paires de pôles)
Ns = M.Ns; p = M.p;
th = 2*pi*(0:Ns-1)/Ns;
Bg = S.Bgap_avg_i(:).';
Bg1 = abs( (2/Ns)*sum(Bg.*exp(-1j*p*th)) );

% Flux fondamental par pôle et f.e.m. induite
tp   = pi*G.Ds/(2*p);
Phi1 = (2/pi)*Bg1*tp*M.L;
lam1 = W.kw1*W.Nph*Phi1;            % flux totalisé fondamental (crête)
E1   = M.w*lam1/sqrt(2);            % f.e.m. efficace
Xm   = E1/max(Im,eps);
Lm   = Xm/M.w;

%  CONTRÔLE INDÉPENDANT DE Bg1. Ci-dessus, Bg1 vient du FLUX par dent
%  divisé par le pas — une reconstruction, qui suppose l'induction
%  uniforme sur le pas. L'opérateur donne accès au champ PONCTUEL au
%  mi-entrefer : on en extrait le même fondamental, sans cette hypothèse.
%  L'écart entre les deux mesure ce que la moyenne par dent efface.
%
%  MESURE. Il dépend fortement du point de fonctionnement :
%      Im = 0,2 A (référence non saturée) : 0,03245 / 0,02918  -> -10,1 %
%      pleine charge  s = 0,0188           : 0,9093  / 0,8933   ->  -1,8 %
%  Au point de charge les deux estimateurs se rejoignent, et c'est la
%  reconstruction flux/pas qui reste la plus proche de l'EF (-1,2 % contre
%  -2,9 % pour le champ ponctuel, référence 0,920 T). R.Bg1 demeure donc la
%  grandeur rapportée ; R.Bg1_field est un CONTRÔLE, pas un remplaçant.
%  L'écart de 10 % à faible excitation reste à expliquer.
%  Le test porte sur 'expand' : seul l'opérateur CONDENSÉ expose
%  field(Ut,r,thq) à trois arguments. L'opérateur fin de mec.airgap_fourier
%  attend field(Us,Ur,r,thq) — une signature différente, qu'il ne faut pas
%  appeler à l'aveugle.
if isfield(S,'Usurf') && isfield(AG,'expand') && isfield(AG,'field')
    thq = linspace(0, 2*pi, 4*Ns+1); thq(end) = [];
    Br  = AG.field(S.Usurf, 0.5*(G.Rs+G.Rr), thq);
    Br  = Br(:).';
    R.Bg1_field = abs( (2/numel(thq))*sum(Br.*exp(-1j*p*thq)) );
else
    R.Bg1_field = NaN;
end

R.Xm = Xm; R.Lm = Lm; R.E1 = E1; R.Bg1 = Bg1; R.Phi1 = Phi1;
R.S = S; R.Im = Im;
end
