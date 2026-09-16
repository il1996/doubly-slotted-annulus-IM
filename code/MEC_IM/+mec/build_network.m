function net = build_network(M, G)
%BUILD_NETWORK  Topologie fixe du réseau de réluctances (fer).
%
%   net = mec.build_network(M,G) construit la partie INVARIANTE du réseau
%   (indépendante de la position rotorique et de l'excitation) : nœuds et
%   branches de fer (dents et culasses stator/rotor). Les branches
%   d'entrefer, variables avec la position, sont ajoutées à la résolution
%   (mec.solve_network).
%
%   Numérotation des nœuds (potentiel scalaire magnétique) :
%       [1 .. Ns]                 culasse stator (un nœud derrière chaque dent)
%       [Ns+1 .. 2Ns]             becs de dents stator (face entrefer)
%       [2Ns+1 .. 2Ns+Nr]         becs de dents rotor  (face entrefer)
%       [2Ns+Nr+1 .. 2Ns+2Nr]     culasse rotor
%
%   Branches de fer (non linéaires, relation H(B)*l) :
%       * Ns segments d'anneau de culasse stator  (Y_i - Y_{i+1}, périodique)
%       * Ns dents stator                          (Y_i - Ts_i), source FMM
%       * Nr dents rotor                           (Tr_j - Yr_j), source FMM
%       * Nr segments d'anneau de culasse rotor    (Yr_j - Yr_{j+1}, périod.)
%
%   Chaque branche porte : nœuds (a,b), longueur l, section A, drapeau
%   'kind' de région (pour la carte de B et les pertes), et un indice de
%   source FMM (isrc_s pour dents stator, isrc_r pour dents rotor).
%
%   Voir aussi : mec.solve_network, mec.geometry.

Ns = M.Ns;  Nr = M.Nr;

% Indices de nœuds
YS = @(i) i;                  % culasse stator
TS = @(i) Ns + i;             % bec dent stator
TR = @(j) 2*Ns + j;           % bec dent rotor
YR = @(j) 2*Ns + Nr + j;      % culasse rotor
net.Nnodes = 2*Ns + 2*Nr;
net.idx.YS = YS; net.idx.TS = TS; net.idx.TR = TR; net.idx.YR = YR;
net.Ns = Ns; net.Nr = Nr;

a = []; b = []; l = []; A = []; kind = {}; srcS = []; srcR = [];

wS = @(i) mod(i-1,Ns)+1;
wR = @(j) mod(j-1,Nr)+1;

% --- Anneau de culasse stator ---
for i = 1:Ns
    a(end+1)=YS(i); b(end+1)=YS(wS(i+1)); %#ok<AGROW>
    l(end+1)=G.l_ys; A(end+1)=G.A_ys; kind{end+1}='ys'; %#ok<AGROW>
    srcS(end+1)=0; srcR(end+1)=0; %#ok<AGROW>
end
% --- Dents stator (avec source FMM de bobinage) ---
for i = 1:Ns
    a(end+1)=YS(i); b(end+1)=TS(i); %#ok<AGROW>
    l(end+1)=G.l_ts; A(end+1)=G.A_ts; kind{end+1}='ts'; %#ok<AGROW>
    srcS(end+1)=i; srcR(end+1)=0; %#ok<AGROW>
end
% --- Dents rotor (avec source FMM de barre) ---
for j = 1:Nr
    a(end+1)=TR(j); b(end+1)=YR(j); %#ok<AGROW>
    l(end+1)=G.l_tr; A(end+1)=G.A_tr; kind{end+1}='tr'; %#ok<AGROW>
    srcS(end+1)=0; srcR(end+1)=j; %#ok<AGROW>
end
% --- Anneau de culasse rotor ---
for j = 1:Nr
    a(end+1)=YR(j); b(end+1)=YR(wR(j+1)); %#ok<AGROW>
    l(end+1)=G.l_yr; A(end+1)=G.A_yr; kind{end+1}='yr'; %#ok<AGROW>
    srcS(end+1)=0; srcR(end+1)=0; %#ok<AGROW>
end

net.iron.a = a(:); net.iron.b = b(:);
net.iron.l = l(:); net.iron.A = A(:);
net.iron.kind = kind(:);
net.iron.srcS = srcS(:); net.iron.srcR = srcR(:);
net.iron.n = numel(a);

% Nœud de référence (potentiel imposé nul) : une culasse rotor
net.ref = YR(1);

end
