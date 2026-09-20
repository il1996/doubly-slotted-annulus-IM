function FIG3_B0G_SWEEP(jsonfile, outdir)
%FIG3_B0G_SWEEP  Figure 3 : rapport d'encochage contre le rapport
%                ouverture / entrefer, l'entrefer seul variant.
%
%   Ecart a la reference par elements finis (formulation A, moyenne sur
%   douze positions rotor, recalculee a chaque entrefer) de :
%     - le produit de Carter, forme conforme exacte ;
%     - l'operateur avec la condition Phi_O = 0 ;
%     - l'operateur couple aux cavites d'encoche ;
%   barres : etendue min-max de la reference sur la position.
%
%   Pavage (33, 16), chapeau symetrique, N_h = 8192, b_0 = 2,000 mm aux
%   deux surfaces, R_s fixe, R_r = R_s - g.
%
%   AUCUNE VALEUR N'EST ECRITE DANS CE FICHIER. Tout est relu de
%     outputs/python/prod_sweep_b0g_results.json
%   produit par code/python/prod_sweep_b0g.py (gardes G1-G5), la ligne
%   b0/g = 12 promue au maillage divise par prod_sweep_b0g_meshrow.py.
%   Le champ 'fe' de chaque entrefer porte le maillage retenu ; a
%   b0/g = 12 et 16 le maillage de production est conserve a cote sous
%   'fe_production_lc' et n'est pas trace.
%
%   Le script imprime les nombres traces, a comparer ligne a ligne avec
%   les Tables 1 et 2 de outputs/python/prod_sweep_b0g_out.txt.
%
%   Le script cherche le transcript, dans l'ordre :
%     1. le chemin donne en premier argument ;
%     2. ../../outputs/python/ (place dans code/MEC_IM/ de l'archive) ;
%     3. le dossier du script lui-meme ;
%     4. outputs/python/ dans chaque dossier parent, jusqu'a la racine.
%
%   Sorties dans le dossier de figures : fig3_b0g_sweep.fig (figure MATLAB
%   rouvrable), .pdf (vectoriel) et .png (1350 x 800 px, soit 3,375 x 2,000
%   pouces a 400 ppp, la largeur d'une colonne de l'article).
%
%   Usage :  FIG3_B0G_SWEEP                       % depuis code/MEC_IM/
%            FIG3_B0G_SWEEP(cheminDuJson)
%            FIG3_B0G_SWEEP(cheminDuJson, dossierDeSortie)
%
%   Auteurs : I. Laouar, A. Boukadoum, N. Mezhoud, S. Lekhchine.

%% ---------------------------------------------------------------- options
USE_LATEX = true;      % false : interprete TeX de MATLAB, compatible Octave

here = fileparts(mfilename('fullpath'));
JSONNAME = 'prod_sweep_b0g_results.json';

% ---- 1 a 4 : ou chercher le transcript -----------------------------------
cand = {};
if nargin >= 1 && ~isempty(jsonfile)
    if exist(jsonfile, 'dir') == 7
        cand{end+1} = fullfile(jsonfile, JSONNAME);              % un dossier
    else
        cand{end+1} = jsonfile;                                  % un fichier
    end
end
cand{end+1} = fullfile(here, '..', '..', 'outputs', 'python', JSONNAME);
cand{end+1} = fullfile(here, JSONNAME);
up = here;
for k = 1:6                                   % remonte jusqu'a la racine
    parent = fileparts(up);
    if isempty(parent) || strcmp(parent, up), break; end
    up = parent;
    cand{end+1} = fullfile(up, 'outputs', 'python', JSONNAME);   %#ok<AGROW>
end

jsonfile = '';
for k = 1:numel(cand)
    if exist(cand{k}, 'file') == 2, jsonfile = cand{k}; break; end
end
if isempty(jsonfile)
    msg = sprintf(['transcript %s introuvable. Cherche dans :\n'], JSONNAME);
    for k = 1:numel(cand), msg = [msg sprintf('   %s\n', cand{k})]; end %#ok<AGROW>
    msg = [msg sprintf(['\nDonnez son chemin en argument, par exemple :\n' ...
           '   FIG3_B0G_SWEEP(''...\\archive_v1.1.0\\outputs\\python\\%s'')\n'], JSONNAME)];
    error('FIG3:json', '%s', msg);
end

% ---- dossier de sortie ---------------------------------------------------
if nargin < 2 || isempty(outdir)
    figdir = fullfile(here, '..', 'article', 'figures');         % depuis code/MEC_IM/
    if exist(figdir, 'dir') == 7, outdir = figdir; else, outdir = here; end
end
if exist(outdir, 'dir') ~= 7, mkdir(outdir); end

%% ---------------------------------------------------------------- donnees
S = jsondecode(fileread(jsonfile));

% le point nominal et les entrefers du balayage, dans un ordre quelconque :
% jsondecode prefixe les cles numeriques ('2' -> 'x2'), on ne s'y fie pas.
names = [{'nominal'} ; fieldnames(S.sweep)];
n = numel(names);
[x, feMean, lo, hi, carterEx, carterAp, neu, cav, lcUsed] = deal(zeros(1, n));

for i = 1:n
    if strcmp(names{i}, 'nominal'), e = S.nominal; else, e = S.sweep.(names{i}); end
    fe = e.fe;                                  % maillage retenu pour cet entrefer
    x(i)        = e.b0_over_g;
    feMean(i)   = fe.mean;
    lo(i)       = (fe.min        / fe.mean - 1) * 100;   % <= 0
    hi(i)       = (fe.max        / fe.mean - 1) * 100;   % >= 0
    carterEx(i) = (e.carter.k_exact / fe.mean - 1) * 100;
    carterAp(i) = (e.carter.k       / fe.mean - 1) * 100;
    neu(i)      = (e.op_neu.mean    / fe.mean - 1) * 100;
    cav(i)      = (e.op_cav.mean    / fe.mean - 1) * 100;
    lcUsed(i)   = fe.lc;
end

[x, ord] = sort(x);
feMean = feMean(ord); lo = lo(ord); hi = hi(ord); lcUsed = lcUsed(ord);
carterEx = carterEx(ord); carterAp = carterAp(ord); neu = neu(ord); cav = cav(ord);
xnom = S.nominal.b0_over_g;

% ---- controle imprime, a confronter au transcript -----------------------
fprintf('\nFIG3_B0G_SWEEP -- relu de %s\n', jsonfile);
fprintf('%7s %8s %10s %10s %10s %9s %9s %9s %8s\n', ...
        'b0/g', 'lc', 'FE mean', 'Carter ap', 'Carter ex', 'Phi_O=0', 'cavites', 'bande -', 'bande +');
for i = 1:n
    fprintf('%7.3f %8.4f %10.6f %+10.3f %+10.3f %+9.3f %+9.3f %+9.3f %+8.3f\n', ...
            x(i), lcUsed(i), feMean(i), carterAp(i), carterEx(i), neu(i), cav(i), lo(i), hi(i));
end
fprintf(['bornes du paragraphe : |Carter ex| <= %.3f %%, |Carter ap| <= %.3f %%, ' ...
         'cavites %+.3f a %+.3f %%, Phi_O=0 %+.3f a %+.3f %%\n\n'], ...
        max(abs(carterEx)), max(abs(carterAp)), min(cav), max(cav), min(neu), max(neu));

%% ---------------------------------------------------------------- couleurs
cOp     = [230 126  34]/255;   % operateur, Phi_O = 0        (#e67e22)
cCav    = [ 30 132  73]/255;   % operateur + cavites         (#1e8449)
cCarter = [ 93 109 126]/255;   % produit de Carter, exact    (#5d6d7e)
cFE     = [ 27  79 114]/255;   % barres de la reference EF   (#1b4f72)
cGrid   = [221 221 221]/255;
cAnnot  = [ 51  51  51]/255;

if USE_LATEX
    itp = 'latex';
    sxl = 'slot opening to gap ratio $b_0/g$ (linear scale)';
    lOp = 'operator, $\Phi_O = 0$, (33, 16)';
else
    itp = 'tex';
    sxl = 'slot opening to gap ratio b_0/g (linear scale)';
    lOp = 'operator, \Phi_O = 0, (33, 16)';
end
syl  = 'deviation from FE mean (%)';
lCav = 'operator + slot cavities, (33, 16)';
lCar = 'Carter''s product, exact conformal form';
lFE  = 'FE reference: mean = 0, bars = min-max over position';

%% ---------------------------------------------------------------- figure
f = figure('Units','inches','Position',[1 1 3.375 2.0], ...
           'Color','w','PaperPositionMode','auto');
ax = axes('Parent',f); hold(ax,'on');
set(ax,'Box','on','TickDir','in','LineWidth',0.6, ...
       'FontName','Helvetica','FontSize',7.5, ...
       'YGrid','on','GridColor',cGrid,'GridAlpha',1,'GridLineStyle','-','Layer','bottom');

% bornes verticales, calculees avant les traces qui s'y appuient
ymin = min([lo, carterEx, cav, neu]) - 1.2;
ymax = max([hi, carterEx, cav, neu]) + 1.2;

% machine de l'article
plot(ax, [xnom xnom], [ymin ymax], ':', 'Color', cAnnot, 'LineWidth', 0.6, ...
     'HandleVisibility','off');

% reference EF : moyenne = 0, barres min-max sur la position
hFE = errorbar(ax, x, zeros(1,n), -lo, hi, 'LineStyle','none', ...
               'Color', cFE, 'LineWidth', 0.8, 'CapSize', 2.5);
plot(ax, [min(x)-1 max(x)+1], [0 0], '-', 'Color', cFE, 'LineWidth', 0.4, ...
     'HandleVisibility','off');

hOp  = plot(ax, x, neu,      '--', 'Color', cOp,     'LineWidth', 0.8, ...
            'Marker','v','MarkerSize',3.5,'MarkerFaceColor',cOp,    'MarkerEdgeColor',cOp);
hCav = plot(ax, x, cav,      '-',  'Color', cCav,    'LineWidth', 0.8, ...
            'Marker','s','MarkerSize',3.5,'MarkerFaceColor',cCav,   'MarkerEdgeColor',cCav);
hCar = plot(ax, x, carterEx, '-.', 'Color', cCarter, 'LineWidth', 0.8, ...
            'Marker','o','MarkerSize',3.5,'MarkerFaceColor',cCarter,'MarkerEdgeColor',cCarter);

xlim(ax, [min(x)-1 max(x)+1]);
ylim(ax, [ymin ymax]);

xt = round(x, 1);                                  % 2  4  8.2  12  16
set(ax, 'XTick', xt, 'XTickLabel', arrayfun(@(v) num2str(v), xt, 'UniformOutput', false));
try
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = 1:(ceil(max(x)) + 1);
catch                                              % Octave : pas de MinorTickValues
    set(ax, 'XMinorTick', 'on');
end

text(ax, xnom + 0.25, ymin + 0.35, 'this machine', ...
     'FontName','Helvetica','FontSize',6.2,'Color',cAnnot, ...
     'HorizontalAlignment','left','VerticalAlignment','bottom');

xlabel(ax, sxl, 'Interpreter', itp, 'FontSize', 8);
ylabel(ax, syl, 'Interpreter', 'tex', 'FontSize', 8);

lg = legend(ax, [hOp hCav hCar hFE], {lOp, lCav, lCar, lFE}, ...
            'Interpreter', itp, 'FontSize', 6.5, 'Box', 'off', ...
            'Location', 'northwest');
set(lg, 'Units', 'normalized');          % on garde la taille calculee,
lgp = get(lg, 'Position');               % on ne deplace que le coin
set(lg, 'Position', [0.17, 0.50, lgp(3), lgp(4)]);

%% ---------------------------------------------------------------- export
base = fullfile(outdir, 'fig3_b0g_sweep');

% .fig : la figure elle-meme, rouvrable et modifiable sous MATLAB
%        (openfig, ou double-clic). Ce format ne s'insere pas dans Word :
%        l'article prend le .png, l'impression le .pdf.
if exist('savefig', 'file') == 2
    savefig(f, [base '.fig']);
else
    hgsave(f, [base '.fig']);                    % Octave
end

if exist('exportgraphics', 'file') == 2 || exist('exportgraphics', 'builtin') == 5
    exportgraphics(f, [base '.pdf'], 'ContentType', 'vector');
    exportgraphics(f, [base '.png'], 'Resolution', 400);   % 1350 x 800 px
else
    set(f, 'PaperUnits','inches', 'PaperSize',[3.375 2.0], ...
           'PaperPosition',[0 0 3.375 2.0]);
    print(f, [base '.pdf'], '-dpdf');
    print(f, [base '.png'], '-dpng', '-r400');
end
fprintf('ecrit : %s.fig, %s.pdf et %s.png\n', base, base, base);
end
