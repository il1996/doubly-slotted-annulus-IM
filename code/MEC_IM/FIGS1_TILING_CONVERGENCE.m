function FIGS1_TILING_CONVERGENCE(outdir)
%FIGS1_TILING_CONVERGENCE  Figure S1 du supplement : convergence du rapport d'encochage en
%                         nombre de colonnes d'ouverture.
%
%   k_C a phi = 0, fer infiniment permeable, troncature harmonique
%   N_h = 8192, machine asynchrone 18,5 kW (48 encoches / 44 barres).
%
%   (a) operateur avec la condition Phi_O = 0 : chapeau symetrique (P1),
%       chapeau asymetrique (P1a, eq. 6) et pavage gradue vers les coins
%       en chapeau asymetrique.
%   (b) operateur couple aux cavites d'encoche (eq. 7) : chapeau
%       symetrique et pavage gradue en chapeau asymetrique.
%
%   Provenance des series (aucune valeur n'est recalculee ici) :
%     chapeau symetrique, Phi_O = 0   tiling_sweep_inf_iron.json, 'p1_<nT>_<nO>'
%     chapeau asymetrique, Phi_O = 0  tiling_sweep_inf_iron.json, 'p1a_<nT>_<nO>'
%     pavage gradue, Phi_O = 0        t4_graded_results.json, 'neumann_graded'
%     cavites, chapeau symetrique     prod_op_results.json, 'cav_sweep'
%     pavage gradue + cavites         t4_graded_results.json, 'cav_graded'
%
%   Reperes horizontaux : limite EF de la condition Phi_O = 0 (1,460),
%   reference EF sur encoches reelles (1,2934), Carter forme conforme
%   exacte (1,267932).
%
%   Les points (n_T, n_O) = (65, 32) n'ont pas ete calcules : les deux
%   courbes concernees s'arretent a n_O = 16.
%
%   Sorties : fig3_tiling_convergence.fig (figure MATLAB rouvrable),
%   .pdf (vectoriel) et .png (400 ppp). Le .fig ne s'insere pas dans Word :
%   le supplement prend le .png.
%
%   Renommee le 20 septembre 2026 : cette figure etait la Fig. 3 du corps ;
%   le balayage b0/g l'a remplacee (FIG3_B0G_SWEEP) et elle est passee au
%   supplement sous le numero S1. Les series et leur provenance sont
%   inchangees.
%
%   Usage :  FIGS1_TILING_CONVERGENCE
%            FIGS1_TILING_CONVERGENCE('chemin/de/sortie')
%   (ancien nom : FIG3_TILING_CONVERGENCE)
%
%   Auteurs : I. Laouar, A. Boukadoum, N. Mezhoud.

%% ---------------------------------------------------------------- options
USE_LATEX = true;      % false : interprete TeX de MATLAB, compatible Octave

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(outdir)
    figdir = fullfile(here, '..', 'article', 'figures');      % depuis code/MEC_IM/
    if exist(figdir, 'dir') == 7, outdir = figdir; else, outdir = here; end
end
if exist(outdir, 'dir') ~= 7, mkdir(outdir); end

%% ---------------------------------------------------------------- donnees
nO_a = [2 4 8 16];             % colonnes d'ouverture, panneau (a)
nO_b = [4 8 16 32];            % colonnes d'ouverture, panneau (b)

% (a) condition Phi_O = 0, chapeau symetrique, n_T = 9 / 17 / 33 / 65
sym_a = [1.696410 1.622713 1.578979 1.554411 ;   % n_T = 9
         1.680080 1.591247 1.549597 1.527523 ;   % n_T = 17
         1.675591 1.581470 1.529845 1.509537 ;   % n_T = 33
         1.674396 1.579112 1.523658 1.496652];   % n_T = 65
% (a) chapeau asymetrique, n_T = 33 / 65
asym_a = [1.471297 1.518358 1.538134 1.547119 ;  % n_T = 33
          1.433997 1.476218 1.493466 1.501132];  % n_T = 65
% (a) pavage gradue en chapeau asymetrique, n_T = 33 / 65
grad_a = [1.397840 1.441333 1.456446 1.465442 ;  % n_T = 33
          1.397124 1.440495 1.455546 1.464902];  % n_T = 65

% (b) cavites d'encoche, chapeau symetrique, n_T = 17 / 33 / 65
sym_b = [1.363139 1.325505 1.309222 1.302386 ;   % n_T = 17
         1.358403 1.317462 1.302394 1.296216 ;   % n_T = 33
         1.357249 1.314950 1.297933      NaN];   % n_T = 65 : (65,32) absent
% (b) pavage gradue + cavites, n_T = 33 / 65
grad_b = [1.278627 1.283510 1.284861 1.290004 ;  % n_T = 33
          1.278203 1.283201 1.284720      NaN];  % n_T = 65 : (65,32) absent

KC_FE_REAL   = 1.2934;         % EF, encoches reelles, phi = 0
KC_FE_NEU    = 1.460;          % limite EF de la condition Phi_O = 0
KC_CARTER_EX = 1.267932;       % Carter, forme conforme exacte

%% ---------------------------------------------------------------- couleurs
cSym_a  = [245 183 177 ; 236 112  99 ; 203  67  53 ; 120  40  31]/255;
cAsym_a = [123  36  28]/255;
cSym_b  = [125 206 160 ; 34 153  84 ; 20  90  50]/255;
cGrad   = [125  60 152]/255;
cFE     = [ 27  79 114]/255;
cNeu    = [192  57  43]/255;
cCarter = [ 93 109 126]/255;

mk_a = {'o','s','^','d'};      % marqueurs par n_T, panneau (a)
mk_b = {'s','^','d'};          % marqueurs par n_T, panneau (b)
nT_sym_a  = [9 17 33 65];      % n_T des series du panneau (a)
nT_pair   = [33 65];           % n_T des series a deux courbes
nT_sym_b  = [17 33 65];        % n_T des series du panneau (b)

if USE_LATEX
    itp = 'latex';  nT = @(v) sprintf('$n_T=%d$',v);
    sxl = 'opening columns $n_O$';
    syl = '$k_C$ at $\varphi = 0$, infinite iron';
    sta = '(a) operator with $\Phi_O = 0$';
else
    itp = 'tex';    nT = @(v) sprintf('n_T = %d',v);
    sxl = 'opening columns n_O';
    syl = 'k_C at \varphi = 0, infinite iron';
    sta = '(a) operator with \Phi_O = 0';
end

%% ---------------------------------------------------------------- figure
f = figure('Units','inches','Position',[1 1 6.8 2.8], ...
           'Color','w','PaperPositionMode','auto');

% ---------------------------------------------------------- panneau (a)
axA = subplot(1,2,1); hold(axA,'on');
set(axA,'Box','on','TickDir','in','LineWidth',0.6,'XScale','log', ...
        'XTick',nO_a,'XTickLabel',arrayfun(@num2str,nO_a,'UniformOutput',false), ...
        'FontName','Helvetica','FontSize',8.5);

plot(axA,nO_a([1 end]),KC_FE_NEU *[1 1],'-.','Color',cNeu,'LineWidth',1,'HandleVisibility','off');
plot(axA,nO_a([1 end]),KC_FE_REAL*[1 1],'-.','Color',cFE ,'LineWidth',1,'HandleVisibility','off');

ha = []; la = {}; k = 0;
for i = 1:4                                       % chapeau symetrique
    k = k+1;
    ha(k) = plot(axA,nO_a,sym_a(i,:),'-','Color',cSym_a(i,:),'LineWidth',1, ...
                 'Marker',mk_a{i},'MarkerSize',3.5, ...
                 'MarkerFaceColor',cSym_a(i,:),'MarkerEdgeColor',cSym_a(i,:));
    la{k} = [nT(nT_sym_a(i)) ', symmetric hat'];
end
for i = 1:2                                       % chapeau asymetrique
    k = k+1;
    ha(k) = plot(axA,nO_a,asym_a(i,:),'--','Color',cAsym_a,'LineWidth',1, ...
                 'Marker',mk_a{i+2},'MarkerSize',3.5, ...
                 'MarkerFaceColor','none','MarkerEdgeColor',cAsym_a);
    la{k} = [nT(nT_pair(i)) ', asymmetric hat'];
end
for i = 1:2                                       % pavage gradue
    k = k+1;
    ha(k) = plot(axA,nO_a,grad_a(i,:),':','Color',cGrad,'LineWidth',1, ...
                 'Marker',mk_a{i+2},'MarkerSize',3.5, ...
                 'MarkerFaceColor','none','MarkerEdgeColor',cGrad);
    la{k} = [nT(nT_pair(i)) ', graded tiling, asym. hat'];
end

text(axA,2.1,KC_FE_NEU +0.008,'FE limit, flux-barrier openings (1.46)', ...
     'FontSize',7,'Color',cNeu,'VerticalAlignment','bottom');
text(axA,2.1,KC_FE_REAL+0.007,'FE, real slots (1.293)', ...
     'FontSize',7,'Color',cFE ,'VerticalAlignment','bottom');

xlim(axA,[nO_a(1) nO_a(end)]); ylim(axA,[1.25 2.05]);
xlabel(axA,sxl,'Interpreter',itp,'FontSize',8.5);
ylabel(axA,syl,'Interpreter',itp,'FontSize',8.5);
title(axA,sta,'Interpreter',itp,'FontSize',8.5);
lgA = legend(axA,ha,la,'Interpreter',itp,'Location','north','FontSize',5.8);
set(lgA,'Box','off');
try, set(lgA,'NumColumns',2); end                 %#ok<TRYNC> R2017b et au-dela

% ---------------------------------------------------------- panneau (b)
axB = subplot(1,2,2); hold(axB,'on');
set(axB,'Box','on','TickDir','in','LineWidth',0.6,'XScale','log', ...
        'XTick',nO_b,'XTickLabel',arrayfun(@num2str,nO_b,'UniformOutput',false), ...
        'FontName','Helvetica','FontSize',8.5);

plot(axB,nO_b([1 end]),KC_FE_REAL  *[1 1],'-.','Color',cFE    ,'LineWidth',1,'HandleVisibility','off');
plot(axB,nO_b([1 end]),KC_CARTER_EX*[1 1],':' ,'Color',cCarter,'LineWidth',1,'HandleVisibility','off');

hb = []; lb = {}; k = 0;
for i = 1:3                                       % chapeau symetrique
    k = k+1;
    hb(k) = plot(axB,nO_b,sym_b(i,:),'-','Color',cSym_b(i,:),'LineWidth',1, ...
                 'Marker',mk_b{i},'MarkerSize',3.5, ...
                 'MarkerFaceColor',cSym_b(i,:),'MarkerEdgeColor',cSym_b(i,:));
    lb{k} = [nT(nT_sym_b(i)) ', symmetric hat'];
end
for i = 1:2                                       % pavage gradue
    k = k+1;
    hb(k) = plot(axB,nO_b,grad_b(i,:),':','Color',cGrad,'LineWidth',1, ...
                 'Marker',mk_b{i+1},'MarkerSize',3.5, ...
                 'MarkerFaceColor','none','MarkerEdgeColor',cGrad);
    lb{k} = [nT(nT_pair(i)) ', graded tiling, asym. hat'];
end

text(axB,4.1,KC_FE_REAL  +0.003,'FE, real slots (1.293)', ...
     'FontSize',7,'Color',cFE    ,'VerticalAlignment','bottom');
text(axB,4.1,KC_CARTER_EX-0.010,'Carter 1.268', ...
     'FontSize',7,'Color',cCarter,'VerticalAlignment','bottom');

xlim(axB,[nO_b(1) nO_b(end)]); ylim(axB,[1.25 1.40]);
xlabel(axB,sxl,'Interpreter',itp,'FontSize',8.5);
ylabel(axB,syl,'Interpreter',itp,'FontSize',8.5);
title(axB,'(b) operator with slot cavities','Interpreter',itp,'FontSize',8.5);
lgB = legend(axB,hb,lb,'Interpreter',itp,'Location','northeast','FontSize',5.8);
set(lgB,'Box','off');

%% ---------------------------------------------------------------- export
% Le nom du fichier reste 'fig3_tiling_convergence' : c'est celui que portent
% l'image publiee, le supplement, README §3 et MANIFEST §9. Seule la figure a
% change de numero (Fig. 3 du corps -> Fig. S1 du supplement), pas le fichier.
base = fullfile(outdir, 'fig3_tiling_convergence');

% .fig : la figure elle-meme, rouvrable et modifiable (openfig ou double-clic)
if exist('savefig','file') == 2
    savefig(f, [base '.fig']);
else
    hgsave(f, [base '.fig']);                    % Octave
end

if exist('exportgraphics','file')                % MATLAB R2020a et au-dela
    exportgraphics(f, [base '.png'], 'Resolution', 400);
    exportgraphics(f, [base '.pdf'], 'ContentType', 'vector');
else
    print(f, '-dpng', '-r400', [base '.png']);
    print(f, '-dpdf', '-painters', [base '.pdf']);
end
fprintf('ecrit : %s.fig, %s.pdf et %s.png\n', base, base, base);
end
