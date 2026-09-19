function FIG2_KC_VS_POSITION
%FIG2_KC_VS_POSITION  Figure 2 : rapport d'encochage contre la position rotor.
%
%   Rapport d'encochage k_C = B_g1^smooth / B_g1^slotted de la machine
%   asynchrone 18,5 kW (48 encoches / 44 barres), a permeabilite du fer
%   infinie, sur un pas d'encoche rotorique tau_r = 8,18 deg mecaniques.
%   Douze positions calculees, la treizieme etant la premiere par
%   periodicite.
%
%   Provenance des cinq series (aucune valeur n'est recalculee ici) :
%     EF encoches reelles      prod_fem_results.json, champ 'pos_real'
%     EF ouvertures Neumann    prod_fem_results.json, champ 'pos_neu'
%     operateur + cavites      prod_op_results.json,  champ 'pos_cav'
%     operateur Phi_O = 0      prod_op_results.json,  champs 'pos_op_33_16'
%                                                     et 'pos_op_17_4'
%   Troncature harmonique N_h = 8192 pour les trois series d'operateur.
%
%   Carter, forme conforme exacte : 1,267932 (forme approchee : 1,266466).
%
%   Usage :  FIG2_KC_VS_POSITION            % affiche et exporte la figure
%
%   Auteurs : I. Laouar, A. Boukadoum, N. Mezhoud.

%% ---------------------------------------------------------------- options
USE_LATEX = true;      % true : rendu LaTeX (comme dans l'article)
                       % false : interprete TeX de MATLAB, compatible Octave

%% ---------------------------------------------------------------- donnees
phi = [0 1 2 3 4 5 6 7 8 9 10 11]/12;      % position rotor, en pas d'encoche

kC.fe     = [1.294435 1.283279 1.268469 1.253700 1.243073 1.239087 ...
             1.243072 1.253699 1.268413 1.283275 1.294455 1.298614];
kC.cav    = [1.302394 1.290566 1.274851 1.259303 1.248102 1.243921 ...
             1.248103 1.259305 1.274854 1.290568 1.302396 1.306829];
kC.feNeu  = [1.457722 1.435301 1.406103 1.377935 1.358062 1.350732 ...
             1.358036 1.377993 1.406058 1.435275 1.457704 1.466232];
kC.op33   = [1.509537 1.482347 1.447185 1.413862 1.390410 1.381901 ...
             1.390410 1.413862 1.447185 1.482347 1.509537 1.519972];
kC.op17   = [1.591247 1.555767 1.510347 1.468322 1.438924 1.428536 ...
             1.438924 1.468322 1.510347 1.555767 1.591247 1.605035];

KC_CARTER_EXACT = 1.267932;

%% ---------------------------------------------------------------- couleurs
col.fe     = [ 27  79 114]/255;    % bleu   : reference EF
col.cav    = [ 30 132  73]/255;    % vert   : operateur + cavites
col.neu    = [192  57  43]/255;    % rouge  : EF, ouvertures en barriere
col.op     = [230 126  34]/255;    % orange : operateur Phi_O = 0 (33,16)
col.op2    = [240 178 122]/255;    % orange clair : idem (17,4)
col.carter = [ 93 109 126]/255;    % gris   : Carter

%% ---------------------------------------------------------------- figure
f = figure('Units','inches','Position',[1 1 4.6 3.0], ...
           'Color','w','PaperPositionMode','auto');
ax = axes('Parent',f,'Box','on','TickDir','in','LineWidth',0.6, ...
          'FontName','Helvetica','FontSize',8.5,'NextPlot','add');

plot(ax,[0 1],KC_CARTER_EXACT*[1 1],'-.','Color',col.carter,'LineWidth',1.0, ...
     'HandleVisibility','off');

h(1) = closed(ax,phi,kC.fe   ,col.fe , '-','o',3.5,1.2);
h(2) = closed(ax,phi,kC.cav  ,col.cav, '-','s',3.0,1.0);
h(3) = closed(ax,phi,kC.feNeu,col.neu,'--','^',3.5,1.2);
h(4) = closed(ax,phi,kC.op33 ,col.op ,'--','v',3.0,1.0);
h(5) = closed(ax,phi,kC.op17 ,col.op2,':' ,'d',2.5,0.9);

text(ax,0.02,KC_CARTER_EXACT+0.006, ...
     'Carter (exact conformal form) 1.2679', ...
     'FontSize',7,'Color',col.carter,'VerticalAlignment','bottom');

xlim(ax,[0 1]); ylim(ax,[1.20 1.65]);
if USE_LATEX
    itp = 'latex';
    sx  = 'rotor position $\varphi/\tau_r$ (one rotor slot pitch)';
    sy  = ['slotting ratio $k_C = B_{g1}^{\mathrm{smooth}}/' ...
           'B_{g1}^{\mathrm{slotted}}$'];
    sl  = {'FE, real slots (reference)', ...
           'operator + slot cavities (33, 16)', ...
           'FE, openings as flux barriers', ...
           'operator, $\Phi_O=0$ (33, 16)', ...
           'operator, $\Phi_O=0$ (17, 4)'};
else
    itp = 'tex';
    sx  = 'rotor position \varphi/\tau_r (one rotor slot pitch)';
    sy  = 'slotting ratio k_C = B_{g1}^{smooth}/B_{g1}^{slotted}';
    sl  = {'FE, real slots (reference)', ...
           'operator + slot cavities (33, 16)', ...
           'FE, openings as flux barriers', ...
           'operator, \Phi_O = 0 (33, 16)', ...
           'operator, \Phi_O = 0 (17, 4)'};
end
xlabel(ax,sx,'Interpreter',itp,'FontSize',8.5);
ylabel(ax,sy,'Interpreter',itp,'FontSize',8.5);

lg = legend(h,sl,'Interpreter',itp,'Location','north','FontSize',6.8);
set(lg,'Box','off');

%% ---------------------------------------------------------------- export
if exist('exportgraphics','file')                % MATLAB R2020a et au-dela
    exportgraphics(f,'fig2_kc_vs_position.png','Resolution',400);
    exportgraphics(f,'fig2_kc_vs_position.pdf','ContentType','vector');
else                                            % versions anterieures
    print(f,'-dpng','-r400','fig2_kc_vs_position.png');
    print(f,'-dpdf','-painters','fig2_kc_vs_position.pdf');
end
end

% ======================================================================
function hh = closed(ax,x,y,c,ls,mk,ms,lw)
%CLOSED  Trace une serie periodique en refermant le pas sur phi/tau_r = 1.
x = [x 1]; y = [y y(1)];
hh = plot(ax,x,y,'LineStyle',ls,'Color',c,'LineWidth',lw, ...
          'Marker',mk,'MarkerSize',ms, ...
          'MarkerFaceColor',c,'MarkerEdgeColor',c);
end
