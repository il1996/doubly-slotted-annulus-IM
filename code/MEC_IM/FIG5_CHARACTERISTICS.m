function FIG5_CHARACTERISTICS
%FIG5_CHARACTERISTICS  Figure 5 : caracteristiques couple-glissement et
%                      courant-glissement de la machine 18,5 kW.
%
%   (a) Couple electromagnetique, (b) courant statorique efficace, de
%   s = 0 a s = 1, machine asynchrone 18,5 kW (48 encoches / 44 barres).
%
%   Quatre jeux superposes :
%     - balayage parametrique EF, 129 glissements resolus a 0,75 s chacun ;
%       la bande s <= 0,125, ou la rugosite point a point reste sous 1,5 %,
%       est distinguee des bandes bruitees (jusqu'a 23 % d'un point au
%       suivant), tracees en gris ;
%     - transitoires EF dedies de 2 s : s = 0,01883 (121,63 N.m ; 19,72 A),
%       rotor bloque s = 1 (104,31 N.m ; 104,4 A) et, pour le courant, le
%       point a vide (8,499 A) ;
%     - reseau avec l'operateur condense, condition Phi_O = 0, (17, 4) ;
%     - reseau avec l'operateur couple aux cavites d'encoche, (33, 16).
%   Vrillage harmonique neutralise dans les deux chaines reseau.
%
%   Donnees : fig5_data.mat
%     sweep_T, sweep_I   balayage EF : colonne 1 glissement, colonne 2
%                        couple (N.m) ou courant (A rms)
%                        -- exportes de Torque Plot 2.tab et
%                           Winding Plot 1.tab (reference ANSYS)
%     sl_op,  T_op,  I_op    reseau Phi_O = 0     -- B10_b1_skewoff.mat
%     sl_cav, T_cav, I_cav   reseau + cavites     -- Z5_sweep_cav.mat
%   Aucune valeur n'est recalculee ici.
%
%   Usage :  FIG5_CHARACTERISTICS
%
%   Auteurs : I. Laouar, A. Boukadoum, N. Mezhoud.

%% ---------------------------------------------------------------- options
USE_LATEX = true;      % false : interprete TeX de MATLAB, compatible Octave

%  Le fichier de donnees est cherche a cote du script, quel que soit le
%  dossier courant de MATLAB.
DATAFILE = fullfile(fileparts(mfilename('fullpath')),'fig5_data.mat');
if exist(DATAFILE,'file') ~= 2
    error('FIG5:donnees', ['fig5_data.mat est introuvable dans %s\n' ...
           'Placez-le a cote de FIG5_CHARACTERISTICS.m.'], ...
          fileparts(mfilename('fullpath')));
end

S_CLEAN = 0.125;       % limite de la bande lisse du balayage parametrique
S_NOM   = 0.01883;     % glissement nominal des transitoires dedies

%% ---------------------------------------------------------------- donnees
D = load(DATAFILE);
sT = D.sweep_T(:,1);  Tsw = D.sweep_T(:,2);      % balayage EF, couple
sI = D.sweep_I(:,1);  Isw = D.sweep_I(:,2);      % balayage EF, courant
clean = sT <= S_CLEAN;

%  transitoires EF dedies
sT_tr = [S_NOM 1.0];        T_tr = [121.63 104.31];
sI_tr = [0 S_NOM 1.0];      I_tr = [8.499 19.72 104.4];

%% ---------------------------------------------------------------- couleurs
cGREY = [149 165 166]/255;     % gris  : bandes bruitees du balayage
cFE   = [ 27  79 114]/255;     % bleu  : elements finis
cNEU  = [192  57  43]/255;     % rouge : reseau, Phi_O = 0
cCAV  = [ 30 132  73]/255;     % vert  : reseau + cavites

if USE_LATEX
    itp = 'latex';
    sxl = 'slip $s$';
    sla = {'FE parametric sweep, noisy bands', ...
           'FE parametric sweep, $s \le 0.125$', ...
           'FE dedicated 2 s transients', ...
           'network, operator $\Phi_O=0$ (17, 4)', ...
           'network, operator + cavities (33, 16)'};
    slb = {'FE parametric sweep','FE dedicated transients', ...
           'network, $\Phi_O=0$ (17, 4)','network, cavities (33, 16)'};
else
    itp = 'tex';
    sxl = 'slip s';
    sla = {'FE parametric sweep, noisy bands', ...
           'FE parametric sweep, s \leq 0.125', ...
           'FE dedicated 2 s transients', ...
           'network, operator \Phi_O = 0 (17, 4)', ...
           'network, operator + cavities (33, 16)'};
    slb = {'FE parametric sweep','FE dedicated transients', ...
           'network, \Phi_O = 0 (17, 4)','network, cavities (33, 16)'};
end

%% ---------------------------------------------------------------- figure
f = figure('Units','inches','Position',[1 1 6.8 2.7], ...
           'Color','w','PaperPositionMode','auto');

% ---------------------------------------------------------- (a) couple
axA = subplot(1,2,1); hold(axA,'on');
set(axA,'Box','on','TickDir','in','LineWidth',0.6, ...
        'FontName','Helvetica','FontSize',8.5);

ha(1) = plot(axA,sT(~clean),Tsw(~clean),'.','Color',cGREY,'MarkerSize',6);
ha(2) = plot(axA,sT( clean),Tsw( clean),'.','Color',cFE  ,'MarkerSize',8);
ha(3) = plot(axA,sT_tr,T_tr,'o','LineStyle','none','MarkerSize',6, ...
             'MarkerEdgeColor',cFE,'MarkerFaceColor','none','LineWidth',1.2);
ha(4) = plot(axA,D.sl_op ,D.T_op ,'-' ,'Color',cNEU,'LineWidth',1.2);
ha(5) = plot(axA,D.sl_cav,D.T_cav,'--','Color',cCAV,'LineWidth',1.2);

xlim(axA,[0 1]); ylim(axA,[0 400]);
xlabel(axA,sxl,'Interpreter',itp,'FontSize',8.5);
ylabel(axA,'electromagnetic torque (N m)','FontSize',8.5);
title(axA,'(a)','FontSize',8.5);
lgA = legend(axA,ha,sla,'Interpreter',itp,'Location','northeast','FontSize',6.2);
set(lgA,'Box','off');

% ---------------------------------------------------------- (b) courant
axB = subplot(1,2,2); hold(axB,'on');
set(axB,'Box','on','TickDir','in','LineWidth',0.6, ...
        'FontName','Helvetica','FontSize',8.5);

hb(1) = plot(axB,sI,Isw,'.','Color',cFE,'MarkerSize',6);
hb(2) = plot(axB,sI_tr,I_tr,'o','LineStyle','none','MarkerSize',6, ...
             'MarkerEdgeColor',cFE,'MarkerFaceColor','none','LineWidth',1.2);
hb(3) = plot(axB,D.sl_op ,D.I_op ,'-' ,'Color',cNEU,'LineWidth',1.2);
hb(4) = plot(axB,D.sl_cav,D.I_cav,'--','Color',cCAV,'LineWidth',1.2);

xlim(axB,[0 1]); ylim(axB,[0 120]);
xlabel(axB,sxl,'Interpreter',itp,'FontSize',8.5);
ylabel(axB,'stator current (A rms)','FontSize',8.5);
title(axB,'(b)','FontSize',8.5);
lgB = legend(axB,hb,slb,'Interpreter',itp,'Location','southeast','FontSize',6.2);
set(lgB,'Box','off');

%% ---------------------------------------------------------------- export
if exist('exportgraphics','file')                % MATLAB R2020a et au-dela
    exportgraphics(f,'fig5_characteristics.png','Resolution',400);
    exportgraphics(f,'fig5_characteristics.pdf','ContentType','vector');
else
    print(f,'-dpng','-r400','fig5_characteristics.png');
    print(f,'-dpdf','-painters','fig5_characteristics.pdf');
end
end
