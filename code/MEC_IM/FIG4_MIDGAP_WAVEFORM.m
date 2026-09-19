function FIG4_MIDGAP_WAVEFORM
%FIG4_MIDGAP_WAVEFORM  Figure 4 : induction radiale au mi-entrefer sur deux
%                      pas d'encoche stator, autour du centre de pole.
%
%   B_r par ampere de courant magnetisant, fer infiniment permeable, rotor
%   a phi = 0, machine asynchrone 18,5 kW (48 encoches / 44 barres).
%   Trois solutions superposees :
%     - elements finis sur la geometrie reelle des encoches (reference) ;
%     - operateur condense couple aux cavites d'encoche, (n_T, n_O) = (33, 16) ;
%     - operateur condense avec la condition Phi_O = 0, meme pavage.
%   Les deux operateurs partagent la carte de la couronne et ne different
%   que par la fermeture du bloc d'ouverture.
%
%   Donnees : field_waveforms.mat, exporte de field_waveforms.npz
%   (archive, outputs/python).  Quatre vecteurs de 16384 points :
%     th        angle mecanique sur le tour complet (rad)
%     Br_fem    elements finis, encoches reelles
%     Br_cav    operateur + cavites
%     Br_op     operateur, Phi_O = 0
%   Aucune valeur n'est recalculee ici.
%
%   Le centre de pole est localise sur l'enveloppe de |B_r| puis arrondi au
%   pas d'encoche stator (7,5 deg) ; le signe est pris positif au centre.
%   Les bandes grises marquent les trois ouvertures stator du champ de vue,
%   de largeur b_s0 = 2,000 mm vue du rayon mi-entrefer R_m = 81,8887 mm.
%
%   Usage :  FIG4_MIDGAP_WAVEFORM
%
%   Auteurs : I. Laouar, A. Boukadoum, N. Mezhoud.

%% ---------------------------------------------------------------- options
USE_LATEX = true;      % false : interprete TeX de MATLAB, compatible Octave

%  Le fichier de donnees est cherche a cote du script, quel que soit le
%  dossier courant de MATLAB.
DATAFILE = fullfile(fileparts(mfilename('fullpath')),'field_waveforms.mat');
if exist(DATAFILE,'file') ~= 2
    error('FIG4:donnees', ['field_waveforms.mat est introuvable dans %s\n' ...
           'Placez-le a cote de FIG4_MIDGAP_WAVEFORM.m.'], ...
          fileparts(mfilename('fullpath')));
end

TAU_S = 7.5;           % pas d'encoche stator (deg mecaniques), 360/48
BS0   = 2.0;           % ouverture d'encoche stator (mm)
RM    = 81.8887;       % rayon mi-entrefer (mm)
NWIN  = 1.55;          % demi-fenetre tracee, en pas d'encoche

%% ---------------------------------------------------------------- donnees
S = load(DATAFILE);
th = S.th(:);  fe = S.Br_fem(:);  cv = S.Br_cav(:);  op = S.Br_op(:);
deg = th*180/pi;

%% ------------------------------------------- centrage sur un centre de pole
%  L'enveloppe de |B_r|, lissee sur 200 echantillons, est maximale au centre
%  du pole ; on arrondit ensuite au pas d'encoche le plus proche.
sm = conv(abs(fe),ones(200,1)/200,'same');
[~,imax] = max(sm);
c = TAU_S*round(deg(imax)/TAU_S);
[~,ic] = min(abs(deg-c));
sgn = sign(fe(ic));                       % signe positif au centre de pole

x = deg - c;                              % angle depuis le centre de pole
x(x >  180) = x(x >  180) - 360;
x(x < -180) = x(x < -180) + 360;
[x,ord] = sort(x);
fe = sgn*fe(ord);  cv = sgn*cv(ord);  op = sgn*op(ord);

sel = x > -NWIN*TAU_S & x < NWIN*TAU_S;
x = x(sel);  fe = fe(sel);  cv = cv(sel);  op = op(sel);

%% ---------------------------------------------------------------- couleurs
cFE  = [ 27  79 114]/255;      % bleu  : elements finis
cCAV = [ 30 132  73]/255;      % vert  : operateur + cavites
cNEU = [192  57  43]/255;      % rouge : operateur, Phi_O = 0
cBAND = [234 237 237]/255;     % gris  : bande d'ouverture stator

if USE_LATEX
    itp = 'latex';
    sxl = ['mechanical angle from the pole centre (degrees), ' ...
           'rotor at $\varphi = 0$'];
    syl = '$B_r$ at mid-gap (T per A of $I_m$)';
    sl  = {'FE, real slot geometry', ...
           'operator + slot cavities (33, 16)', ...
           'operator, $\Phi_O = 0$ (33, 16)'};
else
    itp = 'tex';
    sxl = 'mechanical angle from the pole centre (degrees), rotor at \varphi = 0';
    syl = 'B_r at mid-gap (T per A of I_m)';
    sl  = {'FE, real slot geometry', ...
           'operator + slot cavities (33, 16)', ...
           'operator, \Phi_O = 0 (33, 16)'};
end

%% ---------------------------------------------------------------- figure
f = figure('Units','inches','Position',[1 1 6.4 2.8], ...
           'Color','w','PaperPositionMode','auto');
ax = axes('Parent',f,'Box','on','TickDir','in','LineWidth',0.6, ...
          'FontName','Helvetica','FontSize',8.5,'NextPlot','add');

ymax = 1.15*max(fe);
ymin = min(0,1.1*min(fe));
ylim(ax,[ymin ymax]); xlim(ax,NWIN*TAU_S*[-1 1]);

% bandes des ouvertures stator, tracees en premier plan de fond
w = (BS0/RM)*180/pi;                       % largeur angulaire de l'ouverture
for k = [-1 0 1]
    xc = TAU_S/2 + TAU_S*k;
    fill(ax,xc + w/2*[-1 1 1 -1],[ymin ymin ymax ymax],cBAND, ...
         'EdgeColor','none','HandleVisibility','off');
end
set(ax,'Layer','top');

h(1) = plot(ax,x,fe,'-' ,'Color',cFE ,'LineWidth',1.4);
h(2) = plot(ax,x,cv,'--','Color',cCAV,'LineWidth',1.0);
h(3) = plot(ax,x,op,':' ,'Color',cNEU,'LineWidth',1.0);

text(ax,TAU_S/2,0.96*ymax,'stator opening','HorizontalAlignment','center', ...
     'VerticalAlignment','top','FontSize',6.5,'Color',[127 140 141]/255);

xlabel(ax,sxl,'Interpreter',itp,'FontSize',8.5);
ylabel(ax,syl,'Interpreter',itp,'FontSize',8.5);

lg = legend(ax,h,sl,'Interpreter',itp,'Location','southoutside','FontSize',7);
set(lg,'Box','off');
try, set(lg,'NumColumns',3); end                  %#ok<TRYNC> R2017b et au-dela

%% ---------------------------------------------------------------- export
if exist('exportgraphics','file')                % MATLAB R2020a et au-dela
    exportgraphics(f,'fig4_midgap_waveform_linear.png','Resolution',400);
    exportgraphics(f,'fig4_midgap_waveform_linear.pdf','ContentType','vector');
else
    print(f,'-dpng','-r400','fig4_midgap_waveform_linear.png');
    print(f,'-dpdf','-painters','fig4_midgap_waveform_linear.pdf');
end
end
