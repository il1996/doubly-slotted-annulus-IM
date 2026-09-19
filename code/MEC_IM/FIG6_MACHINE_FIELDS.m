function FIG6_MACHINE_FIELDS
%FIG6_MACHINE_FIELDS  Figure 6 : induction au mi-entrefer sur un quadrant,
%                     reseau contre reference elements finis.
%
%   Composantes radiale B_r et tangentielle B_t au mi-entrefer, de 0 a 90
%   degres mecaniques, machine asynchrone 18,5 kW (48 encoches / 44 barres).
%   Colonne (a) a vide, colonne (b) au point nominal.
%
%   Trois traces par panneau :
%     - reference elements finis, transitoire ANSYS 2 D, une tranche, au
%       dernier instant du calcul ;
%     - reseau avec l'operateur condense, condition Phi_O = 0, (17, 4) ;
%     - reseau avec l'operateur couple aux cavites d'encoche, (33, 16).
%   Les deux champs de reseau sont recales sur la seule phase du
%   fondamental, sans aucun autre decalage ajuste.
%
%   Donnees : fig6_data.mat, exporte des transcripts
%     Z2_fields_avide.txt / Z2_fields_charge.txt      (Phi_O = 0)
%     Z4_fields_cav_avide.txt / Z4_fields_cav_charge.txt  (cavites)
%   Chaque cas porte 1000 points :
%     th_<cas>              angle mecanique (rad)
%     Bre_<cas>, Bte_<cas>  reference EF, radial et tangentiel (T)
%     Brm_<cas>, Btm_<cas>  reseau, Phi_O = 0
%     Brc_<cas>, Btc_<cas>  reseau + cavites
%   La reference EF est la meme dans les deux transcripts, verifie a
%   l'export.  Aucune valeur n'est recalculee ici.
%
%   Le pas de temps de la reference, 1 ms, ne resout pas les harmoniques
%   d'encoche rotorique : la comparaison des rangs eleves porte cette
%   reserve (voir le texte).
%
%   Usage :  FIG6_MACHINE_FIELDS
%
%   Auteurs : I. Laouar, A. Boukadoum, N. Mezhoud.

%% ---------------------------------------------------------------- options
USE_LATEX = true;      % false : interprete TeX de MATLAB, compatible Octave

%  Le fichier de donnees est cherche a cote du script, quel que soit le
%  dossier courant de MATLAB.
DATAFILE = fullfile(fileparts(mfilename('fullpath')),'fig6_data.mat');
if exist(DATAFILE,'file') ~= 2
    error('FIG6:donnees', ['fig6_data.mat est introuvable dans %s\n' ...
           'Placez-le a cote de FIG6_MACHINE_FIELDS.m.'], ...
          fileparts(mfilename('fullpath')));
end

DEG_MAX = 90;          % quadrant trace

%% ---------------------------------------------------------------- donnees
D = load(DATAFILE);
cas   = {'avide','charge'};
titre = {'(a) no load','(b) rated load'};

%% ---------------------------------------------------------------- couleurs
cFE  = [ 27  79 114]/255;      % bleu  : elements finis
cNEU = [192  57  43]/255;      % rouge : reseau, Phi_O = 0
cCAV = [ 30 132  73]/255;      % vert  : reseau + cavites

if USE_LATEX
    itp = 'latex';
    syl = {'$B_r$ at mid-gap (T)','$B_t$ at mid-gap (T)'};
    sl  = {'FE (transient, single slice)', ...
           'network, operator $\Phi_O=0$ (17, 4)', ...
           'network, operator + cavities (33, 16)'};
else
    itp = 'tex';
    syl = {'B_r at mid-gap (T)','B_t at mid-gap (T)'};
    sl  = {'FE (transient, single slice)', ...
           'network, operator \Phi_O = 0 (17, 4)', ...
           'network, operator + cavities (33, 16)'};
end

%% ---------------------------------------------------------------- figure
f = figure('Units','inches','Position',[1 1 6.8 4.4], ...
           'Color','w','PaperPositionMode','auto');

for i = 1:2                                     % colonne : a vide / en charge
    c   = cas{i};
    deg = D.(['th_' c])*180/pi;
    sel = deg <= DEG_MAX;

    %  ligne 1 : composante radiale ; ligne 2 : composante tangentielle
    ye = {D.(['Bre_' c]), D.(['Bte_' c])};      % reference EF
    ym = {D.(['Brm_' c]), D.(['Btm_' c])};      % reseau, Phi_O = 0
    yc = {D.(['Brc_' c]), D.(['Btc_' c])};      % reseau + cavites

    for j = 1:2
        ax = subplot(2,2,(j-1)*2+i); hold(ax,'on');
        set(ax,'Box','on','TickDir','in','LineWidth',0.6, ...
               'FontName','Helvetica','FontSize',8.5);

        h(1) = plot(ax,deg(sel),ye{j}(sel),'-' ,'Color',cFE ,'LineWidth',0.9);
        h(2) = plot(ax,deg(sel),ym{j}(sel),':' ,'Color',cNEU,'LineWidth',0.7);
        h(3) = plot(ax,deg(sel),yc{j}(sel),'--','Color',cCAV,'LineWidth',0.8);

        xlim(ax,[0 DEG_MAX]);
        ylabel(ax,syl{j},'Interpreter',itp,'FontSize',8.5);
        if j == 2
            xlabel(ax,'mechanical angle (degrees)','FontSize',8.5);
        else
            title(ax,titre{i},'FontSize',8.5);
        end
        if i == 1 && j == 1
            lg = legend(ax,h,sl,'Interpreter',itp, ...
                        'Location','southwest','FontSize',6.2);
            set(lg,'Box','off');
        end
    end
end

%% ---------------------------------------------------------------- export
if exist('exportgraphics','file')                % MATLAB R2020a et au-dela
    exportgraphics(f,'fig6_machine_fields.png','Resolution',400);
    exportgraphics(f,'fig6_machine_fields.pdf','ContentType','vector');
else
    print(f,'-dpng','-r400','fig6_machine_fields.png');
    print(f,'-dpdf','-painters','fig6_machine_fields.pdf');
end
end
