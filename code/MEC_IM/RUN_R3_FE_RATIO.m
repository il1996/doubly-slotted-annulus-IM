%% RUN_R3_FE_RATIO - verification EF du rapport lisse/encoche
%
%  BLOC R3 de SPEC_CLAUDE_CODE_v4. Objection centrale de l'Article II :
%  "le modele qui produit la correction de +5,2 % se trompe de 21,6 % sur la
%  grandeur meme que cette correction est censee corriger".
%
%  ETAT : BLOQUE PAR UNE DEPENDANCE EXTERNE.
%  Le bloc exige DEUX resolutions magnetostatiques ANSYS -- geometrie reelle
%  et geometrie lisse equivalente -- sous excitation strictement identique.
%  Le projet IM_18kW_690V.aedt ne contient QU'UN seul design (Setup1) et
%  AUCUNE geometrie lisse. Les quatre occurrences de "smooth" y sont
%  incidentes. Ce script ne peut donc pas produire le rapport EF.
%
%  CE QU'IL PRODUIT MALGRE TOUT :
%   1. la GARDE analytique, calculee : mu0*(2*pi*R*L)/g, valeur que la
%      geometrie lisse DOIT redonner a quelques dixiemes de pour cent ;
%   2. la valeur MEC a confronter, relue de R2 ;
%   3. le traitement des deux .tab, actif des qu'ils existent.
%
%  Aucun chiffre n'est invente : ce qui manque est declare manquant.
clear; clc;
diary('R3_fe_ratio_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G;
mu0=4*pi*1e-7;

fprintf('=== R3 : verification EF du rapport lisse/encoche ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    machine        : MAS 48/44, 18,5 kW\n');
fprintf('    R_s (alesage)  : %.8f m\n',G.Rs);
fprintf('    R_r (rotor)    : %.8f m\n',G.Rr);
fprintf('    entrefer g     : %.8f m = %.4f mm\n',G.Rs-G.Rr,(G.Rs-G.Rr)*1e3);
fprintf('    R mi-entrefer  : %.8f m\n',0.5*(G.Rs+G.Rr));
fprintf('    longueur L     : %.8f m\n',M.L);
fprintf('    projet EF      : C:\\Users\\hp\\Desktop\\ANSYS-\\moteur_18_5\\IM_18kW_690V.aedt\n');

%% ---- 1. LA GARDE ANALYTIQUE ------------------------------------------
%  v4 §R3 : "la geometrie lisse doit redonner, a quelques dixiemes de pour
%  cent, la permeance analytique mu0*(2*pi*R*L)/g. Si elle ne le fait pas,
%  l'extraction est fautive et le rapport ne veut rien dire."
g=G.Rs-G.Rr; Rm=0.5*(G.Rs+G.Rr);
Lam_smooth=mu0*(2*pi*Rm*M.L)/g;
fprintf('\n  ---- 1. GARDE ANALYTIQUE (a verifier sur le run lisse) ----\n');
fprintf('    Lambda_lisse = mu0*(2*pi*R*L)/g\n');
fprintf('                 = %.10e H (unite de permeance)\n',Lam_smooth);
fprintf('                 avec R = %.8f m (mi-entrefer)\n',Rm);
fprintf('    Tolerance exigee : quelques dixiemes de pour cent.\n');
%  Variante au rayon d'alesage, pour lever toute ambiguite de convention
fprintf('    Variante au rayon d''ALESAGE R_s : %.10e H\n', ...
    mu0*(2*pi*G.Rs*M.L)/g);
fprintf('    (l''ecart entre les deux conventions vaut %.4f %% :\n', ...
    100*(G.Rs-Rm)/Rm);
fprintf('     il est INFERIEUR a la tolerance, donc sans consequence.)\n');

%% ---- 2. LA VALEUR MEC A CONFRONTER -----------------------------------
fprintf('\n  ---- 2. VALEUR MEC A CONFRONTER ----\n');
if isfile('R2_kc_tiling.mat')
    S=load('R2_kc_tiling.mat');
    ip=find(S.R(:,1)==S.nTp & S.R(:,2)==S.nOp,1);
    fprintf('    source : R2_kc_tiling.mat (base chapeau, N_h = %d)\n',S.Nh);
    fprintf('    k_C au pavage de production (%d, %d) : %.6f\n', ...
        S.nTp,S.nOp,S.R(ip,6));
    fprintf('    INTERVALLE sur les neuf pavages    : %.6f a %.6f\n', ...
        min(S.R(:,6)),max(S.R(:,6)));
    fprintf('    dispersion                          : %.4f %%\n', ...
        100*(max(S.R(:,6))-min(S.R(:,6)))/mean(S.R(:,6)));
    fprintf('    Carter classique                    : %.6f\n',ctx.AGcarter.kC);
    fprintf(['    NB : R2 etablit que ce rapport n''est PAS convergé en\n' ...
             '    pavage. La confrontation EF doit porter sur l''INTERVALLE,\n' ...
             '    non sur la seule valeur de production.\n']);
else
    fprintf('    R2_kc_tiling.mat ABSENT -- executer R2 d''abord.\n');
end

%% ---- 3. TRAITEMENT DES DEUX .tab, si presents ------------------------
fprintf('\n  ---- 3. RAPPORT EF ----\n');
f_slot=fullfile('R3_fe','slotted_flux.tab');
f_smth=fullfile('R3_fe','smooth_flux.tab');
if isfile(f_slot) && isfile(f_smth)
    rd=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    a=rd(f_slot); b=rd(f_smth);
    Ps=mean(a(:,2)); Pm=mean(b(:,2));
    fprintf('    flux encoche : %.10e\n',Ps);
    fprintf('    flux lisse   : %.10e\n',Pm);
    fprintf('    RAPPORT EF   : %.6f\n',Pm/Ps);
    fprintf('    garde : lisse / analytique = %.6f  (doit valoir 1 a qq 0,1 %%)\n', ...
        Pm/Lam_smooth);
else
    fprintf('    LES DEUX .tab SONT ABSENTS -- le rapport EF n''est PAS produit.\n');
    fprintf('    Attendus : %s\n            et %s\n',f_slot,f_smth);
end

%% ---- 4. SPECIFICATION DES DEUX RUNS ---------------------------------
fprintf(['\n  ---- 4. SPECIFICATION DES DEUX RESOLUTIONS A LANCER ----\n' ...
  '    A tenir IDENTIQUE entre les deux runs :\n' ...
  '      - entrefer mecanique g = %.4f mm, rayons R_s et R_r inchanges ;\n' ...
  '      - excitation : MEME bobine, MEME FMM, MEME courant ;\n' ...
  '      - permeabilite : lineaire OU saturee, mais DECLAREE et IDENTIQUE ;\n' ...
  '      - maillage de reference, critere de convergence, vrillage ;\n' ...
  '      - meme projet IM_18kW_690V.aedt (ne PAS en creer un nouveau :\n' ...
  '        la comparabilite avec le reste de la validation en depend).\n' ...
  '    SEULE difference : encoches supprimees DES DEUX COTES dans le run\n' ...
  '    lisse, alesage et surface rotorique rendus cylindriques.\n' ...
  '    A exporter : flux totalise par la bobine, en .tab, sous\n' ...
  '      MEC_IM/R3_fe/slotted_flux.tab  et  MEC_IM/R3_fe/smooth_flux.tab\n' ...
  '    Ce script formera alors le rapport et verifiera la garde.\n'],g*1e3);

fprintf('\n=== R3 : BLOQUE, en attente des deux resolutions ANSYS ===\n');
diary off;
