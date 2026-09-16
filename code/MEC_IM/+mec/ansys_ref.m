function R = ansys_ref()
%ANSYS_REF  Référence EF (ANSYS Maxwell 2D) du moteur 18,5 kW.
%
%   R = mec.ansys_ref() renvoie les caractéristiques couple-glissement et
%   courant-glissement calculées par éléments finis (ANSYS Maxwell 2D,
%   48/44, 4 pôles, 690 V, 50 Hz), servant de référentiel de validation
%   chiffrée du modèle MEC.
%
%   R.s, R.T   : glissement et couple moyen [N.m]
%   R.s2, R.I  : glissement et courant de phase efficace [A]
%
%   ---------------------------------------------------------------------
%   VRILLAGE — CORRIGÉ le 6 août 2026. L'en-tête annonçait « skew 7,5°/
%   5 tranches ». Le projet en déclare UNE : IM_18kW_690V.aedt porte
%   UseSkewModel=true, SkewAngle='7.5deg' et NumberOfSlices='1', quatre
%   occurrences, les quatre setups. Une seule tranche ne produit aucune
%   rotation relative, donc aucun moyennage axial : la solution de champ
%   est celle d'une SECTION DROITE, NON VRILLÉE. Voir B5_skew_out.txt.
%   ---------------------------------------------------------------------
%   ANCRES — TROIS SONT FAUSSES, MESURÉES le 6 août (A8_verif_ancres_out.txt).
%   Elles sont conservées pour ne pas casser les scripts de diagnostic qui
%   les lisent, mais AUCUN CHIFFRE DES DEUX ARTICLES NE DOIT EN SORTIR.
%
%     I0  = 8,32 A   *** FAUX ***  mesuré 8,4986 A sur
%                    transitoire\a vide\Winding Plot 4.tab, 2001 points,
%                    vitesse établie 1499,942 tr/min. Écart −2,1 %.
%                    B1 emploie 8,49 : correct.
%     Rs  = 0,101 ohm *** FAUX *** B1 et le modèle donnent 0,4450 et
%                    0,4302 ohm. Facteur 4,4. À 19,73 A, l'ancre
%                    donnerait 118 W de pertes Joule statoriques contre
%                    520 W : invraisemblable.
%     Bg1 = 1,01 T   *** FAUX ***  B1 emploie 0,942 T. Écart +7,2 %.
%                    Les deux ne peuvent pas être la même grandeur.
%
%     Xm  = 46 ohm, Rfe = 1740 ohm : employées par B1, non re-mesurables
%                    depuis les .tab disponibles. Statut NON VÉRIFIÉ.
%     T_break = 324,95 N·m à s = 0,105 : NON VÉRIFIABLE comme maximum.
%                    Le maximum brut de Torque Plot 2.tab vaut 336,55 N·m
%                    à s = 0,175 — mais dans la bande 0,135–0,21, dont le
%                    bruit rms atteint 13,2 % (étendue 247,6 à 336,6 N·m).
%                    Le 324,95 appartient, lui, à la bande 0,05–0,15 dont
%                    le bruit vaut 0,64 % : il est cohérent avec la partie
%                    propre de la courbe, mais ce n'est pas le maximum du
%                    fichier. À employer comme ORDRE DE GRANDEUR seulement.
%   ---------------------------------------------------------------------
%   NATURE DE LA TABLE. R.T et R.I sont une TRANSCRIPTION d'un balayage
%   optimetrics à 129 points (Torque Plot 2.tab). Deux bandes y sont
%   bruitées : 0,135–0,21 (13,2 % rms) et 0,90–1,00 (9,9 % rms, écart-type
%   16,2 % sur 21 points, de 56,9 à 126,8 N·m). Les 30 points transcrits
%   ne sont exploitables quantitativement que pour s <= 0,13.
%   NE PAS « corriger » T(0,15) = 297,0 : le fichier y donne 307,77, mais
%   remplacer un échantillon de bruit par un autre n'améliore rien.

% Couple-glissement (extrait bas glissement, régime moteur)
R.s = [0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 ...
       0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.10 ...
       0.105 0.11 0.115 0.12 0.15 0.20 0.30 0.50 0.70 1.00];
R.T = [35.13 68.07 98.98 128.65 156.46 181.59 204.09 224.22 242.81 257.26 ...
       271.38 283.67 293.49 301.92 308.33 313.79 317.48 320.48 323.01 324.79 ...
       324.95 323.40 320.20 315.58 297.0 305.0 231.0 153.8 117.8 98.98];

% Courant-glissement (phase, efficace)
R.s2 = R.s;
R.I  = [9.71 12.87 16.65 20.66 24.65 28.59 32.42 36.02 39.45 42.78 ...
        45.97 48.94 51.70 54.42 57.00 59.16 60.76 60.77 65.12 67.03 ...
        68.88 70.42 71.89 73.39 81.0 87.0 95.1 100.5 102.8 108.22];

% Ancres ponctuelles
%  *** FAUX, mesuré le 6 août — voir l'en-tête. Conservés pour les scripts
%  *** de diagnostic qui les lisent ; interdits aux deux articles.
R.anchor.Bg1 = 1.01;   R.anchor.Bg1_FAUX = true;   % mesuré : 0,942 T
R.anchor.I0  = 8.32;   R.anchor.I0_FAUX  = true;   % mesuré : 8,4986 A
R.anchor.Rs  = 0.101;  R.anchor.Rs_FAUX  = true;   % employé : 0,4450 ohm
%  non vérifiées, mais employées par B1
R.anchor.Xm = 46;     R.anchor.Rfe = 1740;
R.anchor.T_rated = 128.65; R.anchor.s_rated = 0.02;
R.anchor.T_break = 324.95; R.anchor.s_break = 0.105;
end
