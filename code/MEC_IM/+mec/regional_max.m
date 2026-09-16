function [Bts, Bys, Btr, Byr] = regional_max(me, Se)
%REGIONAL_MAX  Maxima d'induction par region du maillage raffine.
%
%   [Bts,Bys,Btr,Byr] = mec.regional_max(me,Se) renvoie, en tesla, le
%   maximum de |B| sur les branches de FER de chacune des quatre regions :
%     Bts : dents stator      Bys : culasse stator
%     Btr : dents rotor       Byr : culasse rotor
%   Ce sont les quatre sondes locales de la Table 17.
%
%   BLOC C1. Etait la fonction LOCALE regional_max de RUN_ARTICLE.m
%   (l. 531-539), donc non appelable ailleurs : le bloc B7 a du la
%   dupliquer. Promue ici.
%
%   RESERVE, ETABLIE PAR B7. Les sondes ainsi obtenues ne sortent PAS de
%   la meme chaine d'entrefer que le reseau de performance :
%   mec.mesh_refined ferme la couronne par mec.airgap_fourier en base P0
%   IMPOSEE (airgap_fourier.m:63, appele sans argument par
%   mesh_refined.m:271) a N_h = 100 par defaut (mesh_refined.m:258), la ou
%   le reseau de performance emploie mec.airgap_dtn_tooth en base P1 a
%   N_h = 8192. Les sondes derivent de jusqu'a 5,55 % sur un facteur 8 de
%   troncature, le maximum sur la dent rotor. Toute publication de ces
%   valeurs doit declarer chaine, base, troncature et derive.
%
%   Voir aussi : mec.mesh_refined, mec.solve_mesh.

nb   = me.gapfirst-1;
isFe = logical(me.iron(1:nb));
a    = me.a(1:nb);
B    = abs(Se.B(1:nb));
Ms   = me.Ms; Ls = me.Ls; Mr = me.Mr; nst = me.nr; nrt = me.nr;

laS = ceil(a/Ms);      laS(a > Ms*Ls) = 0;
aR  = a - Ms*Ls;
laR = ceil(aR/Mr);     laR(a <= Ms*Ls) = 0;

Bts = max(B(isFe & laS>=1 & laS<=nst));
Bys = max(B(isFe & laS>nst));
Btr = max(B(isFe & laR>=1 & laR<=nrt));
Byr = max(B(isFe & laR>nrt));
end
