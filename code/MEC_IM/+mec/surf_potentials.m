function [Us, Ur] = surf_potentials(me, Se)
%SURF_POTENTIALS  Potentiels des noeuds de surface stator et rotor.
%
%   [Us,Ur] = mec.surf_potentials(me,Se) extrait du vecteur solution Se.U
%   les potentiels des Ms noeuds de surface STATOR et des Mr noeuds de
%   surface ROTOR, dans l'ordre ou mec.mesh_refined les a crees
%   (mesh_refined.m:273, champ me.gapF.ids).
%
%   BLOC C1. Etait la fonction LOCALE surfU de RUN_ARTICLE.m (l. 514-516),
%   donc non appelable ailleurs. Promue ici. Renommee surf_potentials :
%   surfU etait ambigu (potentiel ou surface ?).
%
%   Voir aussi : mec.mesh_refined, mec.solve_mesh, mec.inst_currents.

Us = Se.U(me.gapF.ids(1:me.Ms));
Ur = Se.U(me.gapF.ids(me.Ms+1:end));
end
