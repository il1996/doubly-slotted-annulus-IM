function [i3, ibar] = inst_currents(M, W, r)
%INST_CURRENTS  Courants instantanes triphases et de barre au point r.
%
%   [i3,ibar] = mec.inst_currents(M,W,r) construit, a partir du resultat
%   de mec.equivalent_circuit :
%     i3   : 3x1, courants de phase instantanes a l'instant de reference
%            (phase du fondamental complexe r.I1c) ;
%     ibar : Nr x 1, courants de barre instantanes, repartis sur le pas
%            rotorique et cales sur la phase de r.I2c.
%
%   BLOC C1. Cette fonction etait une fonction LOCALE de RUN_ARTICLE.m
%   (l. 508-513), donc non appelable depuis un autre script : le bloc B7 a
%   du la dupliquer pour produire les sondes de la Table 17. Elle est
%   promue ici pour qu'une seule definition existe.
%
%   Voir aussi : mec.equivalent_circuit, mec.mesh_refined, mec.regional_max.

psi1 = angle(r.I1c);
psi2 = angle(r.I2c);
i3   = sqrt(2)*r.I1*[cos(psi1); cos(psi1-2*pi/3); cos(psi1+2*pi/3)];
Ibar = r.I2*(2*M.m*W.kw1*W.Nph)/M.Nr;
ibar = sqrt(2)*Ibar*cos(M.p*(2*pi*(0:M.Nr-1)'/M.Nr) - psi2);
end
