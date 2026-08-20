%% test_slotleak  -  Validation du banc FEM de fuite d'encoche (mec.fem_slot_leakage)
%
%  Sur une encoche RECTANGULAIRE (largeur b, hauteur h, courant uniforme),
%  la permeance de fuite du corps vaut EXACTEMENT lambda = h/(3b).
%  Le banc FEM (encoche seule, parois Neumann = fer infiniment permeable,
%  A=0 en haut) doit la retrouver a ~0,1 % pres.
%
%  L'isthme (hauteur h0, sans courant) contribue EXACTEMENT h0/b (champ
%  uniforme au-dessus de la barre) : on le soustrait analytiquement pour
%  isoler le corps. Test sur 4 rapports d'aspect h/b.
%
%  RAPPEL DU PIEGE (documente dans mec.fem_slot_leakage) : un banc avec
%  A=0 sur TOUT le contour est FAUX (le flux boucle a l'interieur et
%  traverse l'encoche 2x) -> lambda surestime de +9 a +60 %.

clear; clc;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== Validation banc FEM encoche : lambda vs h/(3b) ===\n\n');

b  = 5e-3;                       % largeur d'encoche [m]
h0 = 0.5*b;                      % isthme (sans courant), soustrait h0/b exact
fprintf('%8s | %10s %10s %8s | %7s\n','h/b','lam_FEM','h/(3b)','err %','nelem');

maxerr=0;
for ratio=[1 2 4 8]
    h=ratio*b;
    % machine factice : encoche rotor rectangulaire br1=br2=br0=b
    Mf=struct('br0',b,'hr0',h0,'br1',b,'br2',b,'hr1',h);
    Gf=struct('taur',3*b);
    id=mec.fem_slot_leakage(Mf,Gf,struct('side','rotor','Hmax',b/20));
    lam_body = id.lambda - h0/b;           % corps seul (isthme exact retire)
    lam_ref  = h/(3*b);
    err=100*(lam_body-lam_ref)/lam_ref; maxerr=max(maxerr,abs(err));
    fprintf('%8.1f | %10.4f %10.4f %+7.2f | %7d\n',ratio,lam_body,lam_ref,err,id.nelem);
end

fprintf('\nErreur max = %.2f %%  -> %s\n',maxerr,verdict(maxerr<0.5));

function s=verdict(ok), if ok, s='OK (banc valide)'; else, s='ECHEC'; end, end
