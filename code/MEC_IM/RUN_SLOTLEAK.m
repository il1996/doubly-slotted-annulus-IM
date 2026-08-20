%% RUN_SLOTLEAK  -  Identification FEM de la permeance de fuite d'encoche
%
%  Le diagnostic RUN_LEAKAGE montre que la reactance de fuite du modele est
%  8-13 % TROP FAIBLE au fort glissement. Le maillon faible identifie est le
%  terme empirique 0.4*(br1/br0) ajoute dans lam_r_tip (mec.leakage).
%  On le remplace ici par une IDENTIFICATION FEM (meme demarche que
%  mec.fem_airgap_ident pour sigma0 : mesurer au lieu de supposer).
%
%  Banc : encoche seule, parois = fer infiniment permeable (Neumann),
%  A=0 en haut. Permeance par l'energie : lambda = 2W/(mu0*I^2).
%  VALIDE a 0,0-0,1 % contre h/(3b) sur encoche rectangulaire (test_slotleak).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== Identification FEM de la permeance de fuite d''encoche ===\n\n');

M=mec.machine_18_5kW(); G=mec.geometry(M);

fprintf('%8s | %9s %9s %8s | %7s\n','cote','lam_FEM','lam_MEC','ecart%','nelem');
res=struct();
for side={'rotor','stator'}
    s=side{1};
    id=mec.fem_slot_leakage(M,G,struct('side',s));
    res.(s)=id;
    fprintf('%8s | %9.4f %9.4f %7.1f | %7d\n', s, id.lambda, id.lambda_ana, ...
        100*(id.lambda_ana-id.lambda)/id.lambda, id.nelem);
end

fprintf('\n--- Detail rotor ---\n');
lam_body = M.hr1/(3*0.5*(M.br1+M.br2));
lam_ist  = M.hr0/M.br0;
lam_emp  = 0.4*(M.br1/M.br0);
fprintf('analytique = corps %.3f + isthme %.3f + terme EMPIRIQUE %.3f = %.3f\n',...
    lam_body, lam_ist, lam_emp, lam_body+lam_ist+lam_emp);
fprintf('FEM        = %.3f\n', res.rotor.lambda);
fprintf('=> le terme empirique 0.4*(br1/br0) vaut en realite %.3f (FEM - corps - isthme)\n',...
    res.rotor.lambda - lam_body - lam_ist);
