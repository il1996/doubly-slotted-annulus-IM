%% RUN_LEAKAGE  -  Saturation des chemins de fuite : diagnostic AVANT modelisation
%
%  Hypothese a tester : la saturation des chemins de fuite (becs de dents,
%  pont d'encoche rotorique) expliquerait le residu au fort glissement.
%
%  ATTENTION AU SENS : la saturation REDUIT la permeance de fuite -> Xsigma
%  BAISSE -> couple et courant AUGMENTENT. Or le modele les SURESTIME deja.
%  Ce script tranche quantitativement avant d'ecrire le moindre modele :
%    A) on remonte Rr et Xsigma IMPLIQUES par l'EF (a partir de ses I et T) ;
%    B) on les compare au modele ;
%    C) on verifie si les becs rotoriques saturent reellement (B_bec).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== Chemins de fuite : diagnostic au fort glissement ===\n\n');

M=mec.machine_18_5kW(); ctx=mec.build_context(M);
G=ctx.G; W=ctx.W; ref=mec.ansys_ref();
mu0=4*pi*1e-7; Uph=M.Uph; Rs=ctx.Rs;

%% A/B) Rr et Xsigma impliques par l'EF vs modele
%  A fort glissement le courant magnetisant est faible : I2 ~ I1, d'ou
%     Rr/s = T*omega/(m*p*I1^2)      et    Xsigma = sqrt((Uph/I1)^2-(Rs+Rr/s)^2)
fprintf('--- A/B) Parametres impliques par l''EF vs modele ---\n');
fprintf('%6s | %8s %8s %7s | %8s %8s %7s\n',...
    's','Rr_EF','Rr_MEC','ecart%','Xsig_EF','Xsig_MEC','ecart%');
for s=[0.2 0.5 1.0]
    Ife=interp1(ref.s2,ref.I,s); Tfe=interp1(ref.s,ref.T,s);
    RrS_fe = Tfe*M.w/(M.m*M.p*Ife^2);          % = Rr/s implique par l'EF
    Rr_fe  = s*RrS_fe;
    Xs_fe  = sqrt(max((Uph/Ife)^2 - (Rs+RrS_fe)^2, 0));
    % modele
    r=mec.equivalent_circuit(ctx,s,ctx.Xm0);
    Cg=mec.cage(M,G,W,ctx.Lk,s);
    Rr_me=Cg.Rr;
    Xs_me = sqrt(max((Uph/r.I1)^2 - (Rs+Rr_me/s)^2, 0));
    fprintf('%6.2f | %8.4f %8.4f %6.1f | %8.3f %8.3f %6.1f\n',...
        s, Rr_fe, Rr_me, 100*(Rr_me-Rr_fe)/Rr_fe, ...
           Xs_fe, Xs_me, 100*(Xs_me-Xs_fe)/Xs_fe);
end

%% C) Les becs rotoriques saturent-ils vraiment ?
fprintf('\n--- C) Induction dans le bec/pont rotorique (saturation active ?) ---\n');
fprintf('%6s | %10s %10s | %8s\n','s','Ibar[A crete]','Phi_bec[mWb]','B_bec[T]');
for s=[0.2 0.5 1.0]
    r=mec.equivalent_circuit(ctx,s,ctx.Xm0);
    Ibar = r.I2*(2*M.m*W.kw1*W.Nph)/M.Nr;      % courant de barre efficace
    Ipk  = sqrt(2)*Ibar;
    % flux de fuite traversant l'ouverture d'encoche (air, br0) sur la hauteur hr0
    Phi  = mu0*Ipk*(M.hr0*M.L)/M.br0;
    Btip = Phi/(M.hr0*M.L);                    % induction dans le fer du bec
    fprintf('%6.2f | %10.0f %10.3f | %8.2f\n', s, Ipk, Phi*1e3, Btip);
end

fprintf('\n--- Conclusion ---\n');
fprintf('Si Xsig_MEC < Xsig_EF : la fuite du modele est TROP FAIBLE.\n');
fprintf('La saturation ne pourrait que la REDUIRE encore -> mauvais sens.\n');
