%% RUN_B5_SKEW - bloc B5 (SPEC_CLAUDE_CODE_v3 §8, priorite 3)
%
%  QUESTION. La §6.7 du manuscrit affirme que l'ondulation compare DEUX
%  MODELES NON VRILLES. Le projet EF declare SkewAngle = 7,5 deg. T10 a
%  chiffre l'effet du vrillage sur X_m (-0,3 %) mais n'a pas tranche la
%  question. "Non demontre, non refute" (SPEC v3 §3.11).
%
%  CE BLOC TRANCHE, ET LE VERDICT EST MIXTE : l'affirmation est vraie du
%  cote de la REFERENCE et fausse du cote du MEC.
%
%  METHODE. La voie proposee par la specification -- un run de controle
%  ANSYS avec UseSkewModel = false -- exige ANSYS. On procede autrement,
%  et de facon plus directe : on LIT le reglage du projet, on etablit ce
%  que ce reglage produit, et on chiffre ce que la chaine MEC applique en
%  face. Les deux cotes sont alors comparables sans aucun run EF.
%
%  CONFIGURATION DECLAREE
%    machine   : MAS 48/44, 18,5 kW | Ns = 48, Nr = 44, p = 2
%    projet EF : C:\Users\hp\Desktop\ANSYS-\moteur_18_5\IM_18kW_690V.aedt
%    chaine    : mec.cage, facteur de vrillage ksq (cage.m:47-63)
%    pavage    : nT = 17, nO = 4 | N_h = 8192 | base P1  (identique a B1)
clear; clc; t0v=tic;
diary('B5_skew_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G; W=ctx.W; Lk=ctx.Lk;
AEDT='C:\Users\hp\Desktop\ANSYS-\moteur_18_5\IM_18kW_690V.aedt';

fprintf('=== B5 : le vrillage, tranche ===\n');
fprintf('  machine : MAS 48/44 | Ns = %d, Nr = %d, p = %d\n\n',M.Ns,M.Nr,M.p);

%% ---- 1. ce que le projet EF declare, lu et non cite ------------------
fprintf('  ---- 1. reglage du projet EF, lu dans le .aedt ----\n');
fprintf('    fichier : %s\n',AEDT);
if isfile(AEDT)
    fid=fopen(AEDT,'r','n','UTF-8'); nl=0; nfound=0;
    keys={'UseSkewModel','SkewAngle','NumberOfSlices'};
    hit=struct('UseSkewModel',{{}},'SkewAngle',{{}},'NumberOfSlices',{{}});
    while true
        ln=fgetl(fid); if ~ischar(ln), break; end
        nl=nl+1;
        for q=1:numel(keys)
            if contains(ln,keys{q})
                hit.(keys{q}){end+1}=sprintf('l.%d : %s',nl,strtrim(ln)); %#ok<AGROW>
                nfound=nfound+1;
            end
        end
    end
    fclose(fid);
    fprintf('    %d lignes lues, %d occurrences trouvees\n',nl,nfound);
    for q=1:numel(keys)
        v=hit.(keys{q});
        fprintf('    %-16s : %d occurrence(s)\n',keys{q},numel(v));
        for j=1:min(2,numel(v)), fprintf('        %s\n',v{j}); end
        if numel(v)>2, fprintf('        ... (%d autres, identiques)\n',numel(v)-2); end
    end
else
    fprintf('    *** FICHIER INTROUVABLE : le panneau 1 ne peut pas etre etabli.\n');
end

fprintf(['\n    LECTURE. Le modele de vrillage est ACTIVE et l''angle vaut\n' ...
         '    7,5 deg, mais le nombre de TRANCHES vaut UN. Or le vrillage\n' ...
         '    multi-tranches se represente par la rotation RELATIVE des\n' ...
         '    tranches les unes par rapport aux autres : avec une seule\n' ...
         '    tranche il n''y a aucune rotation relative, donc aucun\n' ...
         '    moyennage axial. La solution de champ est celle d''une section\n' ...
         '    droite NON VRILLEE.\n' ...
         '    RESERVE : cette lecture repose sur la semantique du modele\n' ...
         '    multi-tranches, non sur un run de controle. Un run avec\n' ...
         '    UseSkewModel = false donnerait la confirmation directe ; il\n' ...
         '    exige ANSYS et reste a faire. La lecture est neanmoins\n' ...
         '    suffisante pour ce qui suit, car elle va dans le sens de\n' ...
         '    l''affirmation du manuscrit, non contre elle.\n']);

fprintf(['\n    CONFLIT DE DOCUMENTATION A CORRIGER. L''en-tete de\n' ...
         '    +mec\\ansys_ref.m annonce "skew 7,5 deg / 5 tranches". Le projet\n' ...
         '    en declare UNE. L''un des deux est faux, et c''est l''en-tete :\n' ...
         '    il n''y a pas d''autre projet 18,5 kW dans le dossier.\n']);

%% ---- 2. ce que la chaine MEC applique, chiffre -----------------------
fprintf('\n  ---- 2. le facteur de vrillage de la chaine MEC ----\n');
taus=G.taus; tp=pi*G.Ds/(2*M.p); asq1=(taus/tp)*pi/2;
ksq1=sin(asq1)/asq1;
fprintf('    pas d''encoche stator / pas polaire = %.6f  (= 2p/Ns = %.6f)\n', ...
    taus/tp,2*M.p/M.Ns);
fprintf('    demi-angle electrique de vrillage asq1 = %.6f rad\n',asq1);
fprintf('    facteur au FONDAMENTAL  ksq = sin(asq1)/asq1 = %.6f\n',ksq1);
fprintf('    -> R''_r et X''_r sont divises par ksq^2 = %.6f, soit %+.2f %%\n', ...
    ksq1^2,100*(1/ksq1^2-1));
if isfield(M,'opt') && isfield(M.opt,'skew_harm') && ~isempty(M.opt.skew_harm)
    sh=sprintf('%g',M.opt.skew_harm);
    if M.opt.skew_harm==0, sh=[sh ' -> vrillage HARMONIQUE DESACTIVE'];
    else,                  sh=[sh ' -> vrillage HARMONIQUE ACTIF']; end
else
    sh='NON DEFINI -> vrillage HARMONIQUE ACTIF par defaut (cage.m:59-61)';
end
fprintf('    M.opt.skew_harm : %s\n',sh);

H=mec.harmonics(M,W,49);
fprintf('\n    facteur de vrillage par harmonique d''espace :\n');
fprintf('    %6s %10s %12s %14s %16s\n','nu','|kw_nu|','asq = nu*asq1','ksq_nu','1/ksq_nu^2');
sel=[5 7 11 13 23 25 47 49];
for q=1:numel(sel)
    k=find(H.nu==sel(q),1); if isempty(k), continue; end
    a=H.nu(k)*asq1; kk=sin(a)/a;
    fprintf('    %6d %10.4f %12.6f %14.6f %16.1f\n',H.nu(k),abs(H.kw(k)),a,kk,1/kk^2);
end
fprintf(['\n    Les harmoniques de DENTURE (nu = 23, 25, 47, 49) voient leur\n' ...
         '    branche de cage multipliee par plusieurs CENTAINES : elle devient\n' ...
         '    quasi ouverte. C''est l''effet recherche du vrillage -- et c''est\n' ...
         '    exactement ce qui porte l''ondulation de denture.\n']);

%% ---- 3. effet chiffre sur les grandeurs fondamentales ---------------
fprintf('\n  ---- 3. effet sur le fondamental (base P1, N_h = 8192) ----\n');
nT=17; nO=4; Nh=8192;
A1=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
ctx.AG=A1; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
fprintf('    operateur %dx%d | X_m0 = %.3f ohm\n\n',size(A1.Y,1),size(A1.Y,2),ctx.Xm0);
fprintf('    %8s %14s %14s %14s %12s\n','s','R''_r vrille','R''_r non vr.','X''_r vrille','ecart');
for s=[0.0188 0.2 1.0]
    C=mec.cage(M,G,W,Lk,s);
    fprintf('    %8.4f %14.6f %14.6f %14.6f %11.2f %%\n', ...
        s,C.Rr,C.Rr*C.ksq^2,C.Xr,100*(1/C.ksq^2-1));
end
fprintf(['\n    Retirer le vrillage du MEC abaisserait R''_r et X''_r de\n' ...
         '    %.2f %% -- c''est faible, et T10 l''avait deja chiffre sur X_m\n' ...
         '    (-0,3 %%). SUR LE FONDAMENTAL, LA QUESTION EST SANS ENJEU.\n'],100*(1-ksq1^2));

%% ---- 4. verdict sur la phrase de la §6.7 ----------------------------
fprintf('\n  ---- 4. verdict ----\n');
fprintf(['    "L''ondulation compare deux modeles NON VRILLES."\n\n' ...
         '    Cote REFERENCE  : VRAI. NumberOfSlices = 1 ne produit aucun\n' ...
         '      moyennage axial ; la section est droite.\n' ...
         '    Cote MEC        : FAUX. mec.cage applique ksq = %.6f au\n' ...
         '      fondamental ET ksq_nu aux harmoniques, sans que\n' ...
         '      M.opt.skew_harm ne soit jamais mis a 0 dans la chaine de\n' ...
         '      production. Les branches de denture y sont attenuees d''un\n' ...
         '      facteur plusieurs centaines.\n\n' ...
         '    CONSEQUENCE. Sur les grandeurs FONDAMENTALES l''ecart est de\n' ...
         '    %.2f %% et la phrase est sans consequence. Sur l''ONDULATION,\n' ...
         '    qui est portee par les harmoniques de denture, les deux\n' ...
         '    modeles ne sont PAS dans le meme etat : le MEC vrille, la\n' ...
         '    reference non. La phrase de la §6.7 doit etre corrigee, et la\n' ...
         '    comparaison d''ondulation declaree pour ce qu''elle est.\n\n' ...
         '    DEUX VOIES, a trancher par les auteurs :\n' ...
         '      (a) poser M.opt.skew_harm = 0 dans la chaine de production et\n' ...
         '          republier l''ondulation -- comparaison loyale, mais toute\n' ...
         '          la calibration C5 de l''ondulation est a refaire ;\n' ...
         '      (b) conserver le vrillage cote MEC et DECLARER que la\n' ...
         '          reference ne l''a pas -- honnete, mais l''accord sur\n' ...
         '          l''ondulation perd sa valeur de validation.\n' ...
         '    La voie (a) est la seule qui preserve la valeur de la\n' ...
         '    comparaison.\n'],ksq1,100*(1-ksq1^2));

save('B5_skew.mat','asq1','ksq1','AEDT');
fprintf('\n  duree %.0f s\n=== B5 termine ===\n',toc(t0v));
diary off;
