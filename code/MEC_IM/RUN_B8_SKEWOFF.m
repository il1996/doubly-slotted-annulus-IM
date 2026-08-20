%% RUN_B8_SKEWOFF - exposition de l'Article II au facteur de vrillage
%
%  DECISION 2 (note du 6 aout 2026). La reference EF est une tranche
%  DROITE : IM_18kW_690V.aedt declare UseSkewModel=true, SkewAngle=7,5 deg,
%  mais NumberOfSlices='1' (quatre occurrences, les quatre setups). Une
%  seule tranche ne produit aucune rotation relative, donc aucun moyennage
%  axial : la solution de champ est celle d'une section non vrillee.
%
%  Or mec.cage applique ksq_nu = sin(nu*asq1)/(nu*asq1) aux branches
%  HARMONIQUES, et M.opt.skew_harm n'est jamais defini dans la chaine de
%  production : le vrillage harmonique est donc ACTIF par defaut
%  (cage.m:59-61). A nu = 23 : 1/ksq^2 = 532.
%
%  CE QUI EST EXPOSE. Contrairement a ce qu'annoncait le bloc B5, la carte
%  d'ondulation ne l'est PAS : RUN_TRANSIENT_1SLICE impose des courants de
%  barre purement fondamentaux (l. 50-58) et n'appelle aucune branche
%  harmonique. C'est le SCHEMA EQUIVALENT qui l'est : equivalent_circuit
%  met les branches harmoniques EN SERIE dans le circuit statorique
%  (l. 40-51, 63-71) et elles alimentent Tem (l. 114), Pcu_r (l. 116),
%  res.Ibar (l. 155) et, par Ztot, I1, E1 et cos phi.
%
%  Donc : tab:im_ec, tab:im_tests et tab:price.
%
%  RESERVE A DECLARER. Au fondamental l'interrupteur ne fait rien --
%  cage.m:61 donne asq = asq1 dans les deux cas a nu = 1 -- et ksq =
%  0,997147 reste applique a R'_r et X'_r, soit +0,57 % que la reference
%  n'a pas. Ce residu est irreductible sans toucher au code.
%
%  SENS ATTENDU : defavorable. Avec vrillage les branches harmoniques sont
%  quasi ouvertes et Tem_h ~ 0 ; sans lui elles deviennent actives et
%  Tem_h < 0 (equivalent_circuit.m:110). Le MEC est deja bas de 5,7 % sur
%  le couple : la comparaison loyale va creuser l'ecart. C'est le resultat
%  a etablir, pas a eviter.
%
%  CONFIGURATION : identique a B1 -- pavage nT=17 nO=4, N_h=8192, base P1.
%  REFERENCES EF : famille TRANSITOIRE (decision 1) --
%    en charge 121,63 N.m / 19,73 A ; a vide 8,49 A ; calage 104,31 N.m.
%
%  Sortie : B8_skewoff_out.txt
clear; clc; t0=tic;
diary('B8_skewoff_out.txt'); diary on;

nT=17; nO=4; Nh=8192; s_ch=0.0188;

fprintf('=== B8 : exposition au facteur de vrillage ===\n');
fprintf('  pavage nT=%d nO=%d | N_h = %d | base P1\n',nT,nO,Nh);
fprintf('  references EF : famille TRANSITOIRE (decision 1 du 6 aout)\n');
fprintf('    en charge 121.63 N.m / 19.73 A | a vide 8.49 A | calage 104.31 N.m\n\n');

%% ---- les deux configurations ----------------------------------------
cfg = {struct('name','vrillage ACTIF  (production)','skew',[]), ...
       struct('name','vrillage HARM. NEUTRALISE  ','skew',0)};
R = cell(2,1);

for c = 1:2
    M = mec.machine_18_5kW();
    if ~isempty(cfg{c}.skew), M.opt.skew_harm = cfg{c}.skew; end
    ctx = mec.build_context(M); G = ctx.G;
    % l'operateur ne depend pas de la cage : meme AG dans les deux cas
    A1 = mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
    ctx.AG = A1; ctx.Xm0 = mec.magnetizing(ctx,0.2).Xm;
    % securite : la valeur vue par mec.cage vient de ctx.M
    if ~isempty(cfg{c}.skew), ctx.M.opt.skew_harm = cfg{c}.skew; end

    r  = mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);   % en charge
    r0 = mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);   % a vide
    rb = mec.equivalent_circuit(ctx,1.0 ,ctx.Xm0);   % rotor bloque

    kbar = (2*M.m*ctx.W.kw1*ctx.W.Nph)/M.Nr;
    S = struct();
    S.Xm0    = ctx.Xm0;      S.Xm     = r.Xm;
    S.Rr     = r.Rr;         S.Xr     = r.Xr;
    S.Tem    = r.Tem;        S.Tem1   = r.Tem1;   S.Tem_h = r.Tem_h;
    S.I1     = r.I1;         S.Ibar   = r.Ibar;
    S.cosphi = r.cosphi;     S.eta    = r.eta;    S.Pfe   = r.Pfe;
    S.Pcu_r  = r.Pcu_r;      S.Pcu_r_h= r.Pcu_r_h;
    S.I0     = r0.I1;        S.E1     = r0.E1;
    S.Tb     = rb.Tem;       S.I1b    = rb.I1;
    S.Ibarb  = rb.I2*kbar;   S.Rrb    = rb.Rr;
    S.ksq    = r.Xm*0;                          % place tenue, rempli ci-dessous
    Cg = mec.cage(M,G,ctx.W,ctx.Lk,s_ch);  S.ksq = Cg.ksq;
    S.Th     = r.T_h(:);     S.nu     = r.nu(:);
    R{c} = S;
    fprintf('  %s : fait\n',cfg{c}.name);
end
fprintf('\n');

A = R{1}; B = R{2};
d = @(a,b) 100*(b-a)/max(abs(a),eps);

%% ---- 1. ce que le vrillage harmonique deplace ------------------------
fprintf('  ---- 1. EFFET DE LA NEUTRALISATION, GRANDEUR PAR GRANDEUR ----\n');
fprintf('  %-26s %14s %14s %10s\n','grandeur','vrillage ON','vrillage OFF','variation');
p2=@(l,a,b)fprintf('  %-26s %14.4f %14.4f %9.2f %%\n',l,a,b,d(a,b));
p2('couple net [N.m]',        A.Tem,  B.Tem);
p2('  dont fondamental',      A.Tem1, B.Tem1);
p2('  dont parasite',         A.Tem_h,B.Tem_h);
p2('courant I1 [A]',          A.I1,   B.I1);
p2('courant de barre [A]',    A.Ibar, B.Ibar);
p2('cos(phi)',                A.cosphi,B.cosphi);
p2('rendement',               A.eta,  B.eta);
p2('pertes Joule rotor [W]',  A.Pcu_r,B.Pcu_r);
p2('  dont harmoniques',      A.Pcu_r_h,B.Pcu_r_h);
p2('pertes fer [W]',          A.Pfe,  B.Pfe);
p2('courant a vide [A]',      A.I0,   B.I0);
p2('f.e.m. entrefer [V]',     A.E1,   B.E1);
p2('couple au calage [N.m]',  A.Tb,   B.Tb);
p2('I1 au calage [A]',        A.I1b,  B.I1b);
p2('I barre au calage [A]',   A.Ibarb,B.Ibarb);
fprintf('  (X_m0 = %.3f ohm dans les deux cas : la cage n''entre pas\n',A.Xm0);
fprintf('   dans la branche magnetisante -- controle)\n');

%% ---- 2. couples parasites, harmonique par harmonique -----------------
fprintf('\n  ---- 2. COUPLES PARASITES PAR HARMONIQUE [N.m] ----\n');
fprintf('  %6s %16s %16s\n','nu','vrillage ON','vrillage OFF');
for k=1:numel(A.nu)
    fprintf('  %6d %16.5f %16.5f\n',A.nu(k),A.Th(k),B.Th(k));
end
fprintf('  %6s %16.5f %16.5f\n','somme',sum(A.Th),sum(B.Th));

%% ---- 3. les deux colonnes contre la reference EF ---------------------
fprintf('\n  ---- 3. ECARTS A LA REFERENCE (famille transitoire) ----\n');
fprintf('  %-26s %12s %12s %12s\n','grandeur','vrill. ON','vrill. OFF','EF');
e3=@(l,a,b,f)fprintf('  %-26s %11.1f %% %11.1f %% %12.4f\n', ...
    l,100*(a-f)/f,100*(b-f)/f,f);
e3('couple en charge [N.m]', A.Tem, B.Tem, 121.63);
e3('courant I1 [A]',         A.I1,  B.I1,  19.73);
e3('courant de barre [A]',   A.Ibar,B.Ibar,324.7);
e3('courant a vide [A]',     A.I0,  B.I0,  8.49);
e3('f.e.m. entrefer [V]',    A.E1,  B.E1,  382.1);
e3('couple au calage [N.m]', A.Tb,  B.Tb,  104.31);
e3('I1 au calage [A]',       A.I1b, B.I1b, 108.21);
e3('I barre au calage [A]',  A.Ibarb,B.Ibarb,1929);

%% ---- 4. le residu irreductible du fondamental ------------------------
fprintf('\n  ---- 4. RESIDU DU FONDAMENTAL ----\n');
fprintf('  ksq(nu=1) = %.6f dans LES DEUX configurations (cage.m:61)\n',A.ksq);
fprintf('  R''_r et X''_r restent divises par ksq^2 = %.6f, soit +%.2f %%\n', ...
    A.ksq^2,100*(1/A.ksq^2-1));
fprintf('  que la reference, non vrillee, n''a pas. IRREDUCTIBLE sans\n');
fprintf('  modifier cage.m. A DECLARER dans la note de tab:im_ec.\n');

%% ---- 5. verdict ------------------------------------------------------
fprintf('\n  ---- 5. VERDICT ----\n');
dT = d(A.Tem,B.Tem);
fprintf('  Neutraliser le vrillage harmonique deplace le couple net de\n');
fprintf('  %.2f %% et l''ecart a la reference de %.1f %% a %.1f %%.\n', ...
    dT,100*(A.Tem-121.63)/121.63,100*(B.Tem-121.63)/121.63);
fprintf('  La comparaison LOYALE est la colonne OFF : la reference est une\n');
fprintf('  tranche droite. Si l''ecart y est plus grand, c''est le resultat,\n');
fprintf('  et il se publie tel quel.\n');

save('B8_skewoff.mat','A','B');
fprintf('\n  duree %.0f s\n=== B8 termine ===\n',toc(t0));
diary off;
