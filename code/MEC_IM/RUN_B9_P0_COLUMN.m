%% RUN_B9_P0_COLUMN - la colonne P0 de tab:balance, en une execution
%
%  POURQUOI CE BLOC. tab:balance de l'Article II publie quatre lignes en
%  DEUX colonnes -- base P0 (publiee) et base P1 (convergee) :
%
%      couple au calage             -4,6 %   ->  -1,6 %
%      B_t rms en charge (dentaire) -12,6 %  ->  -9,3 %
%      B_g1 a vide                    --     ->  +0,3 %
%      B_g1 en charge               -1,2 %   ->  -2,4 %
%
%  La colonne P1 sort de B1_im_p1_out.txt. La colonne P0 ne sort de RIEN :
%  elle est reprise du manuscrit v2, donc D'UN TABLEAU. La regle
%  d'engagement n. 1 l'interdit, et le controle le confirme -- le manuscrit
%  donne 0,114 T en charge, soit -13,0 %, la ou tab:balance ecrit -12,6 %.
%
%  Ce script produit les quatre grandeurs en base P0, dans la configuration
%  EXACTE de la colonne publiee, et les quatre memes en base P1 comme
%  controle croise de B1.
%
%  CONFIGURATION P0 : pavage nT=6, nO=2, N_h = 3088, base 'p0' -- c'est la
%  configuration du manuscrit, celle que RUN_B3_PRIX emploie deja pour sa
%  ligne "P0 (publie)". Elle donne X_m0 = 64,78 ohm.
%  CONFIGURATION P1 : pavage nT=17, nO=4, N_h = 8192, base 'p1', X_m0 =
%  61,015 ohm -- celle de B1 et B2.
%
%  REFERENCES EF : famille TRANSITOIRE (decision 1 du 6 aout 2026).
%    calage 104,31 N.m (moyenne du regime etabli du transitoire rotor
%    bloque, 2001 points) -- et NON 98,98 N.m, qui est un point isole d'un
%    balayage dont la dispersion locale atteint 16 % pour s >= 0,90.
%    B_g1 a vide 0,942 T | B_g1 en charge 0,920 T | B_t rms en charge 0,131 T
%
%  Sortie : B9_p0_column_out.txt
clear; clc; t0=tic;
diary('B9_p0_column_out.txt'); diary on;

s_ch = 0.0188;
CFG = { struct('lab','P0 (publiee)','nT', 6,'nO',2,'Nh',3088,'basis','p0'), ...
        struct('lab','P1 (convergee)','nT',17,'nO',4,'Nh',8192,'basis','p1') };

fprintf('=== B9 : colonne P0 de tab:balance, et controle croise P1 ===\n');
fprintf('  references EF, famille TRANSITOIRE (decision 1) :\n');
fprintf('    couple au calage  104.31 N.m   (transitoire rotor bloque)\n');
fprintf('    B_g1 a vide         0.942 T\n');
fprintf('    B_g1 en charge      0.920 T\n');
fprintf('    B_t rms en charge   0.131 T\n\n');

Q = cell(2,1);
for c = 1:2
    C = CFG{c};
    M = mec.machine_18_5kW(); ctx = mec.build_context(M); G = ctx.G;
    A = mec.airgap_dtn_tooth(M,G,0,C.nT,C.nO,C.Nh,C.basis);
    ctx.AG = A; ctx.Xm0 = mec.magnetizing(ctx,0.2).Xm;
    fprintf('  --- %s : pavage nT=%d nO=%d | N_h = %d | base %s\n', ...
        C.lab,C.nT,C.nO,C.Nh,upper(C.basis));
    fprintf('      operateur %dx%d | X_m0 = %.3f ohm\n',size(A.Y,1),size(A.Y,2),ctx.Xm0);

    r  = mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
    r0 = mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
    rb = mec.equivalent_circuit(ctx,1.0 ,ctx.Xm0);

    % --- champs a mi-entrefer, meme protocole que B1 §4 -----------------
    Rm = 0.5*(G.Rs+G.Rr); thq = linspace(0,2*pi,2001); thq(end) = [];
    Rnl = mec.magnetizing(ctx,r0.Im);
    [Br0,~] = A.field(Rnl.S.Usurf,Rm,thq); Br0 = Br0(:).';

    p1 = angle(r.I1c);
    i3 = sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
    Fs = ctx.W.slotMMF(i3);
    p2 = angle(r.I2c); tb = 2*pi*(0:M.Nr-1)/M.Nr;
    Fu = cumsum(cos(M.p*tb-p2).'); Fu = Fu-mean(Fu);
    c1 = (2/M.Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
    Fr = -Fu*((3/2)*(4/pi)*(ctx.W.kw1*ctx.W.Nph/(2*M.p))*sqrt(2)*abs(r.I2c)/abs(c1));
    SL = mec.solve_network(ctx.net,G,ctx.BH,A,Fs,Fr,M.opt);
    [Br1,Bt1] = A.field(SL.Usurf,Rm,thq); Br1 = Br1(:).'; Bt1 = Bt1(:).';

    ff = @(y) abs((2/numel(y))*sum(y.*exp(-1j*M.p*thq)));
    S = struct('Xm0',ctx.Xm0,'Tb',rb.Tem, ...
               'Bt_ch',sqrt(mean(Bt1.^2)),'Bg1_0',ff(Br0),'Bg1_ch',ff(Br1));
    Q{c} = S;
    fprintf('      fait (%.0f s cumulees)\n',toc(t0));
end
fprintf('\n');

P0 = Q{1}; P1 = Q{2};

%% ---- tab:balance, les deux colonnes ---------------------------------
fprintf('  ---- tab:balance : LES QUATRE LIGNES, DEUX BASES ----\n');
fprintf('  %-30s %12s %12s %10s\n','grandeur','P0','P1','EF');
v=@(l,a,b,f)fprintf('  %-30s %12.4f %12.4f %10.4f\n',l,a,b,f);
v('couple au calage [N.m]',   P0.Tb,    P1.Tb,    104.31);
v('B_t rms en charge [T]',    P0.Bt_ch, P1.Bt_ch, 0.131);
v('B_g1 a vide [T]',          P0.Bg1_0, P1.Bg1_0, 0.942);
v('B_g1 en charge [T]',       P0.Bg1_ch,P1.Bg1_ch,0.920);

fprintf('\n  ---- ECARTS, FORMES SUR LES VALEURS PLEINES ----\n');
fprintf('  %-30s %12s %12s\n','grandeur','P0','P1');
e=@(l,a,b,f)fprintf('  %-30s %11.2f %% %11.2f %%\n',l,100*(a-f)/f,100*(b-f)/f);
e('couple au calage',   P0.Tb,    P1.Tb,    104.31);
e('B_t rms en charge',  P0.Bt_ch, P1.Bt_ch, 0.131);
e('B_g1 a vide',        P0.Bg1_0, P1.Bg1_0, 0.942);
e('B_g1 en charge',     P0.Bg1_ch,P1.Bg1_ch,0.920);

fprintf('\n  X_m0 : P0 %.3f ohm | P1 %.3f ohm  (ecart %.2f %%)\n', ...
    P0.Xm0,P1.Xm0,100*(P1.Xm0-P0.Xm0)/P0.Xm0);

%% ---- controles -------------------------------------------------------
fprintf('\n  ---- CONTROLES ----\n');
fprintf('  1. La colonne P1 doit reproduire B1_im_p1_out.txt :\n');
fprintf('     couple calage 102.6023 | B_t charge 0.1188 | B_g1 vide 0.9452\n');
fprintf('     | B_g1 charge 0.8983.  Ecarts observes :\n');
ref = [102.6023 0.1188 0.9452 0.8983];
got = [P1.Tb P1.Bt_ch P1.Bg1_0 P1.Bg1_ch];
nm  = {'couple calage','B_t charge','B_g1 vide','B_g1 charge'};
for k=1:4
    fprintf('       %-14s %12.4f  vs %10.4f  ->  %+.3f %%\n', ...
        nm{k},got(k),ref(k),100*(got(k)-ref(k))/ref(k));
end
fprintf('     Un ecart superieur a 0,05 %% signale une difference de\n');
fprintf('     configuration entre ce script et B1 : la trouver avant de\n');
fprintf('     publier la colonne P0.\n');
fprintf('  2. X_m0 en P0 doit valoir 64,78 ohm (valeur du manuscrit).\n');
fprintf('     Sinon la configuration P0 n''est pas celle qui a ete publiee.\n');

save('B9_p0_column.mat','P0','P1');
fprintf('\n  duree %.0f s\n=== B9 termine ===\n',toc(t0));
diary off;
