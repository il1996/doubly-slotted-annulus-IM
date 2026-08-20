%% RUN_B3_PRIX - le prix de la base P1, chiffre et non estime
%
%  B3. Le rapport estimait +23 % sur le courant a vide et ~-11,5 % sur X_m
%  sature. Il faut les valeurs CALCULEES.
%
%  On construit donc le contexte en base P1 au N_h converge (B2 : k_C stable
%  a 0,62 %, X_m encoche 61,015 ohm) et on relit :
%    - courant a vide et f.e.m. d'entrefer
%    - X_m sature au point nominal s = 0,0188
%    - les ecarts a la reference EF
%  Comparaison a TROIS colonnes : Carter (etat anterieur), P0 (etat publie),
%  P1 (etat converge).
clear; clc; t0=tic;
diary('B3_prix_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx0=mec.build_context(M); G=ctx0.G;
nT=17; nO=4; Nh=8192; s_ch=0.0188;
%  References EF
I0F=8.49; E1F=382.1; XmF=46.0;

fprintf('=== B3 : le prix de la base P1, chiffre ===\n');
fprintf('  pavage nT=%d nO=%d | N_h = %d | point nominal s = %.4f\n\n',nT,nO,Nh,s_ch);

%% ---- trois contextes -------------------------------------------------
C=struct(); nm={'CARTER','P0 (publie)','P1 (converge)'};
%  (1) Carter : modele de comparaison
C(1).ctx=ctx0; C(1).ctx.AG=ctx0.AGcarter;
C(1).ctx.Xm0=mec.magnetizing(C(1).ctx,0.2).Xm;
%  (2) P0 au N_h du manuscrit (3088, pavage 6/2) -- etat publie
A0=mec.airgap_dtn_tooth(M,G,0,6,2,3088,'p0');
C(2).ctx=ctx0; C(2).ctx.AG=A0; C(2).ctx.Xm0=mec.magnetizing(C(2).ctx,0.2).Xm;
%  (3) P1 converge
A1=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
C(3).ctx=ctx0; C(3).ctx.AG=A1; C(3).ctx.Xm0=mec.magnetizing(C(3).ctx,0.2).Xm;

%% ---- grandeurs -------------------------------------------------------
V=nan(3,5);
for k=1:3
    r0=mec.equivalent_circuit(C(k).ctx,1e-4,C(k).ctx.Xm0);
    rn=mec.equivalent_circuit(C(k).ctx,s_ch,C(k).ctx.Xm0);
    V(k,:)=[C(k).ctx.Xm0, rn.Xm, r0.I1, r0.E1, rn.Tem];
end

fprintf('  %-16s %11s %11s %10s %10s %11s\n', ...
    'modele','Xm0 (ohm)','Xm charge','I0 (A)','E1 (V)','couple (N.m)');
for k=1:3
    fprintf('  %-16s %11.3f %11.3f %10.3f %10.2f %11.3f\n',nm{k},V(k,:));
end
fprintf('  %-16s %11s %11.1f %10.2f %10.1f %11s\n','REFERENCE EF','--',XmF,I0F,E1F,'121.63');

fprintf('\n  --- ecarts a la reference EF ---\n');
fprintf('  %-16s %11s %10s %10s %11s\n','modele','Xm charge','I0','E1','couple');
TF=121.63;
for k=1:3
    fprintf('  %-16s %10.1f %% %9.1f %% %9.1f %% %10.1f %%\n',nm{k}, ...
        100*(V(k,2)-XmF)/XmF,100*(V(k,3)-I0F)/I0F, ...
        100*(V(k,4)-E1F)/E1F,100*(V(k,5)-TF)/TF);
end

%% ---- le prix, isole ---------------------------------------------------
fprintf('\n  --- LE PRIX : P0 publie -> P1 converge ---\n');
lb={'Xm0 non sature','Xm sature charge','courant a vide','f.e.m. entrefer','couple'};
for j=1:5
    fprintf('  %-18s %10.3f -> %10.3f   %+7.2f %%\n', ...
        lb{j},V(2,j),V(3,j),100*(V(3,j)-V(2,j))/V(2,j));
end
fprintf(['\n  LECTURE. La base P1 n''ameliore PAS l''accord : elle le degrade\n' ...
         '  sur le courant a vide. Ce qu''elle apporte est que la valeur\n' ...
         '  EXISTE -- en P0 la serie diverge (T17, forme fermee), donc\n' ...
         '  X_m n''y est pas defini. On echange un chiffre flatteur mais\n' ...
         '  non convergent contre un chiffre defavorable mais VRAI.\n']);
save('B3_prix.mat','V','nm');
fprintf('\n  duree %.0f s\n=== B3 termine ===\n',toc(t0));
diary off;
