%% RUN_R8_TABLE2 - R8 : regenerer la Table 2 de l'Article II sans sa note
%
%  CONSTAT DU RAPPORTEUR. "Le tableau principal du manuscrit est dans une
%  configuration differente de celle du reste du manuscrit, et la difference
%  est absorbee par une note."
%
%  LOCALISATION. Table 2 = tab:price (ArticleII_Carter_IM.tex:783-834). Ses
%  trois colonnes sont calculees VRILLAGE HARMONIQUE ACTIF, contrairement aux
%  Tables 3 (tab:im_tests) et 4 (tab:im_ec). La note absorbante occupe les
%  lignes 808-818 : elle annonce +21,3 % et -10,3 % ici contre +21,6 % et
%  -10,6 % ailleurs, "la difference etant le cage harmonique, non la
%  fermeture".
%
%  CE QU'ON FAIT. On regenere les TROIS fermetures dans la configuration
%  LOYALE -- M.opt.skew_harm = 0, celle du reste du manuscrit -- et on produit
%  aussi la configuration ACTIVE, pour que l'affirmation de la note soit
%  VERIFIEE au lieu d'etre admise. La reference EF est une tranche DROITE
%  (IM_18kW_690V.aedt : NumberOfSlices = 1), donc la colonne loyale est la
%  seule comparable.
%
%  REFERENCES. Lues des .tab quand elles le sont ; declarees avec leur statut
%  quand elles ne le sont pas. mec.ansys_ref declare explicitement X_m = 46
%  ohm "NON RE-MESURABLE depuis les .tab disponibles" : on ne fait donc pas
%  semblant de la relire.
%
%  GARDE. Apres regeneration, les deux ecarts de la colonne convergee doivent
%  coincider avec ceux que le manuscrit cite AILLEURS -- relus du .tex, non
%  transcrits. S'ils coincident, la note n'a plus d'objet ; sinon la
%  regeneration n'a pas ete faite dans la bonne configuration.
clear; clc; t0=tic;
diary('R8_table2_out.txt'); diary on;
nT=17; nO=4; Nh=8192; s_ch=0.0188;
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
TEX=fullfile('..','article','ArticleII_Carter_IM.tex');
rd=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
rmsw=@(A,c,tw)sqrt(mean(A(A(:,1)>=tw,c).^2,'omitnan'));
avgw=@(A,c,tw)mean(A(A(:,1)>=tw,c),'omitnan');

fprintf('=== R8 : Table 2 de l''Article II, configuration loyale ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    machine      : MAS 48/44, 18,5 kW, 690 V (mec.machine_18_5kW)\n');
fprintf('    pavage       : n_T = %d, n_O = %d ; N_h = %d\n',nT,nO,Nh);
fprintf('    point        : s = %.4f\n',s_ch);
fprintf('    vrillage     : harmonique NEUTRALISE (M.opt.skew_harm = 0)\n');
fprintf('                   ET actif, pour verifier la note\n');
fprintf('    reference EF : %s\n',ROOT);

%% ---- references EF, relues ou declarees -------------------------------
D0=fullfile(ROOT,'transitoire','a vide');
DC=fullfile(ROOT,'transitoire','en charge');
tw=1.0;
I0F=rmsw(rd(fullfile(D0,'Winding Plot 4.tab')),2,tw);
E1F=rmsw(rd(fullfile(D0,'Winding Plot 2.tab')),2,tw);
TF_tr=avgw(rd(fullfile(DC,'Plot 1.tab')),3,tw);
fprintf('\n  ---- REFERENCES ----\n');
fprintf('    I0  = %.4f A     RELUE  (transitoire\\a vide\\Winding Plot 4.tab, t > %.1f s)\n',I0F,tw);
fprintf('    E1  = %.2f V   RELUE  (transitoire\\a vide\\Winding Plot 2.tab)\n',E1F);
fprintf('    T   = %.2f N.m RELUE  (transitoire\\en charge\\Plot 1.tab)\n',TF_tr);
XmF=46.0;
fprintf('    X_m = %.1f ohm   *** NON RELUE *** : mec.ansys_ref la declare\n',XmF);
fprintf('          "non re-mesurable depuis les .tab disponibles, statut NON VERIFIE".\n');
fprintf('          Elle est employee telle quelle et son statut est declare.\n');

%% ---- six contextes : trois fermetures x deux reglages de vrillage -----
nm={'Carter','piecewise const. (publ.)','hat (converged)'};
V=nan(3,5,2);            % (fermeture, grandeur, reglage) ; 1 = OFF, 2 = ON
for c=1:2
    for k=1:3
        M=mec.machine_18_5kW();
        if c==1, M.opt.skew_harm=0; end
        ctx=mec.build_context(M); G=ctx.G;
        switch k
            case 1, ctx.AG=ctx.AGcarter;
            case 2, ctx.AG=mec.airgap_dtn_tooth(M,G,0,6,2,3088,'p0');
            case 3, ctx.AG=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
        end
        if c==1, ctx.M.opt.skew_harm=0; end     % la cage lit ctx.M
        ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
        r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
        rn=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
        V(k,:,c)=[ctx.Xm0, rn.Xm, r0.I1, r0.E1, rn.Tem];
    end
    fprintf('  reglage %d/2 : fait (%.0f s)\n',c,toc(t0));
end

%% ---- le tableau, configuration loyale ---------------------------------
lbl={'X_m0 (ohm)','X_m charge (ohm)','I0 (A)','E1 (V)','couple (N.m)'};
fprintf('\n  ---- TABLE 2 REGENEREE : VRILLAGE HARMONIQUE NEUTRALISE ----\n');
fprintf('  %-26s %11s %13s %10s %10s %12s\n','fermeture',lbl{:});
for k=1:3
    fprintf('  %-26s %11.3f %13.3f %10.3f %10.2f %12.3f\n',nm{k},V(k,:,1));
end
fprintf('  %-26s %11s %13.1f %10.4f %10.2f %12.2f\n','REFERENCE EF','---',XmF,I0F,E1F,TF_tr);
fprintf('\n  ecarts a la reference :\n');
fprintf('  %-26s %13s %10s %10s %12s\n','fermeture','X_m charge','I0','E1','couple');
dev=nan(3,4);
for k=1:3
    dev(k,:)=[100*(V(k,2,1)-XmF)/XmF, 100*(V(k,3,1)-I0F)/I0F, ...
              100*(V(k,4,1)-E1F)/E1F, 100*(V(k,5,1)-TF_tr)/TF_tr];
    fprintf('  %-26s %12.1f %% %9.1f %% %9.1f %% %11.1f %%\n',nm{k},dev(k,:));
end

%% ---- la note, verifiee au lieu d'etre admise --------------------------
fprintf('\n  ---- LA NOTE DE LA TABLE 2, VERIFIEE ----\n');
fprintf('  Elle affirme que l''ecart entre les deux reglages est "le cage\n');
fprintf('  harmonique, non la fermeture". Si c''est vrai, le decalage doit\n');
fprintf('  etre le MEME pour les trois fermetures.\n\n');
fprintf('  %-26s %14s %14s %12s\n','fermeture','I0 vrill. ON','I0 vrill. OFF','decalage');
for k=1:3
    dON=100*(V(k,3,2)-I0F)/I0F; dOFF=100*(V(k,3,1)-I0F)/I0F;
    fprintf('  %-26s %12.2f %% %12.2f %% %10.3f pt\n',nm{k},dON,dOFF,dOFF-dON);
end
fprintf('\n  %-26s %14s %14s %12s\n','fermeture','Xm vrill. ON','Xm vrill. OFF','decalage');
for k=1:3
    dON=100*(V(k,2,2)-XmF)/XmF; dOFF=100*(V(k,2,1)-XmF)/XmF;
    fprintf('  %-26s %12.2f %% %12.2f %% %10.3f pt\n',nm{k},dON,dOFF,dOFF-dON);
end

%% ---- GARDE ------------------------------------------------------------
fprintf('\n  ---- GARDE : coincidence avec les valeurs citees ailleurs ----\n');
tx=fileread(TEX);
t1=regexp(tx,'no-load current is high by\s*\\SI\{([\d.]+)\}\{\\percent\}','tokens','once');
if isempty(t1)
    t1=regexp(tx,'\\SI\{(21\.\d+)\}\{\\percent\}','tokens','once');
end
t2=regexp(tx,'magnetising reactance is low by\s*\\SI\{([\d.]+)\}\{\\percent\}','tokens','once');
if isempty(t2)
    t2=regexp(tx,'\\SI\{-?(10\.\d+)\}\{\\percent\}','tokens','once');
end
cI=str2double(t1{1}); cX=str2double(t2{1});
gI=dev(3,2); gX=dev(3,1);
fprintf('    cite ailleurs dans le .tex : I0 %+.1f %% | X_m -%.1f %%\n',cI,cX);
fprintf('    Table 2 regeneree          : I0 %+.1f %% | X_m %+.1f %%\n',gI,gX);
fprintf('    ecart                      : %.3f pt | %.3f pt\n',abs(gI-cI),abs(abs(gX)-cX));
G=abs(gI-cI)<0.1 && abs(abs(gX)-cX)<0.1;
if G
    fprintf('    GARDE PASSEE : les chiffres coincident SANS note explicative.\n');
    fprintf('    La note des lignes 808-818 du .tex n''a plus d''objet.\n');
else
    fprintf('    GARDE ECHOUEE : la regeneration n''est pas dans la bonne\n');
    fprintf('    configuration, ou les valeurs citees ailleurs sont autres.\n');
end

%% ---- panneau LaTeX ----------------------------------------------------
fprintf('\n  ---- PANNEAU LaTeX POUR LA TABLE 2 ----\n\n');
for k=1:3
    if k==3
        fprintf('\\textbf{hat (converged)}& \\textbf{\\num{%.3f}} & \\textbf{\\num{%.3f}} & \\textbf{\\num{%.3f}} & \\textbf{\\num{%.2f}} & \\textbf{\\num{%.3f}}\\\\\n',V(k,:,1));
    else
        fprintf('%-23s& \\num{%.3f} & \\num{%.3f} & \\num{%.3f} & \\num{%.2f} & \\num{%.3f}\\\\\n',nm{k},V(k,:,1));
    end
end
fprintf('\\midrule\nFEA reference          & --- & \\num{%.1f} & \\num{%.4f} & \\num{%.2f} & \\num{%.2f}\\\\\n',XmF,I0F,E1F,TF_tr);
fprintf('\\midrule\n\\multicolumn{6}{l}{\\emph{deviation from the reference}}\\\\\n');
for k=1:3
    fprintf('%-23s& --- & $%+.1f\\%%$ & $%+.1f\\%%$ & $%+.1f\\%%$ & $%+.1f\\%%$\\\\\n',nm{k},dev(k,:));
end

save('R8_table2.mat','V','dev','nm','I0F','E1F','TF_tr','XmF','G');
fprintf('\n  duree %.0f s\n=== R8 termine ===\n',toc(t0));
diary off;
