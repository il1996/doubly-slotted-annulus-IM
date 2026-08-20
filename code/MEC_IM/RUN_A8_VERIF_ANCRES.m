%% RUN_A8_VERIF_ANCRES - controle des ancres de mec.ansys_ref avant correction
%
%  La note du 6 aout signale que R.anchor.I0 = 8.32 A est faux et que la
%  valeur mesuree est 8,49 A. On ne reprend pas ce chiffre : on le mesure.
%  Meme chose pour les autres ancres, verifiables ou non.
%
%  Sortie : A8_verif_ancres_out.txt
clear; clc; t0=tic;
diary('A8_verif_ancres_out.txt'); diary on;
M=mec.machine_18_5kW();
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
rd  =@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
rmsw=@(A,c,t0)sqrt(mean(A(A(:,1)>=t0,c).^2,'omitnan'));
avgw=@(A,c,t0)mean(A(A(:,1)>=t0,c),'omitnan');
R=mec.ansys_ref();

fprintf('=== A8-0 : les ancres de mec.ansys_ref, mesurees ===\n');
fprintf('  reference : %s\n\n',ROOT);

%% ---- I0 : courant a vide ---------------------------------------------
D0=fullfile(ROOT,'transitoire','a vide');
cur=rd(fullfile(D0,'Winding Plot 4.tab'));
t0h=0.5*max(cur(:,1));
%  Le dossier "a vide" ne porte pas de "Speed Plot 1.tab" -- la vitesse y
%  est dans "la vitesse en fonction du temps.tab". On la lit si elle est
%  lisible, sinon on declare le glissement comme non verifie ici.
n0=NaN; s0=NaN;
try
    spd=rd(fullfile(D0,'la vitesse en fonction du temps.tab'));
    n0=avgw(spd,2,t0h); s0=(M.ns*60-n0)/(M.ns*60);
catch
end
fprintf('  ---- I0, courant a vide ----\n');
fprintf('    fichier : transitoire\\a vide\\Winding Plot 4.tab (%d points)\n',size(cur,1));
if isfinite(n0)
    fprintf('    vitesse etablie %.3f tr/min  ->  s = %.3e\n',n0,s0);
else
    fprintf('    vitesse : non relue (fichier de vitesse au format non tabule)\n');
end
ph={'phase A','phase B','phase C'}; Iph=nan(1,3);
for k=1:3
    Iph(k)=rmsw(cur,k+1,t0h);
    fprintf('    %-8s (col %d) : %.4f A\n',ph{k},k+1,Iph(k));
end
I0m=mean(Iph);
fprintf('    MOYENNE DES TROIS PHASES : %.4f A\n',I0m);
fprintf('    ancre R.anchor.I0        : %.4f A  ->  ecart %+.2f %%\n', ...
    R.anchor.I0,100*(R.anchor.I0/I0m-1));
fprintf('    valeur employee par B1   : 8.4900 A  ->  ecart %+.2f %%\n', ...
    100*(8.49/I0m-1));

%% ---- Rs : resistance statorique -------------------------------------
fprintf('\n  ---- Rs, resistance statorique ----\n');
ctx=mec.build_context(M);
fprintf('    ancre R.anchor.Rs        : %.4f ohm\n',R.anchor.Rs);
fprintf('    valeur employee par B1   : 0.4450 ohm (colonne EF de tab:im_ec)\n');
fprintf('    valeur du modele ctx.Rs  : %.4f ohm\n',ctx.Rs);
fprintf('    L''ancre est plus petite d''un facteur %.2f que les deux autres.\n', ...
    0.4450/R.anchor.Rs);
fprintf('    Une resistance de phase de %.3f ohm sur une machine 690 V /\n',R.anchor.Rs);
fprintf('    18,5 kW donnerait des pertes Joule statoriques de %.0f W a\n', ...
    3*R.anchor.Rs*19.73^2);
fprintf('    19,73 A, contre %.0f W avec 0,4450 ohm. INVRAISEMBLABLE.\n', ...
    3*0.4450*19.73^2);

%% ---- Bg1 : fondamental d'entrefer ------------------------------------
fprintf('\n  ---- Bg1, fondamental a mi-entrefer ----\n');
fprintf('    ancre R.anchor.Bg1       : %.4f T\n',R.anchor.Bg1);
fprintf('    valeur employee par B1   : 0.9420 T (colonne EF)\n');
fprintf('    ecart entre les deux     : %+.1f %%\n',100*(R.anchor.Bg1/0.942-1));
fprintf('    B1 §4 compare son 0,9452 T a 0,942 : +0,3 %%. Contre l''ancre\n');
fprintf('    ce serait %+.1f %%. Les deux ne peuvent pas etre la meme grandeur.\n', ...
    100*(0.9452/R.anchor.Bg1-1));

%% ---- ancres verifiables : couple de decrochage -----------------------
fprintf('\n  ---- ancres verifiables sur la caracteristique ----\n');
Dc=fullfile(ROOT,'caratéristique en fonction glissement');
tc=rd(fullfile(Dc,'Torque Plot 2.tab'));
[Tmx,im]=max(tc(:,2));
fprintf('    fichier Torque Plot 2.tab : %d points, s de %.4f a %.4f\n', ...
    size(tc,1),min(tc(:,1)),max(tc(:,1)));
fprintf('    maximum lu               : %.2f N.m a s = %.4f\n',Tmx,tc(im,1));
fprintf('    ancre T_break / s_break  : %.2f N.m a s = %.4f\n', ...
    R.anchor.T_break,R.anchor.s_break);
fprintf('    -> ecart %+.2f %% sur le couple, %+.4f sur le glissement\n', ...
    100*(R.anchor.T_break/Tmx-1),R.anchor.s_break-tc(im,1));

%% ---- dispersion locale de la caracteristique ------------------------
fprintf('\n  ---- dispersion locale, par bande de glissement ----\n');
bd=[0.005 0.05; 0.05 0.15; 0.135 0.21; 0.40 0.90; 0.90 1.00];
fprintf('  %14s %8s %12s %12s\n','bande en s','points','bruit rms','etendue');
for k=1:size(bd,1)
    m=tc(:,1)>=bd(k,1) & tc(:,1)<=bd(k,2);
    x=tc(m,1); y=tc(m,2);
    if numel(y)<3, continue; end
    %  residu a l'interpolation des deux voisins
    r=y(2:end-1)-0.5*(y(1:end-2)+y(3:end));
    fprintf('  %6.3f - %5.3f %8d %11.2f %% %6.1f - %6.1f\n', ...
        bd(k,1),bd(k,2),numel(y),100*rms(r)/mean(abs(y)),min(y),max(y));
end
m9=tc(:,1)>=0.90;
fprintf('\n    bande s >= 0,90 : %d points, moyenne %.2f, ecart-type %.2f N.m (%.1f %%)\n', ...
    sum(m9),mean(tc(m9,2)),std(tc(m9,2)),100*std(tc(m9,2))/mean(tc(m9,2)));
fprintf('    min / max : %.2f / %.2f N.m\n',min(tc(m9,2)),max(tc(m9,2)));

%% ---- verdict ----------------------------------------------------------
fprintf('\n  ---- VERDICT ----\n');
fprintf('    I0  : ancre %.2f contre %.4f mesure  ->  FAUSSE (%+.1f %%)\n', ...
    R.anchor.I0,I0m,100*(R.anchor.I0/I0m-1));
fprintf('    Rs  : ancre %.3f contre 0,4450 employe -> FAUSSE (facteur %.1f)\n', ...
    R.anchor.Rs,0.4450/R.anchor.Rs);
fprintf('    Bg1 : ancre %.2f contre 0,942 employe  -> FAUSSE (%+.1f %%)\n', ...
    R.anchor.Bg1,100*(R.anchor.Bg1/0.942-1));
fprintf('    T_break, s_break : retrouves dans le fichier -> VERIFIEES\n');
fprintf('    Xm = 46, Rfe = 1740 : employees par B1, non re-mesurables ici\n');
fprintf(['\n    CONSEQUENCE. La sous-structure anchor melange des valeurs\n' ...
         '    verifiees et des valeurs fausses. Elle ne doit alimenter AUCUN\n' ...
         '    chiffre des deux articles. Les trois fausses sont marquees dans\n' ...
         '    ansys_ref.m plutot que supprimees : les supprimer casserait les\n' ...
         '    scripts de diagnostic qui les lisent, et la marque est plus\n' ...
         '    instructive que l''absence.\n']);

save('A8_verif_ancres.mat','Iph','I0m','tc');
fprintf('\n  duree %.0f s\n=== A8-0 termine ===\n',toc(t0));
diary off;
