%% SMOKE_CARTES - controle des libelles ANGLAIS des cartes 2D
%
%  Verifie (a) que la surcharge des libelles depuis l'appelant fonctionne,
%  (b) que les chaines arrivent intactes dans les objets graphiques,
%  (c) que les PNG s'ecrivent. Reprend le bloc de RUN_ARTICLE section IV.
clear; clc; close all;
M=mec.machine_18_5kW(); ctx=mec.build_context(M);
G=ctx.G; W=ctx.W; BH=ctx.BH;
OUT=fullfile(tempdir,'smoke_cartes'); if ~isfolder(OUT), mkdir(OUT); end

r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
q1=angle(r0.I1c); q2=angle(r0.I2c);
i30=sqrt(2)*r0.I1*[cos(q1);cos(q1-2*pi/3);cos(q1+2*pi/3)];
Ib=r0.I2*(2*M.m*W.kw1*W.Nph)/M.Nr;
ib0=sqrt(2)*Ib*cos(M.p*(2*pi*(0:M.Nr-1)'/M.Nr)-q2);

me0=mec.mesh_refined(M,G,6,i30,0,3,ib0); Se0=mec.solve_mesh(me0,BH,M.opt);

figEN=struct( ...
    'titleB','Flux density |B| and flux lines (actual slot geometry)', ...
    'titleJ','Instantaneous current density (bars + two stator layers)', ...
    'labB','|B| [T]', ...
    'labJ','J [A/m^2]', ...
    'nameB','Flux density map |B|', ...
    'nameJ','Current density map J');
o0=figEN; o0.prefix='No load: '; o0.Jmax=3e6;
mec.draw_cross_section(me0,Se0,M,G,W,i30,ib0,o0);

figs=get(groot,'Children');
fprintf('\n--- libelles effectivement poses sur les figures ---\n');
for k=numel(figs):-1:1
    ax=findobj(figs(k),'Type','axes');
    cb=findobj(figs(k),'Type','colorbar');
    fprintf('fenetre : %s\n',figs(k).Name);
    if ~isempty(ax),  fprintf('  titre  : %s\n',ax(1).Title.String); end
    if ~isempty(cb),  fprintf('  echelle: %s\n',cb(1).Label.String); end
end
saveas(figs(2),fullfile(OUT,'carte_B.png'));
saveas(figs(1),fullfile(OUT,'carte_J.png'));
d=dir(fullfile(OUT,'*.png'));
fprintf('\nPNG ecrits dans %s :\n',OUT);
for k=1:numel(d), fprintf('  %-12s %6.0f ko\n',d(k).name,d(k).bytes/1024); end
fprintf('\n=== SMOKE_CARTES termine ===\n');
