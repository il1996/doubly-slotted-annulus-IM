%% RUN_B7_PROBES - bloc B7 (SPEC_CLAUDE_CODE_v3 §8, priorite 3)
%
%  OBJET. B1 a laisse ouvertes les sondes locales de la Table 17 : elles
%  exigent mec.mesh_refined, "chaine distincte du reseau de performance".
%  B7 les produit -- et etablit en quoi cette chaine est distincte, ce qui
%  est le point le plus important du bloc.
%
%  CE QUE LA LECTURE DU CODE ETABLIT AVANT TOUT CALCUL
%
%   (1) L'ENTREFER N'EST PAS LE MEME OPERATEUR. Le reseau de performance
%       passe par mec.airgap_dtn_tooth (pavage nT/nO, condensation de Schur
%       sur 92 noeuds dentaires). La chaine de sondes passe par
%       mec.airgap_fourier (mesh_refined.m:271), sur la grille de surface
%       remappee aux largeurs physiques. Deux operateurs, deux chaines.
%
%   (2) LA BASE EST P0, ET LE CALLEUR NE PEUT PAS EN CHANGER.
%       airgap_fourier.m:63 : basis = 'p0' par defaut. mesh_refined
%       l'appelle SANS argument de base. Il n'existe aucun crochet
%       M.opt pour la changer. Les sondes sont donc calculees dans la base
%       dont la section 3 de l'Article I demontre qu'elle n'est pas
%       conforme -- et sur des grandeurs LOCALES, celles qui y sont le
%       plus sensibles.
%
%   (3) LA TRONCATURE VAUT 100. mesh_refined.m:258 : NhF = 100, reglable
%       par M.opt.gap_fourier_Nh. Le reseau de performance, lui, tourne a
%       N_h = 8192 (B1/B2). Un facteur 82.
%
%  CE BLOC NE MODIFIE PAS LE CODE DE PRODUCTION. Il balaie la seule
%  poignee que le caller possede -- la troncature -- et mesure si les
%  sondes derivent. C'est le test de la section 3 applique aux grandeurs
%  locales.
%
%  CONFIGURATION DECLAREE
%    machine    : MAS 48/44, 18,5 kW
%    maillage   : mec.mesh_refined(M,G,nc=6,i3,theta=0,nr=3,ibar)
%    entrefer   : mec.airgap_fourier, base P0 (imposee), N_h BALAYE
%    circuit    : mec.equivalent_circuit sur ctx en base P1, N_h = 8192,
%                 pavage nT=17 nO=4 (identique a B1) -- seuls les COURANTS
%                 en viennent ; le maillage garde son propre entrefer
%    reference  : sondes lues des cartes de champ ANSYS (.docx), valeurs
%                 NON re-derivables ici, declarees comme telles
clear; clc; t0v=tic;
diary('B7_probes_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G; W=ctx.W; BH=ctx.BH;
nT=17; nO=4; Nh=8192;
A1=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
ctx.AG=A1; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;

fprintf('=== B7 : sondes locales de la Table 17 ===\n');
fprintf('  machine   : MAS 48/44, 18,5 kW\n');
fprintf('  circuit   : base P1, N_h = %d, pavage nT=%d nO=%d, X_m0 = %.3f ohm\n', ...
    Nh,nT,nO,ctx.Xm0);
fprintf('  maillage  : mec.mesh_refined, nc = 6, nr = 3\n');
fprintf('  entrefer du maillage : mec.airgap_fourier, base P0 IMPOSEE\n\n');

fprintf('  ---- 0. la chaine de sondes n''est pas la chaine de performance ----\n');
fprintf('    reseau de performance : mec.airgap_dtn_tooth, base P1, N_h = %d\n',Nh);
fprintf('    chaine de sondes      : mec.airgap_fourier,   base P0, N_h = 100\n');
fprintf('    airgap_fourier.m:63   : basis = ''p0'' par defaut\n');
fprintf('    mesh_refined.m:271    : appel SANS argument de base\n');
fprintf('    mesh_refined.m:258    : NhF = 100 (M.opt.gap_fourier_Nh)\n');
fprintf('    => rapport de troncature entre les deux chaines : %.0f\n\n',Nh/100);

%% ---- points de fonctionnement (identiques a RUN_ARTICLE §IV) ---------
n1=1469.772776; s1=(M.ns*60-n1)/(M.ns*60);
r0 =mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
r1f=mec.equivalent_circuit(ctx,s1 ,ctx.Xm0);
fprintf('  points : a vide (s = 1e-4) | en charge (n = %.3f tr/min -> s = %.6f)\n', ...
    n1,s1);
fprintf('  courants : I1 a vide %.3f A | I1 en charge %.3f A\n\n',r0.I1,r1f.I1);

%  sondes de reference : LITTERAUX lus des cartes ANSYS (.docx). Elles ne
%  sont pas re-derivables depuis les .tab -- ce sont des maxima lus sur des
%  cartes. Declarees ici comme telles, avec leur provenance.
FEA=struct('lbl',{{'culasse stator','dent stator','dent rotor','culasse rotor'}}, ...
           'vide',[1.837 1.633 1.929 1.674], ...
           'chg' ,[1.896 1.684 1.870 1.577]);
fprintf('  reference : sondes lues des cartes de champ ANSYS (.docx),\n');
fprintf('              valeurs NON re-derivables depuis les .tab.\n\n');

%% ---- 1. balayage de la troncature ------------------------------------
fprintf('  ---- 1. les sondes derivent-elles avec la troncature ? ----\n');
NHL=[100 200 400 800];
P0m=nan(numel(NHL),4); P1m=nan(numel(NHL),4);
for k=1:numel(NHL)
    Mk=M; Mk.opt.gap_fourier_Nh=NHL(k);
    [i30,ib0]=b7_currents(Mk,W,r0);
    [i31,ib1]=b7_currents(Mk,W,r1f);
    me0=mec.mesh_refined(Mk,G,6,i30,0,3,ib0); Se0=mec.solve_mesh(me0,BH,Mk.opt);
    me1=mec.mesh_refined(Mk,G,6,i31,0,3,ib1); Se1=mec.solve_mesh(me1,BH,Mk.opt);
    [a1,a2,a3,a4]=b7_regional(me0,Se0); P0m(k,:)=[a2 a1 a3 a4];
    [b1,b2,b3,b4]=b7_regional(me1,Se1); P1m(k,:)=[b2 b1 b3 b4];
    fprintf('    N_h = %4d fait (%.0f s)\n',NHL(k),toc(t0v));
end

fprintf('\n    -- A VIDE : maxima regionaux MEC [T] --\n');
fprintf('    %-16s %10s %10s %10s %10s %12s %10s\n', ...
    'region','Nh=100','Nh=200','Nh=400','Nh=800','derive','FEA');
for j=1:4
    fprintf('    %-16s %10.4f %10.4f %10.4f %10.4f %11.2f %% %10.3f\n', ...
        FEA.lbl{j},P0m(1,j),P0m(2,j),P0m(3,j),P0m(4,j), ...
        100*(P0m(end,j)/P0m(1,j)-1),FEA.vide(j));
end
fprintf('\n    -- EN CHARGE : maxima regionaux MEC [T] --\n');
fprintf('    %-16s %10s %10s %10s %10s %12s %10s\n', ...
    'region','Nh=100','Nh=200','Nh=400','Nh=800','derive','FEA');
for j=1:4
    fprintf('    %-16s %10.4f %10.4f %10.4f %10.4f %11.2f %% %10.3f\n', ...
        FEA.lbl{j},P1m(1,j),P1m(2,j),P1m(3,j),P1m(4,j), ...
        100*(P1m(end,j)/P1m(1,j)-1),FEA.chg(j));
end

%% ---- 2. la Table 17 telle qu'elle serait publiee ---------------------
fprintf('\n  ---- 2. Table 17 au reglage de production (N_h = 100) ----\n');
fprintf('    %-16s %10s %10s %10s | %10s %10s %10s\n', ...
    'region','MEC vide','FEA vide','ecart','MEC chg','FEA chg','ecart');
for j=1:4
    fprintf('    %-16s %10.4f %10.3f %9.1f %% | %10.4f %10.3f %9.1f %%\n', ...
        FEA.lbl{j},P0m(1,j),FEA.vide(j),100*(P0m(1,j)/FEA.vide(j)-1), ...
        P1m(1,j),FEA.chg(j),100*(P1m(1,j)/FEA.chg(j)-1));
end

%% ---- 3. verdict ------------------------------------------------------
dmax=max(abs([100*(P0m(end,:)./P0m(1,:)-1), 100*(P1m(end,:)./P1m(1,:)-1)]));
fprintf('\n  ---- 3. verdict ----\n');
fprintf('    derive maximale des huit sondes sur un facteur 8 de troncature :\n');
fprintf('    %.2f %%\n\n',dmax);
if dmax>2
    fprintf(['    LES SONDES DERIVENT. La Table 17 est donc publiee a une\n' ...
             '    troncature (N_h = 100) qui n''est pas une valeur convergee,\n' ...
             '    dans une base P0 que la section 3 de l''Article I montre non\n' ...
             '    conforme, et sur des grandeurs LOCALES -- le cas le plus\n' ...
             '    defavorable des trois. Elle ne peut pas etre publiee sans\n' ...
             '    que sa configuration ne soit declaree ET sa derive donnee.\n']);
else
    fprintf(['    LES SONDES NE DERIVENT PAS a cette echelle de troncature.\n' ...
             '    Cela ne valide pas la base P0 : cela signifie que, comme sur\n' ...
             '    le bore du PMSM (Article I §3.5), la chaine opere dans un\n' ...
             '    regime ou le defaut ne se manifeste pas. La configuration\n' ...
             '    doit neanmoins etre DECLAREE dans la legende de la Table 17,\n' ...
             '    puisqu''elle differe de celle de toutes les autres tables.\n']);
end
fprintf(['\n    DANS LES DEUX CAS, deux points sont a porter au manuscrit :\n' ...
         '      - la Table 17 ne sort PAS de la meme chaine que les autres\n' ...
         '        tables de l''Article II (operateur, base et troncature\n' ...
         '        differents) : "une grandeur, une chaine" impose de le dire ;\n' ...
         '      - les huit valeurs de reference sont des maxima LUS SUR DES\n' ...
         '        CARTES, non des sorties de .tab : elles ne sont pas\n' ...
         '        re-derivables et doivent etre declarees comme telles.\n' ...
         '    Le bloc C1 devra en outre promouvoir inst_currents / surfU /\n' ...
         '    regional_max dans +mec : ce sont aujourd''hui des fonctions\n' ...
         '    LOCALES de RUN_ARTICLE.m, donc non appelables ailleurs, ce qui\n' ...
         '    force toute reprise a les dupliquer -- comme ce bloc a du le\n' ...
         '    faire.\n']);

save('B7_probes.mat','NHL','P0m','P1m','FEA','s1');
fprintf('\n  duree %.0f s\n=== B7 termine ===\n',toc(t0v));
diary off;

% ======================================================================
%  Reprises A L'IDENTIQUE des fonctions locales de RUN_ARTICLE.m
%  (l. 508-513 et 531-539). Elles y sont LOCALES a un script, donc non
%  appelables : la duplication est subie, pas choisie. Bloc C1.
function [i3,ib]=b7_currents(M,W,r)
    psi1=angle(r.I1c); psi2=angle(r.I2c);
    i3=sqrt(2)*r.I1*[cos(psi1);cos(psi1-2*pi/3);cos(psi1+2*pi/3)];
    Ibar=r.I2*(2*M.m*W.kw1*W.Nph)/M.Nr;
    ib=sqrt(2)*Ibar*cos(M.p*(2*pi*(0:M.Nr-1)'/M.Nr)-psi2);
end
function [Bts,Bys,Btr,Byr]=b7_regional(me,Se)
    nb=me.gapfirst-1; isFe=logical(me.iron(1:nb));
    a=me.a(1:nb); B=abs(Se.B(1:nb));
    Ms=me.Ms; Ls=me.Ls; Mr=me.Mr; nst=me.nr; nrt=me.nr;
    laS=ceil(a/Ms); laS(a>Ms*Ls)=0;
    aR=a-Ms*Ls; laR=ceil(aR/Mr); laR(a<=Ms*Ls)=0;
    Bts=max(B(isFe&laS>=1&laS<=nst)); Bys=max(B(isFe&laS>nst));
    Btr=max(B(isFe&laR>=1&laR<=nrt)); Byr=max(B(isFe&laR>nrt));
end
