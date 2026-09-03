%% RUN_R9_ARCHIVE - R9 : assembler l'archive de depot public
%
%  CE QUE LA v4 DEMANDE. Une archive contenant le manifeste date, les scripts
%  de production des deux articles, les sorties diary de chaque bloc, un
%  README donnant l'ordre d'execution et les versions, et une licence.
%
%  CE QUE CE SCRIPT FAIT. Il COPIE ; il ne reecrit rien. Le code archive est
%  octet pour octet celui qui a produit les nombres publies -- c'est la
%  condition de la tracabilite, et c'est pourquoi l'adaptation des chemins
%  est confiee a un utilitaire separe (zenodo/code/SET_REFERENCE_PATH.m)
%  plutot qu'appliquee ici.
%
%  IL INCLUT LA REFERENCE. Les exports ANSYS (.tab) pesent 4,7 Mo au total :
%  les embarquer rend l'archive autonome, donc verifiable sans licence de
%  solveur. C'est le point qui separe "le lecteur peut lire le code" de
%  "le lecteur peut refaire le calcul".
%
%  GARDE EXECUTABLE. Chaque sortie .txt et chaque script nommes dans
%  MANIFEST.md doivent se retrouver dans l'archive. Un manifeste qui renvoie
%  a un fichier absent n'est pas un manifeste.
%
%  GARDE NON EXECUTABLE, ET DECLAREE COMME TELLE. La v4 exige qu'un TIERS
%  regenere une grandeur de chaque article avec le seul README. Je ne peux
%  pas passer ce test : j'ai ecrit le code. Il est enonce au §9 du README et
%  reste OUVERT.
clear; clc; t0=tic;
diary('R9_archive_out.txt'); diary on;

HERE = fileparts(mfilename('fullpath'));
if isempty(HERE), HERE = pwd; end
ROOT = fileparts(HERE);                       % ...\MEC
Z    = fullfile(HERE,'zenodo');
%  Trois emplacements, non deux : quelques blocs transverses aux deux
%  articles (temps de calcul, audit des ecarts) vivent dans article/ et non
%  dans un dossier machine. La garde du §6 l'a signale a la premiere
%  execution -- RUN_T9_TEMPS.m et sa sortie manquaient.
%
%  L'ARBORESCENCE EST UN MIROIR, et ce n'est pas cosmetique. Cinq gardes
%  relisent leur manuscrit par un chemin RELATIF -- fullfile('..','article',
%  '*.tex') depuis MEC_BLDC/ ou MEC_IM/. code/ joue donc le role de MEC/ et
%  les trois dossiers gardent leur nom d'origine, sans quoi ces gardes ne
%  peuvent pas s'executer chez un tiers.
SRC  = { fullfile(ROOT,'MEC_BLDC'), 'MEC_BLDC'; ...
         fullfile(ROOT,'MEC_IM'),   'MEC_IM'; ...
         HERE,                      'article' };
TEXF = {'ArticleI_DtN_PMSM.tex','ArticleII_Carter_IM.tex'};
REF  = { 'C:\Users\hp\Desktop\ANSYS résultat 750W',   'ANSYS_750W'; ...
         'C:\Users\hp\Desktop\ANSYS résultat 18.5KW', 'ANSYS_18_5kW' };
NOTES= {'MANIFEST.md','R1_NOTE.md','R2_NOTE.md','R3_NOTE.md','R4_NOTE.md', ...
        'R5_NOTE.md','R6_NOTE.md','R7_NOTE.md','R8_NOTE.md', ...
        'RAPPORT_SPEC_v2.md','RAPPORT_FINAL_HANDOFF.md'};

fprintf('=== R9 : assemblage de l''archive de depot ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    destination  : %s\n',Z);
fprintf('    MATLAB       : %s (%s)\n',version,version('-release'));
fprintf('    produits     : MATLAB de base seul (mesure, voir README §2)\n');
fprintf('    ANSYS        : Electronics Desktop 2023 R1 (lu dans les .aedt)\n');
fprintf('    code         : COPIE INCHANGE, chemins adaptes par SET_REFERENCE_PATH\n');

%% ---- 1. arborescence --------------------------------------------------
D = {'code','outputs','reference','notes'};
for k=1:numel(D)
    p=fullfile(Z,D{k}); if ~isfolder(p), mkdir(p); end
end
fprintf('\n  ---- 1. ARBORESCENCE ----\n');
for k=1:numel(D), fprintf('    zenodo/%s/\n',D{k}); end

%% ---- 2. code ----------------------------------------------------------
fprintf('\n  ---- 2. CODE (copie inchange) ----\n');
nm=0;
for k=1:size(SRC,1)
    dst=fullfile(Z,'code',SRC{k,2}); if ~isfolder(dst), mkdir(dst); end
    f=dir(fullfile(SRC{k,1},'*.m'));
    for i=1:numel(f), copyfile(fullfile(f(i).folder,f(i).name),dst); end
    n=numel(f);
    pk=fullfile(SRC{k,1},'+mec');                     % paquet du MAS
    if isfolder(pk)
        dpk=fullfile(dst,'+mec'); if ~isfolder(dpk), mkdir(dpk); end
        g=dir(fullfile(pk,'*.m'));
        for i=1:numel(g), copyfile(fullfile(g(i).folder,g(i).name),dpk); end
        n=n+numel(g);
        fprintf('    %-10s %3d fichiers .m (dont %d dans +mec/)\n',SRC{k,2},n,numel(g));
    else
        fprintf('    %-10s %3d fichiers .m\n',SRC{k,2},n);
    end
    nm=nm+n;
end
fprintf('    total %d scripts\n',nm);

%% ---- 3. sorties -------------------------------------------------------
fprintf('\n  ---- 3. SORTIES D''EXECUTION ----\n');
nt=0;
for k=1:size(SRC,1)
    dst=fullfile(Z,'outputs',SRC{k,2}); if ~isfolder(dst), mkdir(dst); end
    f=dir(fullfile(SRC{k,1},'*.txt'));
    %  On n'archive pas le journal de la presente execution : il est ouvert
    %  en ecriture, la copie serait tronquee. Il est ecrit a cote.
    f=f(~strcmpi({f.name},'R9_archive_out.txt'));
    for i=1:numel(f), copyfile(fullfile(f(i).folder,f(i).name),dst); end
    fprintf('    %-10s %3d fichiers .txt\n',SRC{k,2},numel(f));
    nt=nt+numel(f);
end
fprintf('    total %d sorties\n',nt);

%% ---- 4. reference EF --------------------------------------------------
fprintf('\n  ---- 4. REFERENCE ELEMENTS FINIS (.tab seuls) ----\n');
nr=0; szr=0;
for k=1:size(REF,1)
    if ~isfolder(REF{k,1})
        fprintf('    [!] %s INTROUVABLE -- archive incomplete\n',REF{k,1});
        continue;
    end
    f=dir(fullfile(REF{k,1},'**','*.tab'));
    for i=1:numel(f)
        rel=erase(f(i).folder,[REF{k,1} filesep]);
        if strcmp(rel,f(i).folder), rel=''; end       % fichier a la racine
        dst=fullfile(Z,'reference',REF{k,2},rel);
        if ~isfolder(dst), mkdir(dst); end
        copyfile(fullfile(f(i).folder,f(i).name),dst);
        szr=szr+f(i).bytes;
    end
    fprintf('    %-14s %3d fichiers .tab\n',REF{k,2},numel(f));
    nr=nr+numel(f);
end
fprintf('    total %d fichiers, %.2f Mo\n',nr,szr/1e6);

%% ---- 5. notes ---------------------------------------------------------
fprintf('\n  ---- 5. NOTES ET MANIFESTE ----\n');
nn=0;
for k=1:numel(NOTES)
    s=fullfile(HERE,NOTES{k});
    if ~isfile(s), s=fullfile(ROOT,NOTES{k}); end     % MANIFEST.md est a la racine
    if isfile(s)
        copyfile(s,fullfile(Z,'notes')); nn=nn+1;
        fprintf('    %s\n',NOTES{k});
    else
        fprintf('    [!] %s ABSENT\n',NOTES{k});
    end
end
fprintf('    total %d documents\n',nn);

%% ---- 5b. LES DEUX MANUSCRITS, et la decision qu'ils portent ----------
fprintf('\n  ---- 5b. MANUSCRITS ----\n');
nx=0;
for k=1:numel(TEXF)
    s=fullfile(HERE,TEXF{k});
    if isfile(s)
        copyfile(s,fullfile(Z,'code','article')); nx=nx+1;
        d=dir(s); fprintf('    %-28s %6.0f ko\n',TEXF{k},d.bytes/1024);
    else
        fprintf('    [!] %s ABSENT\n',TEXF{k});
    end
end
fprintf('    Ils sont inclus parce que CINQ gardes les relisent pour se\n');
fprintf('    verifier (R5, R6a, R8, A5, B6). Les retirer rend l''archive\n');
fprintf('    incapable de se controler elle-meme.\n');
fprintf('    DECISION D''AUTEUR AVANT DEPOT : verifier que l''accord\n');
fprintf('    editeur autorise leur diffusion. Sinon, les retirer et le\n');
fprintf('    declarer -- README §4 indique alors quels tests tombent.\n');

%% ---- 5c. DEPENDANCES NON SATISFAITES, MESUREES -----------------------
%  Plutot que d'affirmer que l'archive est autonome, on cherche ce qui n'y
%  est pas. Trois familles de references pointent hors archive.
fprintf('\n  ---- 5c. DEPENDANCES NON SATISFAITES ----\n');
PAT = {'ANSYS-\',            'projet .aedt / .aedtresults (non inclus)'; ...
       'conception\New Folder','dossier externe de code analytique (non inclus)'; ...
       'MEC_DtN_paper_v2.tex','manuscrit herite, archive (non inclus)'};
mm=dir(fullfile(Z,'code','**','*.m'));
%  Le scanner porte lui-meme les trois motifs dans sa table : il se
%  detecterait. On l'exclut, et on le dit plutot que de laisser un faux
%  positif dans le livrable.
mm=mm(~strcmpi({mm.name},'RUN_R9_ARCHIVE.m'));
nd=0;
for p=1:size(PAT,1)
    hits={};
    for i=1:numel(mm)
        f=fullfile(mm(i).folder,mm(i).name);
        tx2=fileread(f);
        if contains(tx2,PAT{p,1})
            %  on ne retient que les lignes de CODE, pas les commentaires
            L=regexp(tx2,'\r\n|\n|\r','split');
            j=find(contains(L,PAT{p,1}) & ~startsWith(strtrim(L),'%'));
            if ~isempty(j)
                hits{end+1}=sprintf('%s (L%s)',mm(i).name, ...
                    strjoin(string(j),',')); %#ok<SAGROW>
            end
        end
    end
    fprintf('    %-42s : %d script(s)\n',PAT{p,2},numel(hits));
    for i=1:numel(hits), fprintf('        %s\n',hits{i}); end
    nd=nd+numel(hits);
end
fprintf('    total %d scripts ne peuvent pas s''executer dans l''archive.\n',nd);
fprintf('    Ils sont conserves -- ils documentent des diagnostics reels --\n');
fprintf('    mais le README doit les nommer plutot que les laisser decouvrir.\n');

%% ---- 6. GARDE : le manifeste renvoie-t-il a des fichiers presents ? ----
fprintf('\n  ---- 6. GARDE : completude du manifeste ----\n');
mf=fullfile(Z,'notes','MANIFEST.md');
if ~isfile(mf)
    fprintf('    GARDE IMPOSSIBLE : MANIFEST.md absent de l''archive.\n');
    G=false;
else
    tx=fileread(mf);
    outs=unique(regexp(tx,'[A-Za-z0-9_]+_out\.txt','match'));
    scr =unique(regexp(tx,'RUN_[A-Za-z0-9_]+\.m','match'));
    have=[dir(fullfile(Z,'outputs','**','*.txt'))];
    haveN={have.name};
    hs  =[dir(fullfile(Z,'code','**','*.m'))];
    hsN ={hs.name};
    mo=outs(~ismember(outs,haveN));
    ms=scr(~ismember(scr,hsN));
    fprintf('    sorties citees par le manifeste : %d | absentes : %d\n',numel(outs),numel(mo));
    for i=1:numel(mo), fprintf('      MANQUE %s\n',mo{i}); end
    fprintf('    scripts cites par le manifeste  : %d | absents  : %d\n',numel(scr),numel(ms));
    for i=1:numel(ms), fprintf('      MANQUE %s\n',ms{i}); end
    G=isempty(mo)&&isempty(ms);
    if G, fprintf('    GARDE PASSEE : tout ce que le manifeste nomme est present.\n');
    else, fprintf('    GARDE ECHOUEE : le manifeste renvoie a des fichiers absents.\n'); end
end

%% ---- 7. inventaire ----------------------------------------------------
fprintf('\n  ---- 7. INVENTAIRE ----\n');
a=dir(fullfile(Z,'**','*')); a=a(~[a.isdir]);
fprintf('    %d fichiers, %.2f Mo au total\n',numel(a),sum([a.bytes])/1e6);
for k=1:numel(D)
    b=dir(fullfile(Z,D{k},'**','*')); b=b(~[b.isdir]);
    fprintf('    %-11s %4d fichiers %8.2f Mo\n',[D{k} '/'],numel(b),sum([b.bytes])/1e6);
end
rt=dir(fullfile(Z,'*')); rt=rt(~[rt.isdir]);
for i=1:numel(rt), fprintf('    %-11s %s\n','(racine)',rt(i).name); end

%% ---- 8. CE QUI RESTE OUVERT ------------------------------------------
fprintf('\n  ---- 8. CE QUI RESTE OUVERT ----\n');
fprintf('    La garde d''acceptation de la v4 -- un TIERS regenere une\n');
fprintf('    grandeur de chaque article avec le seul README -- n''est PAS\n');
fprintf('    passee. Elle ne peut pas l''etre par qui a ecrit le code.\n');
fprintf('    Le test a lui soumettre est au README §4 : RUN_R5_NSH (2 s)\n');
fprintf('    et RUN_R8_TABLE2 (95 s), qui impriment chacun leur garde.\n');
fprintf('    Tant qu''il n''a pas ete passe, l''archive est ASSEMBLEE mais\n');
fprintf('    NON VALIDEE comme reproductible.\n');
fprintf('\n    La licence est un MODELE : la ligne de copyright est vide et\n');
fprintf('    le choix appartient a l''auteur. Aucun depot n''a ete effectue.\n');

fprintf('\n  duree %.0f s\n=== R9 termine ===\n',toc(t0));
diary off;
