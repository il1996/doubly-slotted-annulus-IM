function SET_REFERENCE_PATH(mode)
%SET_REFERENCE_PATH  Pointe le code archive vers les references de l'archive.
%
%   POURQUOI CE FICHIER EXISTE. Les scripts sont archives INCHANGES, octet
%   pour octet identiques a ceux qui ont produit les nombres publies. Ils
%   portent donc les chemins absolus de la machine d'origine. Les modifier
%   avant archivage aurait casse la tracabilite ; on les adapte ici, une
%   fois, de facon visible et reversible.
%
%   SET_REFERENCE_PATH            applique et liste chaque fichier modifie
%   SET_REFERENCE_PATH('dryrun')  montre ce qui serait fait, sans rien ecrire
%
%   Ce que la fonction remplace, et rien d'autre :
%     'C:\Users\hp\Desktop\ANSYS résultat 750W'    -> <archive>/reference/ANSYS_750W
%     'C:\Users\hp\Desktop\ANSYS résultat 18.5KW'  -> <archive>/reference/ANSYS_18_5kW
%
%   Les sous-dossiers gardent leur nom d'origine, accents et parentheses
%   compris : les scripts les nomment explicitement.
%
%   Un fichier .bak est ecrit a cote de chaque fichier modifie, pour que
%   l'operation soit annulable.

if nargin<1, mode='apply'; end
dry = strcmpi(mode,'dryrun');

here = fileparts(mfilename('fullpath'));          % <archive>/code
root = fileparts(here);                           % <archive>
ref  = fullfile(root,'reference');

sub = { 'C:\Users\hp\Desktop\ANSYS résultat 750W',   fullfile(ref,'ANSYS_750W'); ...
        'C:\Users\hp\Desktop\ANSYS résultat 18.5KW', fullfile(ref,'ANSYS_18_5kW') };

for k=1:size(sub,1)
    if ~isfolder(sub{k,2})
        warning('SET_REFERENCE_PATH:absent', ...
            'Dossier de reference introuvable : %s',sub{k,2});
    end
end

fprintf('=== SET_REFERENCE_PATH (%s) ===\n',upper(mode));
fprintf('  archive : %s\n',root);
for k=1:size(sub,1)
    fprintf('  %-46s -> %s\n',sub{k,1},sub{k,2});
end
fprintf('\n');

files = dir(fullfile(here,'**','*.m'));
me = [mfilename('fullpath') '.m'];
nf = 0; nl = 0;
for i=1:numel(files)
    f = fullfile(files(i).folder,files(i).name);
    if strcmp(f,me), continue; end
    txt = fileread(f);
    hit = false;
    for k=1:size(sub,1)
        if contains(txt,sub{k,1}), hit = true; end
    end
    if ~hit, continue; end

    lines = regexp(txt,'\r\n|\n|\r','split');
    changed = false(1,numel(lines));
    for j=1:numel(lines)
        for k=1:size(sub,1)
            if contains(lines{j},sub{k,1})
                lines{j} = strrep(lines{j},sub{k,1},sub{k,2});
                changed(j) = true;
            end
        end
    end
    nf = nf + 1; nl = nl + sum(changed);
    rel = erase(f,[root filesep]);
    fprintf('  %s  (%d ligne(s))\n',rel,sum(changed));
    for j=find(changed)
        fprintf('      L%-4d %s\n',j,strtrim(lines{j}));
    end
    if ~dry
        copyfile(f,[f '.bak']);
        fid = fopen(f,'w','n','UTF-8');
        fprintf(fid,'%s\n',lines{1:end-1});
        fprintf(fid,'%s',lines{end});
        fclose(fid);
    end
end

fprintf('\n  %d fichier(s), %d ligne(s)%s\n',nf,nl, ...
    string_if(dry,'  [DRYRUN : rien ecrit]',''));
if ~dry && nf>0
    fprintf('  Un .bak a ete ecrit a cote de chaque fichier modifie.\n');
end
fprintf('=== termine ===\n');
end

function s = string_if(c,a,b)
if c, s=a; else, s=b; end
end
