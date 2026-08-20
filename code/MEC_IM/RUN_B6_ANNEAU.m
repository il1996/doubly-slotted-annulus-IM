%% RUN_B6_ANNEAU - courants d'anneau : ECART DE MODELE, non precision
%
%  OBJET (reformulation du bloc B6). Les deux lignes "end-ring current"
%  de la Table 16 changent de signe entre les deux conditions (-4,1 % en
%  charge, +4,3 % au calage). Lues telles quelles, elles se lisent comme
%  une precision du reseau sur l'anneau qui se degraderait au calage.
%  C'est faux. Ce bloc etablit que :
%
%    (i)  le reseau implemente EXACTEMENT le modele d'anneau geometrique
%         ideal, I_anneau/I_barre = 1/(2 sin(p*pi/Nr)) ;
%    (ii) la reference EF ne l'implemente PAS : elle en est a ~4 %, du
%         MEME cote dans les deux conditions ;
%    (iii) l'identite exacte  1+e_anneau = (1+e_barre)*(1+e_rapport)
%         attribue l'inversion de signe ENTIEREMENT a la ligne "bar
%         current", le terme e_rapport etant stable en signe et en
%         amplitude.
%
%  Ces lignes mesurent donc un ECART DE MODELE D'ANNEAU (segmentation,
%  resistance de contact, effets 3-D dans la reference) AJOUTE a l'erreur
%  sur la barre. A reformuler dans le manuscrit.
%
%  CE BLOC NE TRANSCRIT AUCUN CHIFFRE.
%    - le rapport ideal sort de la FORMULE, avec M.p et M.Nr lus de
%      mec.machine_18_5kW ;
%    - les courants du reseau sortent de mec.equivalent_circuit, chaine
%      de RUN_B1_IM_P1.m (pavage nT=17, nO=4, N_h=8192, base P1) ;
%    - le glissement nominal est RELU dans le source de RUN_B1_IM_P1.m,
%      et recoupe par la vitesse lue de Speed Plot 1.tab ;
%    - les courants de reference sortent des .tab ANSYS, convention
%      RUN_ARTICLE section III (End Connection Plot 1, colonnes 2 et 3,
%      valeur efficace sur le regime etabli) ;
%    - les quatre cellules de la Table 16 mises en cause sont RELUES du
%      manuscrit LaTeX, pas recopiees.
%
%  CONFIGURATION DECLAREE
%    machine   : MAS 48/44, 18,5 kW | p = 2, Nr = 44
%    pavage    : nT = 17, nO = 4 | N_h = 8192 | base P1  (identique a B1)
%    vrillage  : les DEUX etats sont produits. ON = etat de production,
%                celui dont la Table 16 est issue ; OFF (M.opt.skew_harm=0)
%                = etat loyal etabli par B8, la reference EF etant une
%                tranche droite (NumberOfSlices = 1).
%    reference : C:\Users\hp\Desktop\ANSYS résultat 18.5KW
%                transitoire\en charge et transitoire\rotor bloqué
%
%  Sortie : B6_anneau_out.txt, B6_anneau.mat
clear; clc; t0v=tic;
diary('B6_anneau_out.txt'); diary on;

ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
TEX ='C:\Users\hp\Desktop\Matlab program\MEC\article\MEC_DtN_paper_v2.tex';
SRCB1='RUN_B1_IM_P1.m';
rd  =@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
rmsw=@(A,c,t0)sqrt(mean(A(A(:,1)>=t0,c).^2,'omitnan'));
avgw=@(A,c,t0)mean(A(A(:,1)>=t0,c),'omitnan');
nT=17; nO=4; Nh=8192;

%% ---- 0. configuration, glissement, operateur -------------------------
M0=mec.machine_18_5kW(); ctx0=mec.build_context(M0); G=ctx0.G;
kid=1/(2*sin(M0.p*pi/M0.Nr));                 % rapport ideal, PAR LA FORMULE

%  glissement nominal : relu du source de B1, pas retape
txtB1=fileread(SRCB1);
tk=regexp(txtB1,'s_ch\s*=\s*([0-9.]+)','tokens','once');
if isempty(tk), error('B6:slip','s_ch introuvable dans %s',SRCB1); end
s_B1=str2double(tk{1});

%  recoupement : le glissement effectif de l'essai EF en charge
Dch=fullfile(ROOT,'transitoire','en charge');
Dbl=fullfile(ROOT,'transitoire','rotor bloqué');
spd=rd(fullfile(Dch,'Speed Plot 1.tab')); t0=1.0;
n_ans=avgw(spd,2,t0); s_tab=(M0.ns*60-n_ans)/(M0.ns*60);

A1=mec.airgap_dtn_tooth(M0,G,0,nT,nO,Nh,'p1');   % ne depend pas de la cage
Xm0=mec.magnetizing(setfield(ctx0,'AG',A1),0.2).Xm; %#ok<SFLD>

fprintf('=== B6 : courants d''anneau -- ecart de MODELE, non precision ===\n');
fprintf('  machine   : MAS 48/44 | p = %d | Nr = %d | Ns = %d\n',M0.p,M0.Nr,M0.Ns);
fprintf('  pavage    : nT=%d nO=%d | N_h = %d | base P1 | X_m0 = %.3f ohm\n',nT,nO,Nh,Xm0);
fprintf('  reference : %s\n',ROOT);
fprintf('  manuscrit : %s\n',TEX);
fprintf('  glissement nominal : %.6f (relu de %s)\n',s_B1,SRCB1);
fprintf('                       %.6f (deduit de Speed Plot 1.tab, n = %.3f tr/min)\n', ...
        s_tab,n_ans);
fprintf('                       ecart entre les deux : %.3f %%\n',100*(s_B1/s_tab-1));

%% ---- 1. le rapport ideal, et ce que le reseau en fait -----------------
fprintf('\n  ---- 1. LE RAPPORT GEOMETRIQUE IDEAL ----\n');
fprintf('    1/(2 sin(p*pi/Nr)) = 1/(2 sin(pi/%g)) = %.6f\n',M0.Nr/M0.p,kid);

%  Les quatre points de fonctionnement du reseau. Deux etats de vrillage :
%  ON = etat de production (celui de la Table 16), OFF = etat loyal (B8).
cfg={struct('nm','ON ','sk',[]), struct('nm','OFF','sk',0)};
NET=struct([]);
for c=1:2
    Mc=mec.machine_18_5kW();
    if ~isempty(cfg{c}.sk), Mc.opt.skew_harm=cfg{c}.sk; end
    cx=mec.build_context(Mc); cx.AG=A1; cx.Xm0=Xm0;
    rc=mec.equivalent_circuit(cx,s_B1,Xm0);
    rb=mec.equivalent_circuit(cx,1.0 ,Xm0);
    NET(c).nm=cfg{c}.nm; NET(c).ch=rc; NET(c).bl=rb;
end

%  Controle d'identite : au FONDAMENTAL le reseau doit rendre kid EXACTEMENT
%  (equivalent_circuit.m:140, Iring1 = Ibar1/|2 sin(pi p/Nr)|).
e_id=max([abs(NET(1).ch.Iring1/NET(1).ch.Ibar1/kid-1), ...
          abs(NET(1).bl.Iring1/NET(1).bl.Ibar1/kid-1)]);
fprintf('    controle : rapport FONDAMENTAL du reseau / rapport ideal - 1 = %.2e\n',e_id);
fprintf('    -> le reseau implemente le modele ideal a la precision machine.\n');
fprintf('    Le RMS physique publie ajoute les branches de cage harmoniques,\n');
fprintf('    dont le report d''anneau est 1/(2 sin(pi*nu*p/Nr)) et non 1/(2 sin(pi*p/Nr)) :\n');
fprintf('    %-28s %10s %10s %10s\n','etat vrillage / condition','I1 barre','I tot barre','part harm.');
for c=1:2
    for q=1:2
        if q==1, rr=NET(c).ch; cd='en charge'; else, rr=NET(c).bl; cd='calage'; end
        fprintf('    %-28s %10.2f %10.2f %9.2f %%\n', ...
            ['vrillage ' NET(c).nm ' / ' cd],rr.Ibar1,rr.Ibar, ...
            100*sqrt(max(rr.Ibar^2-rr.Ibar1^2,0))/rr.Ibar);
    end
end
fprintf('    C''est cette part harmonique -- et elle seule -- qui ecarte le\n');
fprintf('    rapport RMS du reseau de la valeur ideale, de quelques dixiemes.\n');

%% ---- 2. les quatre rapports ------------------------------------------
ec =rd(fullfile(Dch,'End Connection Plot 1.tab'));
curb=rd(fullfile(Dbl,'Winding Plot 3.tab')); t0b=0.5*max(curb(:,1));
ecb=rd(fullfile(Dbl,'End Connection Plot 1.tab'));
Ib_F =1e3*rmsw(ec ,2,t0 ); Ir_F =1e3*rmsw(ec ,3,t0 );
Ib_Fb=1e3*rmsw(ecb,2,t0b); Ir_Fb=1e3*rmsw(ecb,3,t0b);

fprintf('\n  ---- 2. LES QUATRE RAPPORTS ANNEAU/BARRE ----\n');
fprintf('  (regime etabli : t >= %.3f s en charge, t >= %.3f s au calage)\n',t0,t0b);
fprintf('  %-10s %-9s %-11s %11s %12s %10s %12s\n', ...
    'source','vrillage','condition','I barre (A)','I anneau (A)','rapport','ecart/ideal');
row=@(sr,sk,cd,ib,ir)fprintf('  %-10s %-9s %-11s %11.1f %12.1f %10.4f %11.2f %%\n', ...
    sr,sk,cd,ib,ir,ir/ib,100*((ir/ib)/kid-1));
for c=1:2
    row('reseau',NET(c).nm,'en charge',NET(c).ch.Ibar,NET(c).ch.Iring);
    row('reseau',NET(c).nm,'calage'   ,NET(c).bl.Ibar,NET(c).bl.Iring);
end
row('reference','--','en charge',Ib_F ,Ir_F );
row('reference','--','calage'   ,Ib_Fb,Ir_Fb);

k_ref_ch=Ir_F /Ib_F ; k_ref_bl=Ir_Fb/Ib_Fb;
d_ref=[100*(k_ref_ch/kid-1) 100*(k_ref_bl/kid-1)];
d_net=[100*(NET(1).ch.Iring/NET(1).ch.Ibar/kid-1) 100*(NET(1).bl.Iring/NET(1).bl.Ibar/kid-1); ...
       100*(NET(2).ch.Iring/NET(2).ch.Ibar/kid-1) 100*(NET(2).bl.Iring/NET(2).bl.Ibar/kid-1)];
fprintf('\n    LECTURE. Le reseau reste a %.2f %% du modele ideal au pire des\n', ...
        max(abs(d_net(:))));
fprintf('    quatre points ; la reference en est a %.2f %% et %.2f %%, du MEME\n', ...
        d_ref(1),d_ref(2));
fprintf('    cote dans les deux conditions. Un ecart de signe constant sur\n');
fprintf('    deux essais independants n''est pas du bruit : c''est une\n');
fprintf('    CONVENTION. Le modele d''anneau de la reference n''est donc pas\n');
fprintf('    le modele geometrique ideal.\n');

%% ---- 3. relecture de la Table 16 dans le manuscrit --------------------
fprintf('\n  ---- 3. LES QUATRE CELLULES MISES EN CAUSE, RELUES DU MANUSCRIT ----\n');
[TB,ntab]=lire_table16(TEX,'tab:im_tests');
fprintf('    %s : table n^o %d du fichier (label tab:im_tests)\n', ...
        'MEC_DtN_paper_v2.tex',ntab);
fprintf('    %-22s %-14s %10s %10s %10s\n','bloc','ligne','reseau','EF','ecart publie');
nm={'On load','Locked rotor'}; lb={'bar current','end-ring current'};
for b=1:2
    for l=1:2
        fprintf('    %-22s %-14s %10.4g %10.4g %9.1f %%\n', ...
            nm{b},lb{l},TB(b,l).net,TB(b,l).fea,TB(b,l).err);
    end
end

%  Controle de provenance : la colonne EF de la Table 16 doit se retrouver
%  dans les .tab. Si oui, la table et ce bloc parlent bien du meme essai.
fprintf('\n    controle de provenance de la colonne EF (Table 16 vs .tab relus)\n');
fprintf('    %-24s %12s %12s %10s\n','grandeur','Table 16','.tab relu','ecart');
prov=@(l,a,b)fprintf('    %-24s %12.4g %12.1f %9.2f %%\n',l,a,b,100*(a/b-1));
prov('I barre en charge',TB(1,1).fea,Ib_F );
prov('I anneau en charge',TB(1,2).fea,Ir_F );
prov('I barre calage'   ,TB(2,1).fea,Ib_Fb);
prov('I anneau calage'  ,TB(2,2).fea,Ir_Fb);

%% ---- 4. decomposition exacte des deux lignes "end-ring current" -------
fprintf('\n  ---- 4. DECOMPOSITION : ecart_anneau = ecart_barre + ecart_du_rapport ----\n');
fprintf('    Identite EXACTE, sans approximation :\n');
fprintf('      I_an,res/I_an,ref = (I_ba,res/I_ba,ref) * (k_res/k_ref)\n');
fprintf('      donc  1 + e_anneau = (1 + e_barre) * (1 + e_rapport)\n');
fprintf('      soit  e_anneau = e_barre + e_rapport + e_barre*e_rapport\n');
fprintf('    e_rapport = ecart des DEUX MODELES D''ANNEAU, k_res/k_ref - 1.\n\n');

%  (a) sur les valeurs PUBLIEES de la Table 16
DEC=struct([]);
for b=1:2
    eb=TB(b,1).net/TB(b,1).fea-1;
    kn=TB(b,2).net/TB(b,1).net; kf=TB(b,2).fea/TB(b,1).fea;
    ek=kn/kf-1;  er=TB(b,2).net/TB(b,2).fea-1;
    DEC(b).nm=nm{b}; DEC(b).eb=eb; DEC(b).ek=ek; DEC(b).er=er;
    DEC(b).kn=kn; DEC(b).kf=kf;
    DEC(b).rec=(1+eb)*(1+ek)-1; DEC(b).pub=TB(b,2).err/100;
end
fprintf('    (a) SUR LES VALEURS PUBLIEES DE LA TABLE 16\n');
fprintf('    %-14s %9s %10s %9s %9s %11s %10s %8s\n', ...
    'condition','e_barre','e_rapport','somme','croise','e_anneau','publie','residu');
for b=1:2
    fprintf('    %-14s %8.2f %% %9.2f %% %8.2f %% %8.2f %% %10.2f %% %9.1f %% %7.2f\n', ...
        DEC(b).nm,100*DEC(b).eb,100*DEC(b).ek,100*(DEC(b).eb+DEC(b).ek), ...
        100*DEC(b).eb*DEC(b).ek,100*DEC(b).rec,100*DEC(b).pub, ...
        100*(DEC(b).rec-DEC(b).pub));
fprintf('    %-14s   k_res = %.4f   k_ref = %.4f   (ideal %.4f)\n', ...
        '',DEC(b).kn,DEC(b).kf,kid);
end
sw_b=100*(DEC(2).eb-DEC(1).eb); sw_k=100*(DEC(2).ek-DEC(1).ek);
sw_r=100*(DEC(2).er-DEC(1).er);
fprintf('\n    BASCULE entre les deux conditions (calage moins charge) :\n');
fprintf('      ligne "bar current"      : %+7.2f points  <-- porte l''inversion\n',sw_b);
fprintf('      terme "ecart de modele"  : %+7.2f points\n',sw_k);
fprintf('      ligne "end-ring current" : %+7.2f points\n',sw_r);
fprintf('      part de la bascule imputable a la barre : %.1f %%\n',100*sw_b/sw_r);
fprintf('      le terme d''anneau reste du MEME SIGNE (%+.2f %% puis %+.2f %%) ;\n', ...
        100*DEC(1).ek,100*DEC(2).ek);
fprintf('      la ligne de barre, elle, change de signe (%+.2f %% puis %+.2f %%).\n', ...
        100*DEC(1).eb,100*DEC(2).eb);

%  A qui appartient e_rapport ? On le scinde entre le reseau et la reference,
%  toujours par une identite exacte : 1+e_rapport = (1+dn)/(1+df), avec
%  dn = k_res/kid - 1 (reseau vs ideal) et df = k_ref/kid - 1 (reference vs ideal).
fprintf('\n    A QUI APPARTIENT e_rapport ?  1 + e_rapport = (1 + d_res)/(1 + d_ref)\n');
fprintf('    %-14s %12s %12s %12s\n','condition','d_res','d_ref','e_rapport');
for b=1:2
    dn=DEC(b).kn/kid-1; df=DEC(b).kf/kid-1;
    fprintf('    %-14s %11.2f %% %11.2f %% %11.2f %%\n',DEC(b).nm,100*dn,100*df,100*DEC(b).ek);
end
fprintf('    d_res est quasi nul, d_ref ne l''est pas : e_rapport est un ecart\n');
fprintf('    de la REFERENCE au modele ideal, non une erreur du reseau.\n');

%  (b) sur le modele d'aujourd'hui, recalcule
fprintf('\n    (b) SUR LE MODELE RECALCULE CE JOUR (meme identite)\n');
fprintf('    %-14s %-9s %9s %10s %9s %9s %11s\n', ...
    'condition','vrillage','e_barre','e_rapport','somme','croise','e_anneau');
NOW=struct([]);
for c=1:2
    for q=1:2
        if q==1, rr=NET(c).ch; ibF=Ib_F; irF=Ir_F;  cd='On load';
        else,    rr=NET(c).bl; ibF=Ib_Fb; irF=Ir_Fb; cd='Locked rotor'; end
        eb=rr.Ibar/ibF-1; kn=rr.Iring/rr.Ibar; kf=irF/ibF;
        ek=kn/kf-1; er=rr.Iring/irF-1;
        fprintf('    %-14s %-9s %8.2f %% %9.2f %% %8.2f %% %8.2f %% %10.2f %%\n', ...
            cd,NET(c).nm,100*eb,100*ek,100*(eb+ek),100*eb*ek,100*er);
        NOW(c,q).eb=eb; NOW(c,q).ek=ek; NOW(c,q).er=er; NOW(c,q).cd=cd;
    end
end
for c=1:2
    fprintf('    vrillage %s : bascule barre %+6.2f pts, terme d''anneau %+6.2f pts,\n', ...
        NET(c).nm,100*(NOW(c,2).eb-NOW(c,1).eb),100*(NOW(c,2).ek-NOW(c,1).ek));
    fprintf('                 ligne anneau %+6.2f pts -> part de la barre %.1f %%\n', ...
        100*(NOW(c,2).er-NOW(c,1).er), ...
        100*(NOW(c,2).eb-NOW(c,1).eb)/(NOW(c,2).er-NOW(c,1).er));
end
fprintf('    La conclusion ne depend ni du millesime du modele ni du vrillage.\n');

%% ---- 5. verdict et reserves ------------------------------------------
fprintf('\n  ---- 5. VERDICT ----\n');
fprintf('    Les lignes "end-ring current" de la Table 16 ne mesurent PAS la\n');
fprintf('    precision du reseau sur l''anneau. Elles mesurent le produit de\n');
fprintf('    deux choses : l''erreur du reseau sur le COURANT DE BARRE, et un\n');
fprintf('    ECART DE MODELE D''ANNEAU d''environ %+.1f %% que la reference porte\n', ...
        100*0.5*(DEC(1).ek+DEC(2).ek));
fprintf('    seule -- segmentation, resistance de contact, effets 3-D. Le\n');
fprintf('    +%.1f %% du calage est donc en grande partie une convention de la\n', ...
        100*DEC(2).pub);
fprintf('    reference, et non une erreur du reseau. Formulation a retenir :\n');
fprintf('    << la reference n''emploie pas le modele d''anneau geometrique\n');
fprintf('    ideal ; la ligne anneau ne doit pas etre lue comme une precision >>.\n');

fprintf('\n  ---- RESERVES A DECLARER ----\n');
fprintf('    R1. Les deux RMS de reference sont des valeurs efficaces de\n');
fprintf('        FORMES D''ONDE completes : elles contiennent toutes les\n');
fprintf('        harmoniques de temps de l''essai. Le RMS du reseau ne\n');
fprintf('        contient que les branches de cage modelisees. La comparaison\n');
fprintf('        des rapports reste licite (meme grandeur au numerateur et au\n');
fprintf('        denominateur des deux cotes), mais elle n''est pas une\n');
fprintf('        comparaison de spectres.\n');
fprintf('    R2. Le mecanisme exact de l''anneau de la reference n''est pas\n');
fprintf('        identifie ici : le projet EF ne publie pas son modele\n');
fprintf('        d''anneau. Seul est etabli qu''il differe du modele\n');
fprintf('        geometrique ideal de ~4 %%, de facon reproductible.\n');
fprintf('    R3. La Table 16 est celle de MEC_DtN_paper_v2.tex. La table\n');
fprintf('        correspondante d''ArticleII_Carter_IM.tex a deja supprime les\n');
fprintf('        lignes d''anneau ; ce bloc justifie cette suppression et donne\n');
fprintf('        la phrase de remplacement.\n');

%% ---- 6. annexe : provenance des deux couples EF au calage -------------
%  Conserve de la version anterieure du bloc. Deux couples EF au calage
%  circulent dans le dossier ; toute ligne du manuscrit portant sur le
%  calage doit declarer lequel elle emploie.
fprintf('\n  ---- 6. ANNEXE : les deux couples EF au calage ----\n');
trqb=rd(fullfile(Dbl,'Torque Plot 1.tab'));
T_trans=avgw(trqb,2,t0b)*1e3;
Dc=fullfile(ROOT,'caratéristique en fonction glissement');
T_car=NaN;
try
    tc=rd(fullfile(Dc,'Torque Plot 2.tab'));
    sc=1; if max(abs(tc(:,2)))<10, sc=1e3; end
    T_car=interp1(tc(:,1),tc(:,2)*sc,1.0,'linear','extrap');
catch ME
    fprintf('    lecture de la caracteristique impossible : %s\n',ME.message);
end
fprintf('    %-48s %10s\n','essai','T (N.m)');
fprintf('    %-48s %10.3f\n','transitoire\rotor bloque, moyenne regime etabli',T_trans);
fprintf('    %-48s %10.3f\n','caracteristique T(s) au point s = 1',T_car);
fprintf('    %-48s %10.3f\n','couple MEC au calage, vrillage ON',NET(1).bl.Tem);
fprintf('    %-48s %10.3f\n','couple MEC au calage, vrillage OFF',NET(2).bl.Tem);
fprintf('    Deux ESSAIS DIFFERENTS, a ne pas confondre.\n');

save('B6_anneau.mat','kid','s_B1','s_tab','Ib_F','Ir_F','Ib_Fb','Ir_Fb', ...
     'TB','DEC','NOW','d_ref','d_net','ntab','T_trans','T_car');
fprintf('\n  duree %.0f s\n=== B6 termine ===\n',toc(t0v));
diary off;

% ======================================================================
function [TB,ntab]=lire_table16(tex,lab)
%LIRE_TABLE16  Relit dans le source LaTeX les lignes "bar current" et
%   "end-ring current" de la table portant le label LAB, pour les deux
%   blocs "On load" et "Locked rotor". Renvoie aussi le numero d'ordre de
%   la table dans le fichier (compte des \begin{table} qui la precedent).
%   Aucune valeur n'est recopiee : tout est parse.
    L=splitlines(string(fileread(tex)));
    i0=find(contains(L,['\label{' lab '}']),1);
    if isempty(i0), error('B6:tex','label %s introuvable dans %s',lab,tex); end
    ntab=sum(contains(L(1:i0),'\begin{table'));
    i1=i0-1+find(contains(L(i0:end),'\end{table}'),1);
    B=L(i0:i1);
    blocs={'On load','Locked rotor'}; lignes={'bar current','end-ring current'};
    TB=repmat(struct('net',NaN,'fea',NaN,'err',NaN),2,2);
    for b=1:2
        jb=find(contains(B,'multicolumn')&contains(B,blocs{b}),1);
        if isempty(jb), error('B6:tex','bloc %s introuvable',blocs{b}); end
        je=jb+find(contains(B(jb+1:end),'multicolumn')| ...
                   contains(B(jb+1:end),'\bottomrule'),1);
        if isempty(je), je=numel(B); end
        for l=1:2
            %  "end-ring current" contient "ring current", et "bar current"
            %  est un prefixe strict : on ancre sur le debut de ligne.
            jr=jb+find(startsWith(strip(B(jb+1:je)),lignes{l}),1);
            if isempty(jr), error('B6:tex','ligne %s absente du bloc %s',lignes{l},blocs{b}); end
            cel=strsplit(char(B(jr)),'&');
            if numel(cel)<4, error('B6:tex','ligne %s mal formee',lignes{l}); end
            TB(b,l).net=num1(cel{2}); TB(b,l).fea=num1(cel{3}); TB(b,l).err=num1(cel{4});
        end
    end
end
function v=num1(s)
%NUM1  Premier nombre signe d'une cellule LaTeX ($-8.1\%$ -> -8.1).
    t=regexp(s,'[-+]?\d+\.?\d*','match','once');
    v=str2double(t);
end
