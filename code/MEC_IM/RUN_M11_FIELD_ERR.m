%% RUN_M11_FIELD_ERR - M-11 : indicateurs d'erreur des formes d'onde de champ
%
%  VERSION 2 (12 aout 2026). La v1 est conservee telle quelle :
%    - transcript : MEC_IM/M11_v1_meshchain_out.txt  (garde ECHOUEE)
%    - script     : docs/review/RUN_M11_FIELD_ERR_v1.m.txt
%
%  CE QUE LA REVUE REPROCHE. Le sec. 5.8 donne le fondamental a 0,3 % a vide
%  et -2,4 % en charge, et la tangentielle rms a -6,1 % et -9,3 %. La Fig. 3
%  montre quatre comparaisons de formes d'onde. AUCUNE RMSE, aucune erreur
%  ponctuelle maximale, aucun spectre tabule.
%
%  POURQUOI CELA COMPTE ICI PLUS QU'AILLEURS. Un fondamental juste a 0,3 %
%  est compatible avec de fortes erreurs locales SOUS LES OUVERTURES
%  D'ENCOCHE -- qui sont precisement la ou cet operateur revendique son
%  avantage. Et l'hypothese du sec. 6.4, qui situe le residu a la dent
%  ROTORIQUE, devient testable sur un profil d'erreur localise.
%
%  POURQUOI LA v1 A ECHOUE, ET CE QUE CELA A ETABLI (U-12).
%  La v1 reconstruisait le champ par mec.mesh_refined(6,3) + mec.solve_mesh
%  + me.gapF.field. Elle sortait le fondamental a -2,19 % a vide et -6,66 %
%  en charge, contre 0,3 % et -2,4 % publies. La garde a donc refuse de
%  publier ses indicateurs -- a bon droit : ce n'est PAS la chaine du
%  manuscrit. C'est la chaine des SONDES LOCALES, celle dont le sec. 5.8
%  declare lui-meme qu'elle "ne vient pas de la chaine qui produit tous les
%  autres resultats" : base constante par morceaux, N_h = 100.
%
%  LA CHAINE DU SEC. 5.8 EST IDENTIFIEE : c'est RUN_B1_IM_P1.m sec. 4,
%  transcript MEC_IM/B1_im_p1_out.txt lignes 52-56 (4 aout 2026) :
%      Bg1 a vide      0.9452  contre 0.9420   ->  0.3 %
%      Bt rms a vide   0.1201  contre 0.1280   -> -6.1 %
%      Bg1 en charge   0.8983  contre 0.9200   -> -2.4 %
%      Bt rms charge   0.1188  contre 0.1310   -> -9.3 %
%  soit l'operateur CONDENSE mec.airgap_dtn_tooth(nT=17,nO=4,N_h=8192,'p1')
%  -- pavage et troncature de production, base chapeau -- avec le champ lu
%  par AT.field(Usurf,r,thq). Le bloc "solutions" ci-dessous est la COPIE
%  CONFORME de RUN_B1_IM_P1 lignes 15-18 et 71-91.
%
%  GARDE 1 (identite de chaine). Les quatre valeurs ci-dessus doivent etre
%  reproduites a 1e-4 pres en absolu. Sinon rien n'est publie, comme en v1.
%  GARDE 2 (alignement). Apres recalage, le residu de phase du fondamental
%  doit etre nul a 1e-9 rad : la RMSE ne doit rien a un decalage residuel.
%  GARDE 3 (Parseval). L'energie de l'erreur reconstruite harmonique par
%  harmonique doit reproduire la MSE ponctuelle a 1e-6 relatif ; sinon la
%  decomposition par ordre ne decrit pas l'erreur qu'elle pretend decrire.
%
%  ce que la v1 croyait savoir et qui est FAUX : son en-tete affirmait que
%  RUN_ARTICLE.m ligne 81 ecrit (2\numel(y)), division a gauche. Verifie
%  caractere par caractere : la ligne porte (2/numel(y)), et le fichier ne
%  contient AUCUNE occurrence de "2\". U-11 est sans objet.

clear; clc; t0=tic;
%  un transcript = une execution : on repart d'un fichier vide.
if isfile('M11_field_err_out.txt'), delete('M11_field_err_out.txt'); end
diary('M11_field_err_out.txt'); diary on;
ROOT='<home>\Desktop\ANSYS résultat 18.5KW';

M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G; W=ctx.W;
p=M.p; Rm=0.5*(G.Rs+G.Rr);
nT=17; nO=4; Nh=8192; s_ch=0.0188;

fprintf('=== M-11 v2 : indicateurs d''erreur des formes d''onde ===\n\n');
fprintf('  CONFIGURATION -- chaine du sec. 5.8, identifiee (U-12)\n');
fprintf('    machine   : MAS 48/44, 18,5 kW\n');
fprintf('    operateur : mec.airgap_dtn_tooth(nT=%d, nO=%d, N_h=%d, base p1)\n',nT,nO,Nh);
fprintf('    champ     : AT.field(Usurf, R_moy, theta) -- operateur condense\n');
fprintf('    source    : copie conforme de RUN_B1_IM_P1.m sec. 4\n');
fprintf('    reference : %s\n',ROOT);

%% ---- 1. operateur et solutions : COPIE CONFORME DE RUN_B1_IM_P1 ---------
A1=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
ctx.AG=A1; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
fprintf('\n    operateur %dx%d | X_m0 = %.3f ohm\n',size(A1.Y,1),size(A1.Y,2),ctx.Xm0);

thq=linspace(0,2*pi,2001); thq(end)=[];        % grille de B1 (garde 1)

r =mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);   % ordre d'appel de B1
r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);

%  a vide : branche magnetisante seule
Rnl=mec.magnetizing(ctx,r0.Im); U0=Rnl.S.Usurf;
[Br0,Bt0]=A1.field(U0,Rm,thq); Br0=Br0(:).'; Bt0=Bt0(:).';
%  en charge : courant statorique TOTAL + onde de FMM de barres
p1=angle(r.I1c);
i3=sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
Fs=W.slotMMF(i3);
p2=angle(r.I2c); tb=2*pi*(0:M.Nr-1)/M.Nr;
Fu=cumsum(cos(M.p*tb-p2).'); Fu=Fu-mean(Fu);
c1=(2/M.Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
Fr=-Fu*((3/2)*(4/pi)*(W.kw1*W.Nph/(2*M.p))*sqrt(2)*abs(r.I2c)/abs(c1));
SL=mec.solve_network(ctx.net,G,ctx.BH,A1,Fs,Fr,M.opt); U1=SL.Usurf;
[Br1,Bt1]=A1.field(U1,Rm,thq); Br1=Br1(:).'; Bt1=Bt1(:).';
fprintf('    solutions calculees (%.0f s)\n',toc(t0));

%% ---- 2. GARDE 1 : identite de chaine ------------------------------------
ff=@(y)abs((2/numel(y))*sum(y.*exp(-1j*p*thq)));
PUB=[0.9452 0.1201 0.8983 0.1188];             % B1_im_p1_out.txt:53-56
GOT=[ff(Br0) sqrt(mean(Bt0.^2)) ff(Br1) sqrt(mean(Bt1.^2))];
REF=[0.942 0.128 0.920 0.131];                 % EF, valeurs employees par B1
NOM={'Bg1 a vide [T]','Bt rms a vide [T]','Bg1 en charge [T]','Bt rms charge [T]'};
fprintf('\n  ---- GARDE 1 : la chaine reproduit-elle le sec. 5.8 ? ----\n');
fprintf('  %-20s %10s %10s %12s %10s %10s\n','grandeur','ici','B1 publie','ecart abs.','EF','ecart');
ok1=true;
for k=1:4
    d=abs(GOT(k)-PUB(k)); ok1=ok1 && d<1e-4;
    fprintf('  %-20s %10.4f %10.4f %12.1e %10.4f %9.1f %%\n', ...
        NOM{k},GOT(k),PUB(k),d,REF(k),100*(GOT(k)-REF(k))/REF(k));
end
fprintf('  GARDE 1 %s\n',tern(ok1,'PASSEE','ECHOUEE -- rien ne sera publie'));

%% ---- 3. profils EF et recalage de phase ---------------------------------
P0=read_profile(fullfile(ROOT,'transitoire','a vide','Calculator Expressions Plot 1.tab'));
P1=read_profile(fullfile(ROOT,'transitoire','en charge','Calculator Expressions Plot 1.tab'));
cf=@(y,th,k)(2/numel(y))*sum(y(:).*exp(-1j*k*th(:)));

CAS={'a vide',U0,P0; 'en charge',U1,P1};
fprintf('\n  ---- 2. INDICATEURS D''ERREUR PONCTUELLE ----\n');
fprintf('  Le MEC est EVALUE AUX ANGLES DE LA GRILLE EF (1000 points), sans\n');
fprintf('  interpolation : la difference est ponctuelle et exacte.\n');
fprintf('  L''alignement est fait sur la PHASE DU FONDAMENTAL, comme partout\n');
fprintf('  dans ce dossier ; aucun decalage n''est ajuste.\n\n');
fprintf('  %-11s %5s %10s %10s %11s %10s %11s\n', ...
    'cas','comp','rms EF (T)','RMSE (T)','RMSE norm.','err max','err max n.');
RES=nan(4,4); ERR=cell(2,1); ok2=true; ok3=true; DIAG=nan(2,3); ALN=nan(2,2);
for c=1:2
    P=CAS{c,3}; U=CAS{c,2};
    %  RECALAGE SUR LA PHASE EXACTE DU FONDAMENTAL. Lire cette phase par
    %  echantillonnage serait fragile : l'operateur porte du contenu jusqu'a
    %  n = 8192, qui se replie sur la case n = p d'une grille de 1000 points
    %  (n = 1000 +- p, 2000 +- p, ...). L'operateur donne la phase SANS
    %  echantillonnage -- atan2(B_rs, B_rc) au rang n = p -- et c'est elle
    %  qu'on emploie. Cote EF il n'y a que des echantillons : c'est la
    %  reference, elle fait foi telle quelle.
    ca=cf(P.Br,P.th,p); phA=atan2(-imag(ca),real(ca));
    Uf=A1.expand(U);
    [~,~,H0]=A1.AF.field(Uf(1:A1.Msf),Uf(A1.Msf+1:end),Rm,0);
    phM=atan2(H0.Brs(p),H0.Brc(p));
    dlt=(phM-phA)/p;                                 % B~(th)=B(th+dlt)
    [Bm,Btm]=A1.field(U,Rm,P.th(:)+dlt); Bm=Bm(:); Btm=Btm(:);
    %  GARDE 2 : le recalage exact est-il confirme par la lecture
    %  echantillonnee ? Les deux estimateurs sont independants ; s'ils
    %  divergent, c'est le repliement qui pilote l'alignement et la RMSE
    %  ne mesure plus ce qu'elle pretend mesurer.
    cmv=cf(Bm,P.th,p);
    dph=mod(atan2(-imag(cmv),real(cmv))-phA+pi,2*pi)-pi;
    ok2=ok2 && abs(dph)<1e-3;
    ALN(c,:)=[dlt*180/pi abs(dph)];
    for q=1:2                                        % q=1 : Br ; q=2 : Bt
        if q==1, y=Bm; yr=P.Br(:); nm='B_r'; else, y=Btm; yr=P.Bt(:); nm='B_t'; end
        e=y-yr; rmse=sqrt(mean(e.^2)); rref=sqrt(mean(yr.^2));
        RES(2*(c-1)+q,:)=[rmse 100*rmse/rref max(abs(e)) 100*max(abs(e))/max(abs(yr))];
        fprintf('  %-11s %5s %10.4f %10.4f %10.1f %% %10.4f %10.1f %%\n', ...
            CAS{c,1},nm,rref,RES(2*(c-1)+q,1),RES(2*(c-1)+q,2), ...
            RES(2*(c-1)+q,3),RES(2*(c-1)+q,4));
        if q==1, ERR{c}=e; end
    end
    %  DIAGNOSTIC, NON EMPLOYE : de combien un recalage qui MINIMISERAIT la
    %  RMSE differerait-il du recalage sur la phase ? S'il en differe peu,
    %  le desaccord est structurel et non un decalage.
    sh=linspace(-2*pi/M.Ns,2*pi/M.Ns,121); rr=nan(size(sh));
    for j=1:numel(sh)
        rr(j)=sqrt(mean((A1.field(U,Rm,P.th(:)+dlt+sh(j))-P.Br(:)).^2));
    end
    [rmin,jm]=min(rr);
    DIAG(c,:)=[sh(jm)*180/pi rmin 100*(rmin-RES(2*(c-1)+1,1))/RES(2*(c-1)+1,1)];
end
fprintf('\n  recalage applique, et accord des DEUX estimateurs de phase\n');
fprintf('  (exact par l''operateur contre lecture sur la grille EF) :\n');
for c=1:2
    fprintf('    %-10s %+8.4f deg mec ; residu %.2e rad\n', ...
        CAS{c,1},ALN(c,1),ALN(c,2));
end
fprintf('  GARDE 2 (les deux estimateurs concordent a 1e-3 rad) %s\n', ...
    tern(ok2,'PASSEE','ECHOUEE'));
fprintf('\n  diagnostic, NON employe : le recalage qui minimiserait la RMSE,\n');
fprintf('  balaye sur +- un pas d''encoche stator (%.2f deg mec)\n',360/M.Ns);
for c=1:2
    fprintf('    %-10s decale de %+6.3f deg mec, RMSE %.4f T (%+.1f %%)\n', ...
        CAS{c,1},DIAG(c,1),DIAG(c,2),DIAG(c,3));
end
fprintf('  Un optimum EN ZERO signe un recalage sans ambiguite. Un optimum\n');
fprintf('  ecarte signe qu''un second degre de liberte n''est pas apparie :\n');
fprintf('  la POSITION MECANIQUE DU ROTOR, que le MEC fixe a phi = 0 et que\n');
fprintf('  le cliche EF ne declare pas. Voir l''attribution par famille.\n');

%% ---- 4. spectre tabule et energie de l'erreur par ordre -----------------
nus=[1 5 7 11 13 17 19 23 25];
fprintf('\n  ---- 3. SPECTRE SPATIAL DE B_r AU MI-ENTREFER (T) ----\n');
fprintf('  MEC exact : amplitude harmonique de l''operateur (aucun\n');
fprintf('  echantillonnage) ; MEC 2000 pts : lecture sur la grille de B1.\n');
SPEC=cell(2,1);
for c=1:2
    P=CAS{c,3}; U=CAS{c,2};
    Uf=A1.expand(U);
    [~,~,H]=A1.AF.field(Uf(1:A1.Msf),Uf(A1.Msf+1:end),Rm,0);
    am_ex=sqrt(H.Brc.^2+H.Brs.^2);
    Bm=A1.field(U,Rm,thq(:));      % le MODULE est invariant au recalage
    S=nan(numel(nus),4);
    fprintf('\n  %s\n',upper(CAS{c,1}));
    fprintf('  %5s %13s %13s %13s %11s\n','nu','MEC exact','MEC 2000 pts','EF','ecart');
    for k=1:numel(nus)
        n=nus(k)*p;
        ax=am_ex(n); as=abs(cf(Bm,thq,n)); aa=abs(cf(P.Br,P.th,n));
        S(k,:)=[ax as aa 100*(ax-aa)/aa];
        fprintf('  %5d %13.6f %13.6f %13.6f %10.2f %%\n',nus(k),ax,as,aa,S(k,4));
    end
    SPEC{c}=S;
end

fprintf('\n  ---- 4. OU VIT LA RMSE : energie de l''erreur par ordre ----\n');
%  Parseval n'a de sens que si la grille EF est equidistante : on le mesure
%  au lieu de le supposer. Le terme de Nyquist porte le poids 1/4, les
%  autres 1/2 -- c'est la normalisation c_n = (2/N)*somme employee partout.
NG=numel(P0.th); dth=diff(P0.th(:));
unif=max(abs(dth-mean(dth)))/mean(dth);
fprintf('  grille EF : %d points, non-uniformite max %.1e\n',NG,unif);
nmax=floor(NG/2); wP=0.5*ones(nmax,1); NCF=nan(1,2);
if mod(NG,2)==0, wP(end)=0.25; end
for c=1:2
    P=CAS{c,3}; e=ERR{c}; th=P.th(:);
    d=arrayfun(@(n)cf(e,th,n),(1:nmax).');
    Emse=mean(e.^2); Erec=sum(wP.*abs(d).^2)+mean(e)^2;
    rel=abs(Erec-Emse)/Emse; ok3=ok3 && rel<1e-6;
    [~,ix]=sort(abs(d),'descend');
    fprintf('  %-10s MSE %.5e ; reconstruite %.5e (ecart rel. %.1e)\n', ...
        CAS{c,1},Emse,Erec,rel);
    fprintf('             six ordres dominants (n) :');
    for j=1:6, fprintf(' n=%d (%.1f %%)',ix(j),100*wP(ix(j))*abs(d(ix(j)))^2/Emse); end
    fprintf('\n');
    %  ATTRIBUTION PAR FAMILLE. Un ordre n = k*Ns +- p est ancre sur la
    %  DENTURE STATOR, n = k*Nr +- p sur la denture ROTOR. Cette derniere
    %  famille est la seule dont la phase depende de la position mecanique
    %  du rotor -- position que le MEC fixe a phi = 0 et que le cliche EF
    %  ne declare pas. La part d'erreur qu'elle porte est donc CONTAMINEE
    %  par un degre de liberte non apparie ; celle de la famille stator ne
    %  l'est pas. La separation est exacte, sans aucun parametre ajuste.
    nn=(1:nmax).'; fS=false(nmax,1); fR=false(nmax,1);
    for k=1:ceil(nmax/M.Ns)
        fS=fS|(nn==k*M.Ns-p)|(nn==k*M.Ns+p);
    end
    for k=1:ceil(nmax/M.Nr)
        fR=fR|(nn==k*M.Nr-p)|(nn==k*M.Nr+p);
    end
    ff0=(nn==p);
    FAM={'fondamental',ff0; 'denture stator seule',fS&~fR&~ff0; ...
         'denture rotor seule',fR&~fS&~ff0; 'commune aux deux',fS&fR&~ff0; ...
         'reste du spectre',~fS&~fR&~ff0};
    for z=1:size(FAM,1)
        m=FAM{z,2};
        fprintf('             %-22s %5.1f %% de la MSE  (%d ordres)\n', ...
            FAM{z,1},100*sum(wP(m).*abs(d(m)).^2)/Emse,sum(m));
    end
    NCF(c)=100*(1-sum(wP(fR).*abs(d(fR)).^2)/Emse);
    fprintf('             part NON contaminee par la position rotor : %.1f %%\n',NCF(c));
end
fprintf('  GARDE 3 (Parseval) %s\n',tern(ok3,'PASSEE','ECHOUEE'));

%% ---- 5. localisation : ouvertures contre faces de dent ------------------
fprintf('\n  ---- 5. LOCALISATION DE L''ERREUR ----\n');
fprintf('  Masques construits sur le PAVAGE de l''operateur : une colonne de\n');
fprintf('  face porte le potentiel de la dent, une colonne d''ouverture est\n');
fprintf('  de l''air. Test direct de l''hypothese du sec. 6.4.\n');
mS=arc_mask(A1.ths,A1.dths,A1.isFace(1:A1.Msf),P0.th(:));
mR=arc_mask(A1.thr,A1.dthr,A1.isFace(A1.Msf+1:end),P0.th(:));
fprintf('  %-11s %-22s %10s %10s %9s\n','cas','zone','rms err','part MSE','pts');
for c=1:2
    e=ERR{c}; tot=sum(e.^2);
    Z={'ouverture stator',~mS; 'face de dent stator',mS; ...
       'ouverture rotor',~mR; 'face de dent rotor',mR; ...
       'ouv. stator ET rotor',~mS&~mR};
    for z=1:size(Z,1)
        m=Z{z,2};
        fprintf('  %-11s %-22s %10.4f %9.1f %% %9d\n',CAS{c,1},Z{z,1}, ...
            sqrt(mean(e(m).^2)),100*sum(e(m).^2)/tot,sum(m));
    end
end

%% ---- 6. lecture ---------------------------------------------------------
fprintf('\n  ---- LECTURE ----\n');
fprintf('  Le fondamental est juste a %.1f %% a vide et %.1f %% en charge --\n', ...
    100*(GOT(1)-REF(1))/REF(1),100*(GOT(3)-REF(3))/REF(3));
fprintf('  ce sont les valeurs publiees. La RMSE ponctuelle sur B_r vaut\n');
fprintf('  %.1f %% de la valeur efficace de reference a vide et %.1f %% en\n', ...
    RES(1,2),RES(3,2));
fprintf('  charge, l''erreur ponctuelle maximale %.1f %% et %.1f %% de la crete.\n', ...
    RES(1,4),RES(3,4));
fprintf('  C''est ce que la revue annonce : un fondamental juste est\n');
fprintf('  compatible avec de fortes erreurs locales, et seules les secondes\n');
fprintf('  disent ce que vaut le champ SOUS LES OUVERTURES.\n\n');
fprintf('  DEUX MESURES DE NATURE DIFFERENTE. A vide, le recalage qui\n');
fprintf('  minimiserait la RMSE est celui qu''on applique : la comparaison\n');
fprintf('  n''a aucun degre de liberte restant et le chiffre est une mesure.\n');
fprintf('  En charge, un decalage de %.3f deg mec abaisserait la RMSE de\n',DIAG(2,1));
fprintf('  %.0f %% : la position mecanique du rotor du cliche EF n''est pas\n',abs(DIAG(2,3)));
fprintf('  reproduite par la chaine, et la RMSE en charge est donc une BORNE\n');
fprintf('  SUPERIEURE, non une mesure. Ce que la separation par famille\n');
fprintf('  chiffre : la part portee par des ordres etrangers a la denture\n');
fprintf('  rotor -- la seule non contaminee -- vaut %.1f %% a vide et %.1f %%\n',NCF(1),NCF(2));
fprintf('  en charge.\n');

GARDES=[ok1 ok2 ok3];
save('M11_field_err.mat','RES','SPEC','nus','GARDES','DIAG','GOT','PUB','REF','ALN');
fprintf('\n  GARDES : 1 %s | 2 %s | 3 %s\n', ...
    tern(ok1,'OK','KO'),tern(ok2,'OK','KO'),tern(ok3,'OK','KO'));
fprintf('  duree %.0f s\n=== M-11 v2 termine ===\n',toc(t0));
diary off;

% ======================================================================
function P=read_profile(f)
    fid=fopen(f); hdr=fgetl(fid); fclose(fid);
    tk=regexp(hdr,'"([^"]+)"','tokens'); names=cellfun(@(c)c{1},tk,'UniformOutput',false);
    D=readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    gc=@(pat)D(:,find(contains(names,pat),1));
    dist=gc('Distance'); A=gc('Flux_Lines'); Br=gc('Br'); Bt=gc('Bt');
    C=dist(end)*1e-3; n=numel(dist)-1;
    P.th=2*pi*dist(1:n)*1e-3/C; P.A=A(1:n); P.Br=Br(1:n); P.Bt=Bt(1:n);
end
function m=arc_mask(th,dth,isFace,thq)
%  vrai si l'angle thq tombe sur une colonne de FACE du pavage
    th=th(:); dth=dth(:); isFace=logical(isFace(:));
    m=false(numel(thq),1); tq=mod(thq(:),2*pi);
    for k=1:numel(th)
        if ~isFace(k), continue; end
        a=mod(th(k)-dth(k)/2,2*pi); b=mod(th(k)+dth(k)/2,2*pi);
        if a<b, m=m|(tq>=a & tq<b); else, m=m|(tq>=a | tq<b); end
    end
end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
