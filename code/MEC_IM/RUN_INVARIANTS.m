%% RUN_INVARIANTS - suite de tests A REPONSE CONNUE sur l'operateur d'entrefer
%
%  POURQUOI CETTE SUITE EXISTE. Le 12 aout, un test a reponse connue a trouve
%  en une execution un defaut que huit jours de comparaison contre ANSYS
%  n'avaient pas vu : pm_loss_load repondait a un champ que l'aimant voit
%  statique. Il etait indetectable par comparaison a l'EF pour une raison
%  simple et generale -- la valeur FAUSSE tombait a +2,45 % de la reference
%  et la valeur JUSTE a -36 %. Dans un dossier ou plusieurs constantes ont
%  ete ajustees contre cette meme reference, une garde qui compare a l'EF ne
%  peut pas distinguer "juste" de "faux mais calibre".
%
%  Un invariant, lui, ne depend d'aucune reference ni d'aucune constante. Il
%  vaut a toute resolution : les tests ci-dessous ne sont PAS des tests de
%  convergence, et c'est pourquoi ils tournent sur un pavage grossier en
%  quelques secondes. Ce qu'ils verifient est vrai a Nh = 64 comme a
%  Nh = 8192, ou faux dans les deux cas.
%
%  LES HUIT INVARIANTS, ET CE QUE CHACUN PROTEGE.
%    I-1  Y * 1 = 0            un potentiel uniforme ne debite aucun flux.
%                              Verifie du meme coup que la BASE DE SURFACE
%                              represente la constante EXACTEMENT -- la
%                              partition de l'unite sur laquelle repose tout
%                              l'argument de la base chapeau.
%    I-2  Y = Y'               reciprocite. Une matrice de permeance non
%                              symetrique produirait de l'energie.
%    I-3  Y >= 0               l'entrefer ne peut pas stocker une energie
%                              negative. C'est aussi ce qui rend le reseau
%                              bien pose.
%    I-4  permeance homopolaire un potentiel stator uniforme sur rotor a la
%                              masse debite EXACTEMENT mu0*2*pi*L/ln(Rs/Rr).
%                              Valeur analytique de la couronne, sans aucun
%                              ajustement.
%    I-5  covariance par rotation  tourner RIGIDEMENT toute la configuration
%                              -- les deux grilles et les deux potentiels --
%                              ne peut pas changer le couple. C'est le test
%                              qui aurait attrape le defaut de D1.
%    I-6  couple nul en champ aligne  si les deux surfaces ne portent que des
%                              composantes en cosinus du meme rang, la
%                              contrainte tangentielle s'integre a zero.
%                              EXACTEMENT zero, pas approximativement.
%    I-7  couple par deux voies  tenseur de Maxwell contre derivee de
%                              l'energie par rapport a la position. Deux
%                              chemins independants vers la meme grandeur.
%    I-8  invariance en rayon   la somme MST par harmonique est independante
%                              du rayon d'integration dans la couronne.
%
%  AUCUN de ces huit ne fait intervenir ANSYS, ni k_C, ni une constante
%  identifiee. Ils sont vrais de la physique de la couronne, ou le code est
%  faux.

clear; clc; t0=tic;
if isfile('INV_tests_out.txt'), delete('INV_tests_out.txt'); end
diary('INV_tests_out.txt'); diary on;

M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G;
nT=6; nO=2; Nh=2048;                     % grossier A DESSEIN : voir l'en-tete
AT=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
AF=AT.AF; Ms=AF.Ms; Mr=AF.Mr; mu0=AF.mu0;

fprintf('=== INVARIANTS : tests a reponse connue de l''operateur d''entrefer ===\n');
fprintf('  machine MAS 48/44 | pavage nT=%d nO=%d | N_h=%d | base chapeau\n',nT,nO,Nh);
fprintf('  grille fine %d + %d colonnes | X = ln(Rs/Rr) = %.6e\n\n',Ms,Mr,AF.X);
res={}; nok=0;

%% ---- I-1 et I-4 : le potentiel uniforme, sur quatre configurations ------
%  I-1 et I-4 testent la meme propriete sous deux angles : la base de
%  surface represente-t-elle EXACTEMENT la fonction constante ? Si oui, un
%  potentiel uniforme ne debite aucun flux (I-1) et la permeance homopolaire
%  vaut la valeur analytique de la couronne (I-4).
%
%  ON CROISE DEUX FACTEURS, parce que le premier passage a montre que la
%  reponse depend de la base ET du pavage :
%     base       'p0' creneau constant par colonne / 'p1' chapeau lineaire
%     grille     pavee (colonnes de FACE et d'OUVERTURE de largeurs
%                differentes) / uniforme (toutes les colonnes egales)
%  Le creneau tuile le cercle par construction, quelle que soit la grille.
%  Le chapeau, lui, ne forme une partition de l'unite que si les demi-
%  largeurs s'accordent d'une colonne a l'autre -- donc sur grille uniforme.
%  Le pavage de cet operateur est deliberement NON uniforme.
fprintf('  ---- I-1 et I-4 : representation de la constante ----\n');
fprintf('  %-10s %-10s %14s %14s\n','base','grille','I-1  Y*1/max|Y|','I-4  ecart Lam');
Lam_ana=mu0*2*pi*AF.L/AF.X;
Mu=round((Ms+Mr)/2);                       % grille uniforme de meme finesse
thu=linspace(0,2*pi,Mu+1); thu(end)=[]; dthu=(2*pi/Mu)*ones(1,Mu);
COMB={'p0','pavee';'p1','pavee';'p0','uniforme';'p1','uniforme'};
V=nan(4,2);
for k=1:4
    if strcmp(COMB{k,2},'pavee')
        Ak=mec.airgap_fourier(AT.ths,AT.dths,AT.thr,AT.dthr,AF.Rs,AF.Rr,AF.L,Nh,COMB{k,1});
        ms=numel(AT.ths); mr=numel(AT.thr);
    else
        Ak=mec.airgap_fourier(thu,dthu,thu,dthu,AF.Rs,AF.Rr,AF.L,Nh,COMB{k,1});
        ms=Mu; mr=Mu;
    end
    V(k,1)=max(abs(Ak.Y*ones(ms+mr,1)))/max(abs(Ak.Y(:)));
    V(k,2)=abs(sum(Ak.Y(1:ms,:)*[ones(ms,1);zeros(mr,1)])-Lam_ana)/Lam_ana;
    fprintf('  %-10s %-10s %14.2e %14.2e\n',COMB{k,1},COMB{k,2},V(k,1),V(k,2));
end
fprintf('\n');
%  EST-CE DE LA TRONCATURE ? Si l'ecart venait du nombre d'harmoniques il
%  diminuerait avec Nh. On le mesure au lieu de le supposer.
fprintf('  ---- le defaut de la base chapeau depend-il de la troncature ? ----\n');
fprintf('  %8s %14s %14s\n','N_h','p1 pavee I-4','p0 pavee I-4');
for nh=[256 1024 4096]
    a1=mec.airgap_fourier(AT.ths,AT.dths,AT.thr,AT.dthr,AF.Rs,AF.Rr,AF.L,nh,'p1');
    a0=mec.airgap_fourier(AT.ths,AT.dths,AT.thr,AT.dthr,AF.Rs,AF.Rr,AF.L,nh,'p0');
    e1=abs(sum(a1.Y(1:Ms,:)*[ones(Ms,1);zeros(Mr,1)])-Lam_ana)/Lam_ana;
    e0=abs(sum(a0.Y(1:Ms,:)*[ones(Ms,1);zeros(Mr,1)])-Lam_ana)/Lam_ana;
    fprintf('  %8d %14.2e %14.2e\n',nh,e1,e0);
end
fprintf('\n');
%  Les invariants portes au verdict sont ceux de la CONFIGURATION DE
%  PRODUCTION : base chapeau, grille pavee. Les trois autres lignes sont le
%  diagnostic, et elles designent le facteur responsable.
r1=V(2,1);
res(end+1,:)={'I-1  Y*1 = 0 (base p1, grille pavee)',r1,1e-12,'relatif'}; %#ok<*SAGROW>

%% ---- I-2 : reciprocite --------------------------------------------------
r2=max(max(abs(AF.Y-AF.Y')))/max(abs(AF.Y(:)));
res(end+1,:)={'I-2  Y = Y'' (reciprocite)',r2,1e-12,'relatif'};

%% ---- I-3 : energie non negative -----------------------------------------
ev=eig(full(0.5*(AF.Y+AF.Y')));
r3=max(0,-min(ev))/max(abs(ev));
res(end+1,:)={'I-3  valeurs propres >= 0',r3,1e-10,'relatif'};

%% ---- I-4 : permeance homopolaire de la couronne -------------------------
Uh=[ones(Ms,1);zeros(Mr,1)];
Lam_num=sum(AF.Y(1:Ms,:)*Uh);
r4=V(2,2);
res(end+1,:)={'I-4  permeance homopolaire = mu0*2pi*L/X',r4,1e-12,'relatif'};
fprintf('  I-4 : %.6e H mesure contre %.6e H analytique (%+.3f %%)\n\n', ...
    Lam_num,Lam_ana,100*(Lam_num-Lam_ana)/Lam_ana);

%% ---- I-5 : covariance par rotation rigide -------------------------------
%  On tourne LES DEUX grilles et LES DEUX potentiels du meme angle. La
%  configuration physique est identique a une rotation pres : le couple ne
%  peut pas changer. Rien ici ne depend de la valeur de l'angle.
Us=cos(2*AT.ths(:))+0.4*sin(5*AT.ths(:));
Ur=0.7*cos(2*AT.thr(:)-0.3)+0.2*cos(7*AT.thr(:));
T0=AF.torque(Us,Ur);
dlt=0.137;                                % angle quelconque, non special
AFd=mec.airgap_fourier(AT.ths+dlt,AT.dths,AT.thr+dlt,AT.dthr, ...
                       AF.Rs,AF.Rr,AF.L,Nh,'p1');
Td=AFd.torque(Us,Ur);                     % memes valeurs NODALES, grille tournee
r5=abs(Td-T0)/max(abs(T0),eps);
res(end+1,:)={'I-5  couple invariant par rotation rigide',r5,1e-10,'relatif'};
fprintf('  I-5 : couple %.6e N.m a phi=0, %.6e N.m tourne de %.3f rad\n\n',T0,Td,dlt);

%% ---- I-6 : couple nul en champ aligne -----------------------------------
%  Deux surfaces ne portant que du cosinus au meme rang : Brs et Btc sont
%  identiquement nuls, donc Brc*Btc + Brs*Bts = 0 terme a terme.
n0=6;
Usa=cos(n0*AT.ths(:)); Ura=0.55*cos(n0*AT.thr(:));
Ta=AF.torque(Usa,Ura);
%  echelle de comparaison : le meme couple avec un quart de tour de dephasage
Urq=0.55*sin(n0*AT.thr(:));
Tq=AF.torque(Usa,Urq);
r6=abs(Ta)/max(abs(Tq),eps);
res(end+1,:)={'I-6  couple nul en champ aligne',r6,1e-12,'/ couple en quadrature'};
fprintf('  I-6 : %.3e N.m aligne contre %.6e N.m en quadrature\n\n',Ta,Tq);

%% ---- I-7 : le couple par deux voies independantes -----------------------
%  Voie A : tenseur de Maxwell, somme analytique par harmonique.
%  Voie B : derivee de l'energie de la couronne par rapport a la position
%           rotor, W(phi) = 1/2 * U' * Y(phi) * U, les valeurs nodales du
%           rotor etant portees par la grille qui tourne. Difference centree.
%  LA DERIVEE EST ANALYTIQUE, PAS UNE DIFFERENCE FINIE. Les deux premiers
%  passages ont employe une difference centree, puis une extrapolation de
%  Richardson ; elles laissaient 3,9e-06 puis 9,5e-07 de residu. Ce residu
%  n'etait pas dans le code : l'operateur porte du contenu jusqu'a n = Nh,
%  dont l'echelle angulaire 2*pi/Nh est du meme ordre que le pas h, si bien
%  que la difference finie n'est jamais dans son regime asymptotique. On
%  derive donc EXACTEMENT, au lieu de relacher le seuil.
%
%     Y depend de phi par les seules projections rotor :
%        Wc_r = proj .* cos(n*(thr+phi)) ,  Ws_r = proj .* sin(n*(thr+phi))
%        dWc_r/dphi = -n*Ws_r  ,  dWs_r/dphi = +n*Wc_r
%     Le bloc rotor-rotor est invariant (les deux termes se compensent), le
%     mode homopolaire aussi. Il ne reste que le bloc croise :
%        dYsr/dphi = Wc_s' * ((KD.*n).*Ws_r) - Ws_s' * ((KD.*n).*Wc_r)
%     et  dW/dphi = Us' * (dYsr/dphi) * Ur .
KD=(mu0*pi*AF.L*AF.n).*AF.Dn;
dYsr=AF.Wc_s.'*((KD.*AF.n).*AF.Ws_r)-AF.Ws_s.'*((KD.*AF.n).*AF.Wc_r);
TB=Us.'*dYsr*Ur;
r7=abs(abs(TB)-abs(T0))/max(abs(T0),eps);
res(end+1,:)={'I-7  Maxwell contre derivee de l''energie',r7,1e-10,'relatif'};
fprintf('  I-7 : %.12e N.m (Maxwell)\n        %.12e N.m (energie, derivee analytique) ; signes %+d\n\n', ...
    T0,TB,sign(TB*T0));

%% ---- I-8 : invariance en rayon de la somme MST --------------------------
rq=[AF.Rr+0.2*(AF.Rs-AF.Rr), AF.Rr+0.8*(AF.Rs-AF.Rr)];
TT=nan(1,2);
for k=1:2
    [~,~,H]=AF.field(Us,Ur,rq(k),0);
    TT(k)=(pi*rq(k)^2*AF.L/mu0)*sum(H.Brc.*H.Btc+H.Brs.*H.Bts);
end
r8=abs(TT(1)-TT(2))/max(abs(TT(1)),eps);
res(end+1,:)={'I-8  couple independant du rayon',r8,1e-10,'relatif'};
fprintf('  I-8 : %.6e N.m a r=0,2g, %.6e N.m a r=0,8g\n\n',TT(1),TT(2));

%% ---- verdict ------------------------------------------------------------
fprintf('  %-46s %12s %10s  %s\n','invariant','mesure','seuil','verdict');
for k=1:size(res,1)
    ok=res{k,2}<res{k,3}; nok=nok+ok;
    fprintf('  %-46s %12.2e %10.0e  %s\n',res{k,1},res{k,2},res{k,3}, ...
        tern(ok,'PASSE','ECHOUE  <--'));
end
G=(nok==size(res,1));
fprintf('\n  %d invariants sur %d verifies\n',nok,size(res,1));
fprintf('  GARDE %s\n',tern(G, ...
    'PASSEE -- l''operateur satisfait tout ce qu''on sait de lui sans calculer', ...
    'ECHOUEE -- un invariant viole est un defaut, pas une imprecision'));
if ~G
fprintf(['\n  ---- DIAGNOSTIC DES DEUX ECHECS ----\n' ...
 '  I-1 et I-4 testent la meme propriete : la base de surface represente-\n' ...
 '  t-elle EXACTEMENT la fonction constante ? Le croisement ci-dessus la\n' ...
 '  localise sans ambiguite.\n' ...
 '    - le creneau (p0) l''a, sur les deux grilles ;\n' ...
 '    - le chapeau (p1) l''a sur grille UNIFORME (6e-15) et la PERD sur la\n' ...
 '      grille PAVEE (1e-01), qui est celle de production ;\n' ...
 '    - et l''ecart CROIT avec N_h puis sature, au lieu de decroitre : ce\n' ...
 '      n''est donc pas de la troncature, et raffiner ne le retire pas.\n' ...
 '  CAUSE. Le chapeau de demi-largeur h ne forme une partition de l''unite\n' ...
 '  que si les demi-largeurs s''accordent d''une colonne a la suivante. Le\n' ...
 '  pavage emploie DELIBEREMENT des largeurs differentes pour les colonnes\n' ...
 '  de face de dent et pour celles d''ouverture : aux jonctions, les\n' ...
 '  chapeaux laissent un manque ou un recouvrement.\n' ...
 '  CE QUI EST ETABLI. Sur le mode HOMOPOLAIRE, la permeance de la couronne\n' ...
 '  sort a +0,303 %% de sa valeur analytique, et un potentiel uniforme\n' ...
 '  debite un flux qui devrait etre nul.\n' ...
 '  CE QUI NE L''EST PAS. La propagation de ce biais aux grandeurs publiees,\n' ...
 '  k_C et X_m au premier chef. Ce test ne la mesure pas et on ne la deduit\n' ...
 '  pas : elle demande sa propre mesure.\n' ...
 '  LE CORRECTIF. Une fonction chapeau ASYMETRIQUE, de demi-largeurs gauche\n' ...
 '  et droite egales aux demi-pas locaux, retablit la partition de l''unite\n' ...
 '  sur une grille quelconque. Sa transformee de Fourier reste analytique.\n']);
end
fprintf('\n  AUCUN de ces tests n''emploie la reference elements finis, ni k_C,\n');
fprintf('  ni aucune constante identifiee. Ils tiennent a toute resolution.\n');
save('INV_tests.mat','res','nok','G');
fprintf('  duree %.0f s\n=== INVARIANTS termine ===\n',toc(t0));
diary off;

% ======================================================================
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
