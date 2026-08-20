%% RUN_X3_TRACE_TABLE - table de verification de la section 3 de l'Article I
%
%  POURQUOI CE BLOC. T17 a etabli le bon resultat, mais son panneau 2 le
%  VERIFIE MAL, et la colonne imprimee se refute elle-meme :
%
%      N        Yss(i,i)      /(0.5*ln N)
%    1e3     3.212665e-07    9.301619e-08
%    1e4     6.164801e-07    1.338669e-07
%    1e5     9.193644e-07    1.597100e-07
%    1e6     1.222923e-06    1.770362e-07
%    "Le rapport se stabilise [...] coefficient (4*mu0*L/pi)*(1/2) = 1.3183e-07"
%
%  Le rapport ne se stabilise pas : il croit de facon monotone et a DEPASSE
%  1,3183e-07 des N = 1e5. Deux raisons, toutes deux benignes pour le fond :
%
%   (i) NORMALISATION. Si Yss ~ (4*mu0*L/pi)*(ln N)/2, alors
%       Yss/(0.5*ln N) -> 4*mu0*L/pi = 2.6366e-07, et non sa moitie. La
%       valeur annoncee est celle de la PENTE en ln N, pas celle du rapport.
%
%  (ii) TERME CONSTANT. La forme fermee complete porte une constante
%       additive qui n'est pas petite :
%
%          Yss(i,i) = (2*mu0*L/pi) * [ ln N + gamma + ln|2 sin(d/2)| ] + o(1)
%
%       Avec d = 0,006 rad, gamma + ln|2 sin(d/2)| = -4,539. Diviser par
%       0,5*ln N laisse donc un terme en 1/ln N qui ne s'eteint qu'a la
%       vitesse ou ln N grandit : a N = 1e6 il vaut encore -33 %.
%
%  UN RAPPORTEUR QUI VERIFIE LA COLONNE VERRA CELA. On regenere donc le
%  panneau contre la forme fermee COMPLETE, avec :
%    - la colonne d'ecart relatif, qui doit tendre vers zero ;
%    - un controle de PENTE independant (increment par decade) ;
%    - le meme test a TROIS largeurs de colonne, pour montrer que la
%      constante suit bien ln|2 sin(d/2)| et non un ajustement heureux.
%
%  CONFIGURATION DECLAREE
%    machine    : MAS 48/44, 18,5 kW (mec.machine_18_5kW)
%    grandeur   : terme propre de l'operateur de couronne en base P0,
%                 noyau degenere Cn -> 1 (limite n >> 1/X_g)
%    X_g        : log(Rs/Rr) du contexte, imprime en tete
%    troncature : N = 1e3 a 1e7
%    base P1    : noyau chapeau, queue au-dela de N
clear; clc; t0=tic;
diary('X3_trace_table_out.txt'); diary on;
mu0=4*pi*1e-7;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G;
X=log(G.Rs/G.Rr); L=M.L; gam=0.5772156649015329;
A2=2*mu0*L/pi;                       % pente en ln N de la forme fermee
fprintf('=== X3 : table de verification de la conformite de la trace ===\n');
fprintf('  machine : MAS 48/44, 18,5 kW | L = %.4f m\n',L);
fprintf('  X_g = %.6e | 1/X_g = %.1f\n',X,1/X);
fprintf('  pente predite  2*mu0*L/pi        = %.6e\n',A2);
fprintf('  limite du rapport Yss/(0.5*ln N) = %.6e   (soit 2 x la pente)\n',2*A2);
fprintf('  constante d''Euler-Mascheroni gamma = %.10f\n\n',gam);

S=@(a) -log(abs(2*sin(a/2))+realmin);      % SUM_{n>=1} cos(n a)/n
clos=@(di,dj,dth) (mu0*L/pi)*( S((di-dj)/2+dth) + S((di-dj)/2-dth) ...
                             - S((di+dj)/2+dth) - S((di+dj)/2-dth) );

%% ---- 1. hors diagonale : la forme fermee est exacte ------------------
d1=0.006; dth=0.05;
fprintf('  ---- 1. hors diagonale (di = dj = %.3f rad, dth = %.3f rad) ----\n',d1,dth);
fprintf('  %10s %18s %18s %14s\n','N','somme numerique','forme fermee','ecart rel.');
ref=clos(d1,d1,dth);
for N=[1e3 1e4 1e5 1e6 1e7]
    n=(1:N).';
    num=(4*mu0*L/pi)*sum((1./n).*sin(n*d1/2).^2.*cos(n*dth));
    fprintf('  %10.0e %18.8e %18.8e %13.2e\n',N,num,ref,abs(num-ref)/abs(ref));
end
fprintf('  La serie hors diagonale CONVERGE : la queue n''y vit pas.\n\n');

%% ---- 2. terme propre : forme fermee COMPLETE -------------------------
fprintf('  ---- 2. terme propre (i = j, dth = 0) ----\n');
DS=[0.006 0.012 0.024];
for k=1:numel(DS)
    d=DS(k); C=gam+log(abs(2*sin(d/2)));
    fprintf('  -- largeur de colonne d = %.3f rad : gamma + ln|2 sin(d/2)| = %+.6f --\n',d,C);
    fprintf('  %10s %18s %18s %12s %16s\n', ...
        'N','Yss(i,i) somme','forme fermee','ecart rel.','increment/decade');
    prev=NaN;
    for N=[1e3 1e4 1e5 1e6 1e7]
        n=(1:N).';
        dd=(4*mu0*L/pi)*sum((1./n).*sin(n*d/2).^2);
        cf=A2*(log(N)+C);
        if isnan(prev), inc=NaN; else, inc=dd-prev; end
        prev=dd;
        if isnan(inc)
            fprintf('  %10.0e %18.8e %18.8e %11.2e %16s\n',N,dd,cf,abs(dd-cf)/abs(cf),'--');
        else
            fprintf('  %10.0e %18.8e %18.8e %11.2e %16.8e\n',N,dd,cf,abs(dd-cf)/abs(cf),inc);
        end
    end
    fprintf('  increment par decade predit : 2*mu0*L/pi * ln 10 = %.8e\n\n',A2*log(10));
end
fprintf(['  LECTURE. L''ecart a la forme fermee complete tombe sous 1e-3 des\n' ...
         '  N = 1e5 et continue de decroitre, aux TROIS largeurs de colonne.\n' ...
         '  La constante additive suit ln|2 sin(d/2)| : elle n''est pas un\n' ...
         '  ajustement, c''est le terme que la sommation produit. Le terme\n' ...
         '  propre CROIT SANS LIMITE en ln N -- aucune troncature ne le borne.\n\n']);

%% ---- 3. base P1 : la queue converge absolument -----------------------
fprintf('  ---- 3. base P1 (chapeau) : queue residuelle au-dela de N ----\n');
h=d1;
fprintf('  %10s %20s %14s\n','N','queue P1 restante','x N^2');
for N=[1e2 1e3 1e4 1e5 1e6]
    n=((N+1):(200*N)).';
    q=(mu0*pi*L)*sum(n.*((4./(pi*(n.^2)*h)).*sin(n*h/2).^2).^2);
    fprintf('  %10.0e %20.8e %14.6f\n',N,q,q*N^2);
end
fprintf(['  Le produit queue x N^2 se stabilise : decroissance en 1/N^2,\n' ...
         '  conforme a SUM_{n>N} 1/n^3. La troncature redevient un choix de\n' ...
         '  precision, non une source de biais.\n\n']);

%% ---- 4. pourquoi le PMSM ne derive pas -------------------------------
%  L'EXPLICATION PAR L'EPAISSEUR NE TIENT PAS, ET IL FAUT LE DIRE. On lit
%  parfois que la couronne epaisse du PMSM (aimant de 3,5 mm compris)
%  placerait la machine "sous le seuil" ou la queue s'installe. C'est
%  l'inverse : le noyau degenere pour n >> 1/X, et 1/X vaut 7 pour le
%  PMSM contre 337 pour la MAS. La queue du PMSM commence donc BIEN PLUS
%  TOT. Si l'epaisseur etait le mecanisme, le PMSM deriverait davantage.
%
%  LE MECANISME EST LE COUPLAGE ENTRE TRONCATURE ET PAVAGE. Le terme
%  propre vaut A2*[ln N + gamma + ln|2 sin(d/2)|]. Deux chaines :
%
%    PMSM (cogging_mec) : numax = Nsurf/2 IMPOSE (rang de Y), et les
%      colonnes sont uniformes, d = 2*pi/Nsurf. Alors
%          ln N + ln|2 sin(d/2)| ~ ln(Nsurf/2) + ln(2*pi/Nsurf) = ln(pi),
%      INDEPENDANT de Nsurf. Le terme divergent est STATIONNAIRE par
%      construction : raffiner la surface raffine la troncature d'autant,
%      et les deux logarithmes se compensent exactement.
%
%    MAS (mec.airgap_dtn_tooth) : le pavage est FIXE (nT=17, nO=4, 92
%      noeuds, d ~ 0,00626 rad) et N_h est balaye independamment. Rien ne
%      compense ln N_h : le terme propre croit.
%
%  Ce n'est donc pas une immunite, c'est un VERROUILLAGE. Quiconque
%  deverrouille numax sur le PMSM verra la derive apparaitre. C'est la
%  formulation qu'il faut publier -- et elle se verifie ici.
fprintf('  ---- 4. pourquoi le PMSM ne derive pas : troncature verrouillee ----\n');
Xim=X;
Rs_p=69.356356079956e-3/2; hm=3.5e-3; lag=1e-3; L_p=33e-3;
rmi_p=Rs_p-lag-hm; Xpm=log(Rs_p/rmi_p); A2p=2*mu0*L_p/pi;
fprintf('  epaisseurs logarithmiques :\n');
fprintf('    MAS 48/44  X = ln(Rs/Rr)   = %.6e  -> 1/X = %6.1f\n',Xim,1/Xim);
fprintf('    PMSM 15/14 X = ln(Rs/r_mi) = %.6e  -> 1/X = %6.1f\n',Xpm,1/Xpm);
fprintf('    La queue s''installe pour n >> 1/X : elle commence donc PLUS TOT\n');
fprintf('    sur le PMSM (n ~ %.0f) que sur la MAS (n ~ %.0f). L''epaisseur\n',1/Xpm,1/Xim);
fprintf('    n''explique PAS l''absence de derive -- elle la contredirait.\n\n');

fprintf('  (a) chaine PMSM : numax = Nsurf/2, colonnes d = 2*pi/Nsurf\n');
fprintf('  %8s %8s %12s %14s %18s\n','Nsurf','N','d (rad)','ln N + C(d)','terme propre');
NS=[540 1080 1260 2160 4320 8640]; Tp=nan(size(NS));
for k=1:numel(NS)
    Ns_=NS(k); Nk=floor(Ns_/2); d=2*pi/Ns_;
    br=log(Nk)+gam+log(abs(2*sin(d/2)));
    n=(1:Nk).'; Tp(k)=(4*mu0*L_p/pi)*sum((1./n).*sin(n*d/2).^2);
    fprintf('  %8d %8d %12.6e %14.6f %18.8e\n',Ns_,Nk,d,br,Tp(k));
end
fprintf('    derive du terme propre sur le balayage : %+.2f %%\n', ...
    100*(Tp(end)-Tp(1))/Tp(1));
fprintf('    (ln(pi) + gamma = %.6f -- valeur vers laquelle la colonne\n',log(pi)+gam);
fprintf('     "ln N + C(d)" tend, sans dependre de Nsurf)\n\n');

fprintf('  (b) chaine MAS : pavage FIXE nT=17, nO=4, N_h balaye seul\n');
dIM=0.00626;
fprintf('  %8s %12s %14s %18s\n','N_h','d (rad)','ln N + C(d)','terme propre');
NH=[512 1024 2048 4096 8192]; Ti=nan(size(NH));
for k=1:numel(NH)
    Nk=NH(k); br=log(Nk)+gam+log(abs(2*sin(dIM/2)));
    n=(1:Nk).'; Ti(k)=(4*mu0*L/pi)*sum((1./n).*sin(n*dIM/2).^2);
    fprintf('  %8d %12.6e %14.6f %18.8e\n',Nk,dIM,br,Ti(k));
end
fprintf('    derive du terme propre sur le balayage : %+.2f %%\n', ...
    100*(Ti(end)-Ti(1))/Ti(1));
fprintf(['\n    VERDICT. Le terme divergent est stationnaire a %+.2f %% sur la\n' ...
         '    chaine PMSM et croit de %+.0f %% sur la chaine MAS, pour la meme\n' ...
         '    base P0 et le meme operateur. La difference n''est pas dans la\n' ...
         '    machine : elle est dans le COUPLAGE entre troncature et pavage.\n'], ...
         100*(Tp(end)-Tp(1))/Tp(1),100*(Ti(end)-Ti(1))/Ti(1));

save('X3_trace_table.mat','A2','A2p','X','L','L_p','gam','DS','Xim','Xpm', ...
     'NS','Tp','NH','Ti','dIM');
fprintf('\n  duree %.0f s\n=== X3 termine ===\n',toc(t0));
diary off;
