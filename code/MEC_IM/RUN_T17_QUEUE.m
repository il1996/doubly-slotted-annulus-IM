%% RUN_T17_QUEUE - sommation analytique de la queue harmonique
%
%  T17. Pour n >> 1/X_g le noyau de la couronne degenere : Cn = coth(n*X) -> 1
%  et Dn = 1/sinh(n*X) -> 0 exponentiellement. La queue n'affecte donc que les
%  blocs DIAGONAUX Yss et Yrr, jamais le couplage Ysr.
%
%  Pour la base P0, le terme d'energie entre colonnes i et j vaut
%
%     Yss(i,j) = (4*mu0*L/pi) * SUM_n (1/n) sin(n di/2) sin(n dj/2) cos(n dth)
%
%  avec dth = th_i - th_j. Le produit de trois fonctions circulaires se
%  decompose en QUATRE cosinus :
%
%     sinA sinB cosC = 1/4 [cos(A-B+C) + cos(A-B-C) - cos(A+B+C) - cos(A+B-C)]
%
%  et chacun admet la forme fermee classique
%
%     SUM_{n>=1} cos(n*alpha)/n = -ln|2 sin(alpha/2)|
%
%  d'ou la somme INFINIE en forme fermee, sans aucune troncature :
%
%     Yss_inf(i,j) = (mu0*L/pi) * [ -ln|2sin(a1/2)| - ln|2sin(a2/2)|
%                                   +ln|2sin(a3/2)| + ln|2sin(a4/2)| ]
%     a1 = (di-dj)/2 + dth      a3 = (di+dj)/2 + dth
%     a2 = (di-dj)/2 - dth      a4 = (di+dj)/2 - dth
%
%  CE QUE CELA REVELE. Sur la DIAGONALE (i = j, dth = 0) on a a1 = a2 = 0,
%  donc -ln|2 sin 0| = +infini : le terme propre DIVERGE. La divergence
%  constatee empiriquement en T1/T2 n'est donc pas un defaut de troncature
%  mais une propriete EXACTE de la base P0 -- l'auto-energie d'un potentiel
%  discontinu est infinie. Aucun Nh ne peut la faire converger.
clear; clc;
diary('T17_out.txt'); diary on;
mu0=4*pi*1e-7;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G;
X=log(G.Rs/G.Rr); L=M.L;
fprintf('=== T17 : sommation analytique de la queue ===\n');
fprintf('  X_g = %.4e   1/X_g = %.1f\n\n',X,1/X);

S=@(a) -log(abs(2*sin(a/2))+realmin);      % SUM cos(n a)/n
clos=@(di,dj,dth) (mu0*L/pi)*( S((di-dj)/2+dth) + S((di-dj)/2-dth) ...
                             - S((di+dj)/2+dth) - S((di+dj)/2-dth) );

%% ---- 1. verification de la forme fermee HORS diagonale ---------------
%  On somme numeriquement avec Cn = 1 et on compare a la forme fermee.
d1=0.006; d2=0.006; dth=0.05;             % deux colonnes distinctes
fprintf('  1. Verification hors diagonale (di=dj=%.3f rad, dth=%.3f rad)\n',d1,dth);
fprintf('     %10s %16s %16s %12s\n','N','somme numerique','forme fermee','ecart rel.');
ref=clos(d1,d2,dth);
for N=[1e3 1e4 1e5 1e6]
    n=(1:N).';
    num=(4*mu0*L/pi)*sum((1./n).*sin(n*d1/2).*sin(n*d2/2).*cos(n*dth));
    fprintf('     %10d %16.6e %16.6e %11.2e\n',N,num,ref,abs(num-ref)/abs(ref));
end

%% ---- 2. la diagonale DIVERGE, et en 1/2 ln N -------------------------
fprintf('\n  2. Terme propre (i = j, dth = 0) : croissance en ln N\n');
fprintf('     %10s %16s %14s\n','N','Yss(i,i)','/(0.5*ln N)');
for N=[1e3 1e4 1e5 1e6]
    n=(1:N).';
    dd=(4*mu0*L/pi)*sum((1./n).*sin(n*d1/2).^2);
    fprintf('     %10d %16.6e %14.6e\n',N,dd,dd/(0.5*log(N)));
end
fprintf(['     Le rapport se stabilise : la croissance est bien\n' ...
         '     logarithmique, coefficient (4*mu0*L/pi)*(1/2) = %.4e\n'], ...
         (4*mu0*L/pi)*0.5);
fprintf('     => AUCUN Nh ne fait converger la base P0. T1 mesurait cela.\n');

%% ---- 3. base P1 : la queue converge ABSOLUMENT -----------------------
%  Chapeau : W ~ 1/n^2, donc Kn*|W|^2 ~ 1/n^3. La somme residuelle au-dela
%  de N est bornee par SUM_{n>N} 1/n^3 ~ 1/(2N^2) : elle s'annule vite.
fprintf('\n  3. Base P1 : queue residuelle au-dela de N\n');
h=d1;
fprintf('     %10s %18s %14s\n','N','queue P1 restante','x N^2');
for N=[1e2 1e3 1e4 1e5]
    n=((N+1):(200*N)).';
    q=(mu0*pi*L)*sum(n.*((4./(pi*(n.^2)*h)).*sin(n*h/2).^2).^2);
    fprintf('     %10d %18.6e %14.4f\n',N,q,q*N^2);
end
fprintf(['     Le produit queue x N^2 se stabilise : decroissance en 1/N^2,\n' ...
         '     conforme a SUM_{n>N} 1/n^3. La troncature devient un choix\n' ...
         '     de precision, non une source de biais.\n']);

fprintf(['\n  CONCLUSION T17. La forme fermee explique ANALYTIQUEMENT ce que\n' ...
         '  T1 et T2 avaient mesure : la divergence de la base P0 est portee\n' ...
         '  par le terme propre, elle est exacte et non numerique. Elle\n' ...
         '  justifie a posteriori le passage a P1 (T16), dont la queue est\n' ...
         '  bornee en 1/N^2.\n']);
fprintf('\n=== T17 termine ===\n');
diary off;
