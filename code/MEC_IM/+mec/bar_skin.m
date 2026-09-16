function S = bar_skin(M, G, s, opt)
%BAR_SKIN  Effet de peau d'une barre de forme QUELCONQUE (1D exact).
%
%   S = mec.bar_skin(M,G,s,opt) calcule les facteurs d'effet de peau kR, kX
%   d'une barre rotorique de largeur VARIABLE b(y), en résolvant le problème
%   de diffusion 1D dans l'encoche. Les formules classiques de Field-Emde
%   (utilisées par mec.cage) supposent une barre RECTANGULAIRE ; l'encoche
%   réelle est TRAPÉZOÏDALE (br1 en haut -> br2 au fond), ce qui change kR.
%
%   PHYSIQUE. À fort glissement le courant se refoule vers le HAUT de la
%   barre (côté entrefer, là où le flux de fuite embrassé est minimal). Comme
%   le haut est la partie LARGE du trapèze, le refoulement coûte MOINS cher
%   qu'en rectangulaire => kR(trapèze) < kR(rectangle). Le modèle
%   rectangulaire SURESTIME donc Rr, et donc le couple à fort glissement.
%
%   FORMULATION (exacte en 1D, toute forme b(y)) :
%       Ampère   :  H(y)*b(y) = integrale_0^y J*b dy'      -> u = H*b
%       Faraday  :  rho*dJ/dy = j*omega*mu0*H = j*omega*mu0*u/b
%       d(u)/dy  =  J*b
%   Conditions : u(0)=0 (aucun courant enlacé au fond). Le potentiel vecteur
%   est référencé à 0 au HAUT de l'encoche (y=h) : c'est là que la tension de
%   barre se mesure (le flux au-dessus est extérieur à la barre, traité à part
%   par lambda_r_tip dans mec.leakage). D'où l'impédance INTERNE
%       Z_barre = rho * J(h) * L / I ,   I = u(h)
%   (au continu J est uniforme et l'on retrouve Z = rho*L/A ; au 1er ordre on
%   retrouve exactement L_int = mu0*L*integrale(cumA^2/b)/A^2. Référencer au
%   FOND donnerait L_int/2 — erreur classique.)
%
%   opt.shape : 'trapz' (défaut, géométrie réelle) | 'rect' (validation)
%   opt.N     : nombre de points de marche (défaut 4000)
%
%   Sorties : S.kR, S.kX (rapports AC/DC), S.Z, S.Rdc, S.Ldc, S.xi.
%
%   VALIDATION : avec opt.shape='rect', kR et kX doivent retomber sur les
%   formules exactes de Field-Emde (cf. RUN_BARSKIN).
%
%   Voir aussi : mec.cage, RUN_BARSKIN.

if nargin<4, opt=struct(); end
shape = getdef(opt,'shape','trapz');
N     = getdef(opt,'N',4000);
mu0   = 4*pi*1e-7;

sig = M.al.sigma20/(1 + M.al.alpha*M.al.theta);
rho = 1/sig;
L   = M.L;
h   = M.hr1;                       % hauteur du corps de barre
w   = 2*pi*M.f*abs(s);             % pulsation rotorique

% ---- profil de largeur b(y) : y=0 au FOND, y=h en HAUT (coté entrefer) ----
y = linspace(0,h,N+1).';
switch lower(shape)
    case 'trapz'
        b = M.br2 + (M.br1 - M.br2)*(y/h);        % fond -> haut
    case 'rect'
        Aeq = 0.5*(M.br1+M.br2)*h;                % même aire totale
        b = (Aeq/h)*ones(size(y));
    otherwise
        error('shape inconnue');
end
A_tot = trapz(y,b);
S.Rdc = rho*L/A_tot;

% ---- inductance de fuite DC du corps de barre (omega -> 0) ----
% lambda_dc = integrale_0^h [ (aire sous y)/A_tot ]^2 / b(y) dy
cumA = cumtrapz(y,b);
S.Ldc = mu0*L*trapz(y, (cumA/A_tot).^2 ./ b);

if w < 1e-9
    S.kR=1; S.kX=1; S.Z=S.Rdc + 1i*w*S.Ldc; S.xi=0; return;
end

% ---- marche 1D (RK2) de (J,u) du fond vers le haut ----
dy = h/N;
J = 1+0i; u = 0+0i;                 % J(0)=1 (linéaire : mise à l'échelle après)
c = 1i*w*mu0/rho;
for k=1:N
    bk  = b(k);  bkh = 0.5*(b(k)+b(k+1));  bk1 = b(k+1);
    % k1
    dJ1 = c*u/bk;            du1 = J*bk;
    % k2 (mi-pas)
    Jm  = J + 0.5*dy*dJ1;    um = u + 0.5*dy*du1;
    dJ2 = c*um/bkh;          du2 = Jm*bkh;
    J = J + dy*dJ2;          u = u + dy*du2;   %#ok<AGROW>
end
I = u;                               % courant total de barre  ( = u(h) )
Z = rho*J*L/I;                       % rho*J(h)*L/I  (référence de A au HAUT)

S.Z  = Z;
S.kR = real(Z)/S.Rdc;
S.kX = imag(Z)/(w*S.Ldc);
delta = sqrt(2/(w*mu0*sig));
S.xi  = h/delta;
end

function v=getdef(s,f,d), if isfield(s,f)&&~isempty(s.(f)),v=s.(f);else,v=d;end, end
