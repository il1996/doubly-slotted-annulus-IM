function AF = airgap_fourier(ths, dths, thr, dthr, Rs, Rr, L, Nh, basis)
%   BASE DE SURFACE (argument 'basis', defaut 'p0') -- tache T16.
%     'p0' : potentiel CONSTANT par colonne. Projection en 1/n, donc le terme
%            d'energie Kn*|W|^2 ~ n*(1/n)^2 = 1/n : SERIE HARMONIQUE,
%            logarithmiquement DIVERGENTE. C'est la cause mesuree de la
%            derive de X_m et de k_C avec Nh (taches T1/T2).
%     'p1' : potentiel LINEAIRE par morceaux (fonction chapeau). Le chapeau
%            est la convolution de deux creneaux, sa transformee est donc le
%            CARRE de celle du creneau : projection en 1/n^2, et
%            Kn*|W|^2 ~ n*(1/n^2)^2 = 1/n^3, ABSOLUMENT CONVERGENTE.
%            Les degres de liberte deviennent des valeurs NODALES et non des
%            moyennes de colonne ; les chapeaux forment une partition de
%            l'unite, donc un potentiel uniforme reste uniforme.
%AIRGAP_FOURIER  Couplage d'entrefer par matrice de permeance harmonique (P1).
%
%   AF = mec.airgap_fourier(ths,dths,thr,dthr,Rs,Rr,L,Nh) construit le
%   couplage EXACT de la couronne d'air Rr<=r<=Rs (solution de Laplace du
%   potentiel scalaire en coordonnees polaires) entre les potentiels de
%   surface stator (colonnes d'angles ths, largeurs dths, a r=Rs) et rotor
%   (thr, dthr, a r=Rr), constants par colonne :
%
%     * projection de Fourier EXACTE d'un potentiel constant par colonne :
%         Wc(n,c) = (2/(n*pi))*sin(n*dth_c/2)*cos(n*th_c)   (idem sin)
%     * par harmonique spatial n (equation radiale en t=ln r : Omega''=n^2*Omega) :
%         C_n = coth(n*X),  D_n = 1/sinh(n*X),  X = ln(Rs/Rr)
%       energie de la couronne (les DEUX composantes cos/sin) :
%         W_n = (mu0*pi*n*L/2)*[C_n*(|s|^2+|r|^2) - 2*D_n*(s.r)]
%       (limites verifiees : potentiel un cote -> mu0*pi*R*L/(2g) par
%       harmonique ; meme potentiel des 2 cotes -> energie PUREMENT
%       tangentielle mu0*pi*n^2*g*L/(2R) — c'est le flux zig-zag que le
%       couplage radial discret ne peut pas representer, defaut M6)
%     * mode n=0 : Lambda0 = mu0*2*pi*L/X (flux homopolaire de couronne).
%
%   L'effet Carter, la frange d'encoche, la composante TANGENTIELLE Bt et
%   l'attenuation radiale des harmoniques de denture EMERGENT de la solution :
%   aucun coefficient de Carter, aucune largeur de noyau sigma0.
%
%   Sorties :
%     AF.Y            : (Ms+Mr)x(Ms+Mr), ordre [stator ; rotor], symetrique
%                       semi-definie positive ; flux SORTANT des noeuds
%                       = AF.Y * [Us ; Ur]
%     AF.amp(Us,Ur)   : amplitudes harmoniques de surface (a_s,b_s,a_r,b_r)
%     AF.field(Us,Ur,r,thq) : [Br,Bt] au rayon r aux angles thq
%     AF.torque(Us,Ur): couple par tenseur de Maxwell ANALYTIQUE par
%                       harmonique, T = (pi*r^2*L/mu0)*sum_n(Brc*Btc+Brs*Bts)
%                       (independant de r dans la couronne — exact)
%     AF.n, AF.Cn, AF.Dn, AF.X, AF.Lam0
%
%   Voir aussi : mec.mesh_refined (M.opt.gap_fourier), RUN_FOURIER_GAP.

mu0 = 4*pi*1e-7;
ths=ths(:).'; dths=dths(:).'; thr=thr(:).'; dthr=dthr(:).';
Ms=numel(ths); Mr=numel(thr);
X = log(Rs/Rr);
n = (1:Nh).';

% ---- permeances harmoniques de la couronne ----
Cn = coth(n*X);
Dn = 1./sinh(n*X);
Kn = mu0*pi*L*n;                       % facteur commun (les 2 composantes)

% ---- projections de la base de surface (integrales EXACTES) ----
if nargin<9 || isempty(basis), basis='p0'; end
switch lower(basis)
case 'p0'
    %  creneau de largeur dth centre sur th : a_n = (2/(n*pi))*sin(n*dth/2)
    proj = @(dth) (2./(n*pi)).*sin(n*dth/2);
    w0   = @(dth) dth/(2*pi);
case 'p1'
    %  chapeau de DEMI-largeur h : a_n = (4/(pi*n^2*h))*sin^2(n*h/2).
    %  Verification de coherence : c'est bien le carre du noyau du creneau
    %  divise par h, comme l'exige la convolution de deux creneaux.
    %  Aire du chapeau = h, d'ou le mode homopolaire w0 = h/(2*pi).
    proj = @(h) (4./(pi*(n.^2).*h)).*sin(n.*h/2).^2;
    w0   = @(h) h/(2*pi);
otherwise
    error('mec:airgap_fourier:basis','base inconnue : %s',basis);
end
Wc_s = proj(dths).*cos(n*ths);   % (Nh x Ms)
Ws_s = proj(dths).*sin(n*ths);
Wc_r = proj(dthr).*cos(n*thr);
Ws_r = proj(dthr).*sin(n*thr);

% ---- matrice d'admittance ----
KC = Kn.*Cn;  KD = Kn.*Dn;
Yss = Wc_s.'*(KC.*Wc_s) + Ws_s.'*(KC.*Ws_s);
Yrr = Wc_r.'*(KC.*Wc_r) + Ws_r.'*(KC.*Ws_r);
Ysr = -(Wc_s.'*(KD.*Wc_r) + Ws_s.'*(KD.*Ws_r));
% mode n=0 (moyennes)
Lam0 = mu0*2*pi*L/X;
w0s = w0(dths);  w0r = w0(dthr);
Y = [Yss, Ysr; Ysr.', Yrr] + Lam0*([w0s, -w0r].'*[w0s, -w0r]);

AF.Y=Y; AF.n=n; AF.Cn=Cn; AF.Dn=Dn; AF.X=X; AF.Lam0=Lam0;
AF.Ms=Ms; AF.Mr=Mr; AF.Rs=Rs; AF.Rr=Rr; AF.L=L; AF.mu0=mu0;
AF.Wc_s=Wc_s; AF.Ws_s=Ws_s; AF.Wc_r=Wc_r; AF.Ws_r=Ws_r;
AF.w0s=w0s; AF.w0r=w0r; AF.basis=lower(basis);

AF.amp    = @(Us,Ur) amplitudes(AF,Us,Ur);
AF.field  = @(Us,Ur,r,thq) field_at(AF,Us,Ur,r,thq);
AF.torque = @(Us,Ur) torque_mst(AF,Us,Ur);
end

% ======================================================================
function A = amplitudes(AF,Us,Ur)
    Us=Us(:); Ur=Ur(:);
    A.as = AF.Wc_s*Us;  A.bs = AF.Ws_s*Us;
    A.ar = AF.Wc_r*Ur;  A.br = AF.Ws_r*Ur;
    A.m0s = AF.w0s*Us;  A.m0r = AF.w0r*Ur;
end

function [Br,Bt,H] = field_at(AF,Us,Ur,r,thq)
%  Champ dans la couronne au rayon r : solution exacte par harmonique.
    A = amplitudes(AF,Us,Ur);
    n=AF.n; X=AF.X; mu0=AF.mu0; thq=thq(:).';
    sh = sinh(n*X);
    fs = sinh(n*log(r/AF.Rr))./sh;      % poids du potentiel stator
    fr = sinh(n*log(AF.Rs/r))./sh;      % poids du potentiel rotor
    gs = cosh(n*log(r/AF.Rr))./sh;
    gr = cosh(n*log(AF.Rs/r))./sh;
    % amplitudes de champ par harmonique (H sert au couple/spectres)
    H.Brc = -(mu0*n/r).*(A.as.*gs - A.ar.*gr);
    H.Brs = -(mu0*n/r).*(A.bs.*gs - A.br.*gr);
    H.Btc = -(mu0*n/r).*(A.bs.*fs + A.br.*fr);
    H.Bts = +(mu0*n/r).*(A.as.*fs + A.ar.*fr);
    Br = (H.Brc.'*cos(n*thq) + H.Brs.'*sin(n*thq)).';
    Bt = (H.Btc.'*cos(n*thq) + H.Bts.'*sin(n*thq)).';
    % mode n=0 (radial uniforme)
    Br = Br - mu0*(A.m0s - A.m0r)/(r*X);
end

function T = torque_mst(AF,Us,Ur)
%  Tenseur de Maxwell integre sur un cercle de la couronne — la somme par
%  harmonique est INDEPENDANTE de r (verifiable), on prend le mi-entrefer.
    r = 0.5*(AF.Rs+AF.Rr);
    [~,~,H] = field_at(AF,Us,Ur,r,0);
    T = (pi*r^2*AF.L/AF.mu0) * sum(H.Brc.*H.Btc + H.Brs.*H.Bts);
end
