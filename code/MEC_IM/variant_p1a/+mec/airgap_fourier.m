function AF = airgap_fourier(ths, dths, thr, dthr, Rs, Rr, L, Nh, basis)
%AIRGAP_FOURIER (VARIANTE p1a, HORS DU PAQUET DE PRODUCTION) -- 17 septembre 2026.
%   Copie de +mec/airgap_fourier.m (production, inchange) a laquelle est ajoutee
%   la base 'p1a' (chapeau ASYMETRIQUE, eq. (6) du manuscrit, meme formule que
%   dsop.py l. 153-157) :
%       h_l, h_r = distances aux centres des colonnes voisines (grille ordonnee
%                  le long de l'alesage, periodique 2*pi),
%       T_k(n)   = [ (1 - e^{-i n h_r})/(n^2 h_r) + (1 - e^{+i n h_l})/(n^2 h_l) ] e^{-i n th_k} / pi,
%       Wc = Re T, Ws = -Im T, poids homopolaire w0 = (h_l + h_r)/2 / (2 pi).
%   Pour 'p0' et 'p1' le code est celui de la production : matrices Y identiques
%   bit a bit (MES_R3_VALIDATE, 17/09/2026 : max|dY| = 0) ; pour 'p1a' l'operateur
%   fin est identique a dsop.assemble a 1e-15 pres (meme validation).
%   Ce fichier vit dans variant_p1a/+mec, un dossier de paquet SEPARE que
%   RUN_Z9_BASIS_P1A place en tete du path : il masque la fonction de production
%   sans la modifier. Le marqueur AF.variant_p1a = true en atteste dans la sortie.
%   Ce qu'il faudrait changer pour que la production accepte 'p1a' : ajouter le
%   cas 'p1a' de proj_basis ci-dessous (12 lignes) a +mec/airgap_fourier.m.
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
[Wc_s,Ws_s,w0s] = proj_basis(lower(basis), ths, dths, n);
[Wc_r,Ws_r,w0r] = proj_basis(lower(basis), thr, dthr, n);

% ---- matrice d'admittance ----
KC = Kn.*Cn;  KD = Kn.*Dn;
Yss = Wc_s.'*(KC.*Wc_s) + Ws_s.'*(KC.*Ws_s);
Yrr = Wc_r.'*(KC.*Wc_r) + Ws_r.'*(KC.*Ws_r);
Ysr = -(Wc_s.'*(KD.*Wc_r) + Ws_s.'*(KD.*Ws_r));
% mode n=0 (moyennes)
Lam0 = mu0*2*pi*L/X;
Y = [Yss, Ysr; Ysr.', Yrr] + Lam0*([w0s, -w0r].'*[w0s, -w0r]);

AF.Y=Y; AF.n=n; AF.Cn=Cn; AF.Dn=Dn; AF.X=X; AF.Lam0=Lam0;
AF.Ms=Ms; AF.Mr=Mr; AF.Rs=Rs; AF.Rr=Rr; AF.L=L; AF.mu0=mu0;
AF.Wc_s=Wc_s; AF.Ws_s=Ws_s; AF.Wc_r=Wc_r; AF.Ws_r=Ws_r;
AF.w0s=w0s; AF.w0r=w0r; AF.basis=lower(basis);
AF.variant_p1a=true;                   % marqueur : variante p1a hors production

AF.amp    = @(Us,Ur) amplitudes(AF,Us,Ur);
AF.field  = @(Us,Ur,r,thq) field_at(AF,Us,Ur,r,thq);
AF.torque = @(Us,Ur) torque_mst(AF,Us,Ur);
end

% ======================================================================
function [Wc,Ws,w0] = proj_basis(basis, th, dth, n)
switch basis
case 'p0'
    proj = (2./(n*pi)).*sin(n*dth/2);
    Wc = proj.*cos(n*th); Ws = proj.*sin(n*th); w0 = dth/(2*pi);
case 'p1'
    proj = (4./(pi*(n.^2).*dth)).*sin(n.*dth/2).^2;
    Wc = proj.*cos(n*th); Ws = proj.*sin(n*th); w0 = dth/(2*pi);
case 'p1a'
    thp = circshift(th,1); thn = circshift(th,-1);
    hl = mod(th - thp, 2*pi); hr = mod(thn - th, 2*pi);       % distances aux centres voisins
    Tk = (1 - exp(-1i*n*hr))./((n.^2).*hr) + (1 - exp(1i*n*hl))./((n.^2).*hl);
    T  = Tk.*exp(-1i*n*th)/pi;
    Wc = real(T); Ws = -imag(T); w0 = 0.5*(hl+hr)/(2*pi);
otherwise
    error('mec:airgap_fourier:basis','base inconnue : %s',basis);
end
end

function A = amplitudes(AF,Us,Ur)
    Us=Us(:); Ur=Ur(:);
    A.as = AF.Wc_s*Us;  A.bs = AF.Ws_s*Us;
    A.ar = AF.Wc_r*Ur;  A.br = AF.Ws_r*Ur;
    A.m0s = AF.w0s*Us;  A.m0r = AF.w0r*Ur;
end

function [Br,Bt,H] = field_at(AF,Us,Ur,r,thq)
    A = amplitudes(AF,Us,Ur);
    n=AF.n; X=AF.X; mu0=AF.mu0; thq=thq(:).';
    sh = sinh(n*X);
    fs = sinh(n*log(r/AF.Rr))./sh;
    fr = sinh(n*log(AF.Rs/r))./sh;
    gs = cosh(n*log(r/AF.Rr))./sh;
    gr = cosh(n*log(AF.Rs/r))./sh;
    H.Brc = -(mu0*n/r).*(A.as.*gs - A.ar.*gr);
    H.Brs = -(mu0*n/r).*(A.bs.*gs - A.br.*gr);
    H.Btc = -(mu0*n/r).*(A.bs.*fs + A.br.*fr);
    H.Bts = +(mu0*n/r).*(A.as.*fs + A.ar.*fr);
    Br = (H.Brc.'*cos(n*thq) + H.Brs.'*sin(n*thq)).';
    Bt = (H.Btc.'*cos(n*thq) + H.Bts.'*sin(n*thq)).';
    Br = Br - mu0*(A.m0s - A.m0r)/(r*X);
end

function T = torque_mst(AF,Us,Ur)
    r = 0.5*(AF.Rs+AF.Rr);
    [~,~,H] = field_at(AF,Us,Ur,r,0);
    T = (pi*r^2*AF.L/AF.mu0) * sum(H.Brc.*H.Btc + H.Brs.*H.Bts);
end
