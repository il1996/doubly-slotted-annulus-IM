function res = equivalent_circuit(ctx, s, Xm_init)
%EQUIVALENT_CIRCUIT  Point de fonctionnement complet à un glissement donné.
%
%   res = mec.equivalent_circuit(ctx, s, Xm_init) résout le schéma
%   équivalent monophasé de la machine asynchrone au glissement s, avec :
%     * réactance magnétisante Xm SATURÉE fournie par le réseau MEC
%       (point fixe MEC<->circuit sur le courant magnétisant) ;
%     * résistance de fer Rfe déduite des pertes fer LOCALES du MEC ;
%     * résistance et fuite rotoriques Rr(s),Xr(s) avec effet de peau
%       (mec.cage) ;
%     * fuites différentielles proportionnelles à Xm.
%
%   Schéma :  Uph = (Rs + jXs)*I1 + E1 ,  E1 aux bornes du parallèle
%   { jXm  ||  Rfe  ||  (Rr/s + jXr) }.
%
%   Sorties (res) : I1, I2, Im, E1, cosphi, Pin, Pag, Tem, Pout, rendement,
%   pertes détaillées, Xm/Rfe convergés, Bg1, cartes de champ (res.S).
%
%   Voir aussi : mec.magnetizing, mec.cage, mec.power_balance.

M = ctx.M; G = ctx.G;
m = M.m; Uph = M.Uph; w = M.w;
ws = w/M.p;                       % vitesse mécanique synchrone [rad/s]
Rs = ctx.Rs;
Xs_leak = ctx.Lk.Xs_leak; sigd_s = ctx.Lk.sigd_s;

if nargin < 3 || isempty(Xm_init), Xm_init = ctx.Xm0; end

% Impédance de cage au glissement s (effet de peau)
Cg = mec.cage(M, G, ctx.W, ctx.Lk, s);
Rr = Cg.Rr; Xr = Cg.Xr; sigd_r = Cg.sigd_r;

% ---- Branches HARMONIQUES (couples parasites asynchrones) -------------
%  Chaque harmonique d'espace nu voit son champ a Omega_s/nu, son glissement
%  s_nu = 1 -+ nu*(1-s), et sa propre branche rotorique. Ces branches sont EN
%  SERIE dans le circuit statorique : leur somme REMPLACE la fuite
%  différentielle forfaitaire sigd_s*Xm, qui était une réactance PURE (donc
%  sans perte ni couple). Elles apportent les pertes ET les couples parasites.
%  Report rotorique HARMONIQUE (kw_nu, anneaux en sin^2(pi*nu*p/Nr)).
H = ctx.H;  nH = numel(H.nu);
snu = zeros(nH,1); Rrn = snu; Xrn = snu; frn = snu;
for k = 1:nH
    if H.dir(k) > 0, snu(k) = 1 - H.nu(k)*(1-s);
    else,            snu(k) = 1 + H.nu(k)*(1-s); end
    Ck = mec.cage(M, G, ctx.W, ctx.Lk, max(abs(snu(k)),1e-3), H.nu(k), H.kw(k));
    Rrn(k) = Ck.Rr;  Xrn(k) = Ck.Xr;
    % part d'anneau dans la resistance rotorique de l'harmonique (P2/P3)
    srn = max(sin(pi*H.nu(k)*M.p/M.Nr)^2, 1e-9);
    Req = Ck.Rring/(2*srn);
    frn(k) = Req/(Ck.Rbar + Req);
end

% ---------------- Point fixe MEC <-> circuit --------------------------
Xm = Xm_init;  Rfe = 50*ctx.Xm0;   % init Rfe grand
tolX = M.opt.sat_tol; itmax = M.opt.sat_itmax;
Rm = [];  I1 = 0; E1v = 0; Im = 0; I2 = 0; Zhk = zeros(nH,1);
for it = 1:itmax
    Xs  = Xs_leak;                 % (la fuite différentielle est explicitée)
    Xrp = Xr + sigd_r*Xm;
    Zrot = Rr/s + 1i*Xrp;
    Yp = 1/Zrot + 1/(1i*Xm) + 1/Rfe;
    Zp = 1/Yp;
    % branches harmoniques en série
    Zh = 0;
    for k = 1:nH
        Xmn = Xm*H.sig(k);
        if Xmn <= 0 || abs(snu(k)) < 1e-6, Zhk(k)=0; continue; end
        Zhk(k) = 1/( 1/(1i*Xmn) + 1/(Rrn(k)/snu(k) + 1i*Xrn(k)) );
        Zh = Zh + Zhk(k);
    end
    Ztot = Rs + 1i*Xs + Zh + Zp;
    I1  = Uph/Ztot;
    E1v = I1*Zp;
    Im  = E1v/(1i*Xm);
    I2  = E1v/Zrot;

    % Mise à jour saturée de Xm par le réseau MEC (sur |Im|)
    Rm = mec.magnetizing(ctx, abs(Im));
    Xm_new = Rm.Xm;

    % Pertes fer locales -> Rfe
    PfeS = mec.iron_losses(M, G, Rm.S, M.f);
    Pfe  = PfeS.total;
    Rfe_new = m*abs(E1v)^2/max(Pfe, 1);

    dX = abs(Xm_new - Xm)/max(Xm,eps);
    dR = abs(Rfe_new - Rfe)/max(Rfe,eps);
    Xm  = Xm + 0.6*(Xm_new - Xm);
    Rfe = Rfe + 0.6*(Rfe_new - Rfe);
    if dX < tolX && dR < 1e-2
        break;
    end
end

% ---------------- Grandeurs de sortie ---------------------------------
I1m = abs(I1); I2m = abs(I2); Imm = abs(Im);
cosphi = cos(angle(I1));                 % Uph pris comme référence réelle
Pin  = m*Uph*real(I1);                   % puissance active absorbée
Qin  = m*Uph*(-imag(I1));                % puissance réactive absorbée
Pag  = m*I2m^2*Rr/s;                     % puissance d'entrefer FONDAMENTALE
Tem1 = Pag/ws;                           % couple fondamental

% ---- Harmoniques : puissances d'entrefer, couples parasites, pertes ----
%  T_nu = Pag_nu * nu*p/omega * sens   (le facteur nu est décisif)
%  Verification de cohérence : T_nu*Omega_mec = Pag_nu*(1-s_nu) — exacte.
Pag_h = zeros(nH,1); T_h = zeros(nH,1);
for k = 1:nH
    Pag_h(k) = m*I1m^2*real(Zhk(k));
    T_h(k)   = Pag_h(k)*H.nu(k)*M.p/M.w*H.dir(k);
end
Tem_h   = sum(T_h);                      % couple parasite total (usuellement < 0)
Pcu_r_h = sum(Pag_h.*snu);               % pertes rotoriques harmoniques
Pmech_h = sum(Pag_h.*(1-snu));           % puissance mécanique harmonique

Tem  = Tem1 + Tem_h;                     % couple NET (fondamental + parasites)
Pcu_s = m*Rs*I1m^2;                      % pertes Joule stator
Pcu_r = m*Rr*I2m^2 + Pcu_r_h;            % Joule rotor (fondamental + harmoniques)
Pfe   = m*abs(E1v)^2/Rfe;                % pertes fer
Pmech_gross = Pag*(1-s) + Pmech_h;       % puissance mécanique brute
Pfw   = mec.mech_losses(M, s);           % pertes mécaniques
Padd  = mec.stray_losses(M, I2m);        % pertes supplémentaires en charge (loi I2^2)
Pout  = Pmech_gross - Pfw - Padd;         % puissance utile
n_rpm = M.ns*60*(1-s);                   % vitesse [tr/min]
if Pin > 0
    eta = max(Pout,0)/Pin;
else
    eta = NaN;
end

% ---- P3 : courants de cage PHYSIQUES (fondamental + harmoniques) -------
%  Chaque branche harmonique porte un courant rotorique REEL dans les memes
%  barres, a la frequence |s_nu|*f : les composantes (frequences differentes)
%  s'ajoutent QUADRATIQUEMENT en RMS. Report courant par harmonique :
%      Ibar_nu  = I2_nu * 2*m*kw_nu*Nph/Nr
%      Iring_nu = Ibar_nu / (2*sin(pi*nu*p/Nr))
%  C'est ce que mesure ANSYS a VIDE (27 A/barre a s~0 : denture pure) et que
%  le seul fondamental ne peut pas produire.
Wd = ctx.W;
kbar1  = 2*m*Wd.kw1*Wd.Nph/M.Nr;
Ibar1  = I2m*kbar1;
Iring1 = Ibar1/abs(2*sin(pi*M.p/M.Nr));
sr1  = max(sin(pi*M.p/M.Nr)^2, 1e-9);
Req1 = Cg.Rring/(2*sr1);
fr1  = Req1/(Cg.Rbar + Req1);            % part d'anneau (fondamental)
Ibar_h2 = 0; Iring_h2 = 0; Pring_h = 0; I2n = zeros(nH,1);
for k = 1:nH
    if Zhk(k) == 0, continue; end
    I2n(k) = abs(I1*Zhk(k)/(Rrn(k)/snu(k) + 1i*Xrn(k)));
    Ibn = I2n(k)*2*m*abs(H.kw(k))*Wd.Nph/M.Nr;
    srn = max(abs(2*sin(pi*H.nu(k)*M.p/M.Nr)), 1e-6);
    Ibar_h2  = Ibar_h2  + Ibn^2;
    Iring_h2 = Iring_h2 + (Ibn/srn)^2;
    Pring_h  = Pring_h  + frn(k)*Pag_h(k)*snu(k);
end
res.Ibar1 = Ibar1;  res.Iring1 = Iring1;
res.Ibar  = sqrt(Ibar1^2 + Ibar_h2);     % RMS physique de barre
res.Iring = sqrt(Iring1^2 + Iring_h2);   % RMS physique d'anneau
res.I2n = I2n;
res.Pring = fr1*m*Rr*I2m^2 + Pring_h;    % pertes des anneaux (fond.+harm.)
res.Pbars = Pcu_r - res.Pring;           % pertes des barres
res.kRing = Cg.kRing;

res.s = s;  res.n_rpm = n_rpm;
res.I1 = I1m; res.I2 = I2m; res.Im = Imm;
res.I1c = I1; res.I2c = I2; res.E1 = abs(E1v); res.E1c = E1v;
res.cosphi = cosphi; res.Pin = Pin; res.Qin = Qin;
res.Pag = Pag; res.Tem = Tem;
res.Tem1 = Tem1; res.Tem_h = Tem_h;      % couple fondamental / parasite
res.Pag_h = sum(Pag_h); res.Pcu_r_h = Pcu_r_h;
res.T_h = T_h; res.nu = H.nu; res.snu = snu;
res.Pcu_s = Pcu_s; res.Pcu_r = Pcu_r; res.Pfe = Pfe; res.Pfw = Pfw;
res.Padd = Padd;                         % pertes supplémentaires en charge (stray)
res.Pout = Pout; res.eta = eta;
res.Xm = Xm; res.Rfe = Rfe; res.Rr = Rr; res.Xr = Xr;
res.kR = Cg.kR; res.kX = Cg.kX; res.xi = Cg.xi;
res.Bg1 = Rm.Bg1; res.S = Rm.S; res.iterFP = it;
end
