function out = dq_startup(ctx, opt)
%DQ_STARTUP  Démarrage direct transitoire (modèle dq paramétré par le MEC).
%
%   out = mec.dq_startup(ctx, opt) simule le DÉMARRAGE direct sur le réseau
%   (direct-on-line) du moteur asynchrone par un modèle dynamique dq
%   (référentiel stationnaire) dont TOUS les paramètres proviennent du MEC :
%       Rs           : résistance statorique (mec.stator_resistance)
%       Lm           : inductance magnétisante SATURÉE (réseau MEC)
%       Lls, Llr     : inductances de fuite (mec.leakage, mec.cage)
%       Rr'          : résistance rotorique reportée (mec.cage)
%   couplé à l'équation mécanique :  J dOmega/dt = Tem - TL - B*Omega.
%
%   État x = [psi_ds, psi_qs, psi_dr, psi_qr, Omega_m] (flux + vitesse).
%   Alimentation triphasée équilibrée appliquée à t=0 (démarrage brusque).
%
%   opt (facultatif) :
%       .tend  : durée de simulation [s]        (défaut 0.5)
%       .TL    : couple de charge [N.m]          (défaut M.mech.TL0)
%       .im0   : courant magnétisant de calage de Lm [A] (défaut 8.3)
%
%   Sorties (out) : t, n_rpm, Tem, ia (courant phase A), im (magnétisant),
%   plus les paramètres dq utilisés (Ls,Lr,Lm,Rr,Rs) pour traçabilité.
%
%   LIMITE : Lm constante (valeur saturée au point de fonctionnement) et
%   Rr' constante (basse fréquence) — l'effet de barre profonde pendant le
%   démarrage (kR variable) et la saturation dynamique du flux sont des
%   raffinements documentés. Le run-up et l'enveloppe couple/courant sont
%   bien capturés.
%
%   Voir aussi : mec.equivalent_circuit, mec.saturation_curve.

M = ctx.M;
if nargin < 2, opt = struct(); end
tend = getdef(opt,'tend',0.5);
TL   = getdef(opt,'TL',   M.mech.TL0);
im0  = getdef(opt,'im0',  8.3);
% TL : scalaire (constant) OU table [t, TL] (profil de charge, ex. releve
% ANSYS 'Plot 1.tab' — montee de type ventilateur puis palier) -> interp.
if isscalar(TL)
    TLf = @(t) TL;
else
    TLf = @(t) interp1(TL(:,1), TL(:,2), min(max(t,TL(1,1)),TL(end,1)), 'linear');
end
% Couple d'ONDULATION optionnel opt.Trip = @(theta_m) [N.m] (carte MST
% mono- ou multi-tranche) superpose a l'equation mecanique. Il est mis a
% l'echelle (psi_s/psi_N)^2 : l'ondulation de denture est ~ proportionnelle
% au carre du flux, ce qui evite une ondulation non physique tant que le
% flux n'est pas etabli au demarrage.
Trip = getdef(opt,'Trip',[]);
hasT = ~isempty(Trip);

p = M.p; w = M.w; Vm = sqrt(2)*M.Uph;
Rs = ctx.Rs;

% Paramètres dq issus du MEC
Rm = mec.magnetizing(ctx, im0);   Lm = Rm.Lm;      % magnétisante saturée
Lls = ctx.Lk.Xs_leak/w;                            % fuite stator
Cg  = mec.cage(M, ctx.G, ctx.W, ctx.Lk, 0.03);     % cage à faible glissement
Llr = Cg.Xr/w;  Rr = Cg.Rr;
Ls = Lls + Lm;  Lr = Llr + Lm;  D = Ls*Lr - Lm^2;

J = M.mech.J;  B = M.mech.B;

% Second membre du système d'état
psN = Vm/w;                               % flux statorique nominal (crete)
    function dx = odef(t, x)
        pds=x(1); pqs=x(2); pdr=x(3); pqr=x(4); wm=x(5);
        ids=( Lr*pds - Lm*pdr)/D;
        iqs=( Lr*pqs - Lm*pqr)/D;
        idr=( Ls*pdr - Lm*pds)/D;
        iqr=( Ls*pqr - Lm*pqs)/D;
        Vds=Vm*cos(w*t);  Vqs=Vm*sin(w*t);
        wr=p*wm;
        Tem=1.5*p*(pds*iqs - pqs*ids);
        Tr=0;
        if hasT
            Tr=min((pds^2+pqs^2)/psN^2,1)*Trip(x(6));
        end
        dx=zeros(6,1);
        dx(1)=Vds - Rs*ids;
        dx(2)=Vqs - Rs*iqs;
        dx(3)=-Rr*idr - wr*pqr;
        dx(4)=-Rr*iqr + wr*pdr;
        dx(5)=(Tem + Tr - TLf(t) - B*wm)/J;
        dx(6)=wm;
    end

x0 = zeros(6,1);
mstep = 2e-4; if hasT, mstep = 5e-5; end   % resoudre la periode de denture
oset = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',mstep);
[t,X] = ode45(@odef, [0 tend], x0, oset);

% Reconstruction des sorties
pds=X(:,1); pqs=X(:,2); pdr=X(:,3); pqr=X(:,4); wm=X(:,5);
ids=( Lr*pds - Lm*pdr)/D;   iqs=( Lr*pqs - Lm*pqr)/D;
idr=( Ls*pdr - Lm*pds)/D;   iqr=( Ls*pqr - Lm*pqs)/D;
Tem=1.5*p*(pds.*iqs - pqs.*ids);
imd=ids+idr; imq=iqs+iqr;
Trv=zeros(size(t));
if hasT
    Trv=min((pds.^2+pqs.^2)/psN^2,1).*arrayfun(Trip,X(:,6));
end

out.t     = t;
out.n_rpm = wm*60/(2*pi);
out.Tem   = Tem + Trv;                % couple TOTAL (fondamental + ondulation)
out.Tem1  = Tem;                      % part fondamentale (dq)
out.Trip  = Trv;                      % part d'ondulation (carte MST)
out.theta = X(:,6);
out.ia    = ids;                      % courant de phase A (Clarke amplitude)
out.im    = sqrt(imd.^2+imq.^2)/sqrt(2);   % magnétisant efficace
out.wm    = wm;
out.params = struct('Rs',Rs,'Rr',Rr,'Lm',Lm,'Lls',Lls,'Llr',Llr,'J',J,'B',B);
end

function v = getdef(s,f,d)
if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end
end
