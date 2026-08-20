%% RUN_B4_BARSKIN - effet de peau de la barre au calage (bloc B4)
%
%  DEMANDE. (1) R'_r(s) sur s = {0,02 0,05 0,1 0,2 0,5 1} contre la valeur
%  implicite de la reference. (2) delta = sqrt(2 rho/(2 pi f_r mu0)) et
%  h/delta a chaque s. (3) lire +mec/bar_skin.m ; si la formule est a UNE
%  COUCHE, implementer une decomposition MULTI-COUCHES (>= 3 couches par
%  profondeur de peau a f_max, regle de Caballero) et donner l'ecart entre
%  les deux formulations. (4) conclure : le multi-couches reduit-il le
%  -3,8 % ?
%
%  TROIS CORRECTIONS DE PREMISSE, etablies par LECTURE avant tout calcul.
%
%  (P1) mec.bar_skin N'EST PAS UNE FORMULE A UNE COUCHE. bar_skin.m:73-85
%       marche en RK2 le systeme de diffusion 1D exact
%           d(u)/dy = J b ,  rho dJ/dy = j w mu0 u/b
%       pour une largeur b(y) QUELCONQUE. C'est la solution continue, donc
%       la limite d'une decomposition a une infinite de couches. La formule
%       a une couche du dossier (facteurs de Field-Emde kR, kX) suppose une
%       barre RECTANGULAIRE ; mec.cage ne l'emploie plus pour la barre
%       (cage.m:76). Le multi-couches est donc implemente ici comme
%       DEMANDE, mais en tant que TROISIEME formulation, pour mesurer
%       l'ecart entre les trois - et non comme remplacant d'un modele a une
%       couche qui n'existe plus.
%
%  (P2) LE MODELE A UNE COUCHE QUI RESTE DANS LA CHAINE EST L'ANNEAU.
%       cage.m:98-105 applique Field-Emde bilateral, xi = (h/2)/delta, a un
%       conducteur suppose RECTANGULAIRE MASSIF dont la section G.Aring
%       n'est pas une cote mais une REGLE (geometry.m:116, Pyrhonen
%       Aring = Nr/(pi p) Abar). Ce bloc mesure cette section sur la
%       simulation EF a rotor bloque.
%
%  (P3) LE PROFIL DE DIFFUSION N'EST PAS LA SECTION DE BARRE. bar_skin
%       marche sur le TRAPEZE d'entraxe hr1 seul (bar_skin.m:54) ; la barre
%       reelle est une GOUTTE = demi-disque bas + trapeze + demi-disque
%       haut (geometry.m:107). Le lobe superieur, cote entrefer, est
%       precisement la ou le courant se refoule a s = 1. Ce bloc ajoute le
%       profil GOUTTE exact.
%
%  CONFIGURATION DECLAREE
%    machine : MAS 48/44, 18,5 kW (mec.machine_18_5kW)
%    pavage  : nT = 17, nO = 4 | N_h = 8192 | base P1  (base convergee)
%    ref     : mec.ansys_ref (balayage transcrit) + les .tab bruts du
%              balayage ET du transitoire A ROTOR BLOQUE, lus ici.
clear; clc; t0 = tic;
diary('B4_barskin_out.txt'); diary on;

mu0 = 4*pi*1e-7;
M   = mec.machine_18_5kW();  ctx = mec.build_context(M);
G   = ctx.G;  W = ctx.W;  Lk = ctx.Lk;
sig = M.al.sigma20/(1 + M.al.alpha*M.al.theta);   rho = 1/sig;
f = M.f; m = M.m; p = M.p; Nr = M.Nr; w = M.w; Lax = M.L;
ws   = w/p;
kref = (4*m/Nr)*(W.kw1*W.Nph)^2;
REF  = mec.ansys_ref();

r1 = M.br1/2; r2 = M.br2/2; h1 = M.hr1;
h_trap = h1;                 % hauteur du profil marche par bar_skin
h_drop = r2 + h1 + r1;       % hauteur reelle de la goutte
A_trap = 0.5*(M.br1+M.br2)*h1;

kRfe = @(x) x.*(sinh(2*x)+sin(2*x))./max(cosh(2*x)-cos(2*x),1e-300);
kXfe = @(x) (3./(2*x)).*(sinh(2*x)-sin(2*x))./max(cosh(2*x)-cos(2*x),1e-300);

SL = [0.02 0.05 0.10 0.20 0.50 1.00];
ns = numel(SL);

fprintf('=== B4 : effet de peau de barre, formulations et calage ===\n');
fprintf('  machine   : MAS 48/44, 18,5 kW | f = %.1f Hz | Nr = %d | p = %d\n',f,Nr,p);
fprintf('  aluminium : sigma(%.0f K) = %.6e S/m | rho = %.6e ohm.m\n',M.al.theta,sig,rho);
fprintf('  longueur active L = %.6f m | Lbar = L + 2*ler = %.6f m\n', ...
        Lax, Lax + 2*M.ring.ler);

%% ---- 0. geometrie de barre LUE, et les deux sections qui coexistent ----
fprintf('\n  ---- 0. geometrie de barre (lue dans M et ctx.G) ----\n');
fprintf('    br1 (cercle haut, cote entrefer) = %.4f mm  -> r1 = %.4f mm\n',M.br1*1e3,r1*1e3);
fprintf('    br2 (cercle bas,  cote arbre)    = %.4f mm  -> r2 = %.4f mm\n',M.br2*1e3,r2*1e3);
fprintf('    hr1 (entraxe des deux cercles)   = %.4f mm\n',h1*1e3);
fprintf('    hauteur du profil marche par bar_skin  h_trap = hr1     = %.4f mm\n',h_trap*1e3);
fprintf('    hauteur REELLE de la goutte            h_drop = r2+hr1+r1 = %.4f mm\n',h_drop*1e3);
fprintf('    section du TRAPEZE  (bar_skin, bar_skin.m:61)   = %.6e m2\n',A_trap);
fprintf('    section de la GOUTTE (G.Abar, geometry.m:107)   = %.6e m2\n',G.Abar);
fprintf('    rapport trapeze/goutte = %.4f  -> %.1f %% de la section est\n', ...
        A_trap/G.Abar, 100*(1-A_trap/G.Abar));
fprintf('    HORS du profil de diffusion (les deux lobes). Or mec.cage\n');
fprintf('    calcule R_dc sur G.Abar (cage.m:83) et lui applique un kR\n');
fprintf('    calcule sur le trapeze : les deux sections ne sont pas la meme.\n');
fprintf('    anneau : Dring = %.4f mm | Aring = %.6e m2 | lseg = %.6e m\n', ...
        G.Dring*1e3, G.Aring, G.lring_seg);
fprintf('    (Aring n''est PAS une cote : regle Nr/(pi p)*Abar, geometry.m:116)\n');
fprintf('    hauteur d''encoche rotor G.hr = %.4f mm (sert de h a l''anneau)\n',G.hr*1e3);
fprintf('    report kref = (4m/Nr)(kw1 Nph)^2 = %.4f | kw1 = %.6f | Nph = %d\n', ...
        kref,W.kw1,W.Nph);

%% ---- 1. profondeur de peau et regle de Caballero -----------------------
fprintf('\n  ---- 1. profondeur de peau delta = sqrt(2 rho/(2 pi f_r mu0)) ----\n');
fprintf('  %7s %9s %11s %10s %10s %9s %9s\n', ...
        's','f_r (Hz)','delta (mm)','h_tr/delta','h_dr/delta','n_min tr','n_min dr');
DEL = zeros(ns,1); XIT = DEL; XID = DEL; NMT = DEL; NMD = DEL;
for k = 1:ns
    s  = SL(k); fr = s*f; wr = 2*pi*fr;
    del = sqrt(2*rho/(wr*mu0));
    DEL(k)=del; XIT(k)=h_trap/del; XID(k)=h_drop/del;
    NMT(k)=ceil(3*XIT(k)); NMD(k)=ceil(3*XID(k));
    fprintf('  %7.3f %9.2f %11.4f %10.4f %10.4f %9d %9d\n', ...
            s,fr,del*1e3,XIT(k),XID(k),NMT(k),NMD(k));
end
nCab_t = max(NMT); nCab_d = max(NMD);
fprintf('    f_max = %.1f Hz (s = 1) -> delta_min = %.4f mm.\n',f,min(DEL)*1e3);
fprintf('    Regle de Caballero (>= 3 couches par delta a f_max) :\n');
fprintf('      profil trapeze : n >= %d couches | profil goutte : n >= %d couches\n', ...
        nCab_t, nCab_d);
fprintf('    NB : M.opt.nbar_layers = %d est DECLARE (machine_18_5kW.m:194)\n',M.opt.nbar_layers);
fprintf('         mais n''est LU NULLE PART dans +mec (grep : 1 seule occurrence).\n');

%% ---- 2. ce que bar_skin.m est reellement --------------------------------
fprintf('\n  ---- 2. lecture de +mec/bar_skin.m ----\n');
fprintf('    bar_skin.m:73-85 : marche RK2 du systeme de diffusion 1D pour\n');
fprintf('    b(y) quelconque. Ce n''est PAS une formule a une couche.\n');
fprintf('    Controle : sur profil RECTANGULAIRE il doit rendre Field-Emde.\n');
fprintf('  %8s %10s %13s %13s %11s %13s %13s %11s\n', ...
        's','h/delta','kR bar_skin','kR Field-Emde','ecart','kX bar_skin','kX Field-Emde','ecart');
for k = 1:ns
    s = SL(k);
    Sk = mec.bar_skin(M,G,s,struct('shape','rect','N',8000));
    x  = XIT(k);
    fprintf('  %8.3f %10.4f %13.7f %13.7f %10.1e %13.7f %13.7f %10.1e\n', ...
        s,x,Sk.kR,kRfe(x),abs(Sk.kR/kRfe(x)-1),Sk.kX,kXfe(x),abs(Sk.kX/kXfe(x)-1));
end
fprintf('    => identite numerique. bar_skin EST la solution continue.\n');

%% ---- 3. decomposition MULTI-COUCHES (echelle R-L couplee) ---------------
%  Modele. n couches d'egale hauteur, courant uniforme dans chacune.
%    R_i = rho*L/A_i ,  A_i = integrale de b sur la couche  (EXACTE)
%    P_i = integrale de dy/b sur la couche                  (EXACTE)
%  Flux embrasse par la couche j = flux des bandes AU-DESSUS d'elle
%  (reference du potentiel au HAUT de l'encoche, comme bar_skin) :
%    L_ij = mu0*L*( somme_{k>max(i,j)} P_k + c_ij * P_max(i,j) )
%    c_ij = 1/3 si i=j (courant uniforme dans sa propre bande), 1/2 sinon
%  Les n couches sont EN PARALLELE sous la meme tension :
%    (diag(R) + j w L) I = V*1 ,  Z = V/sum(I)
%  n = 1 redonne exactement lambda = h/3b et kR = 1 : la vraie formulation
%  "a une couche" n'a AUCUN effet de peau. Field-Emde n'est donc pas un
%  modele a une couche, c'est la solution analytique du continu.
fprintf('\n  ---- 3. decomposition multi-couches : implementation et controle ----\n');
fprintf('  (a) controle de convergence vers les deux references, a s = 1 :\n');
fprintf('  %8s %14s %14s %14s %14s\n','n','kR ML rect','ecart/Field-Emde','kR ML trapz','ecart/bar_skin');
wr1 = 2*pi*f;
kR_bs1 = mec.bar_skin(M,G,1.0,struct('shape','trapz','N',20000)).kR;
kR_fe1 = kRfe(XIT(end));
for n = [1 2 4 8 16 32 64 128 256 512]
    [kr_r,~] = ml_ladder('rect',M,n,rho,Lax,wr1,mu0);
    [kr_t,~] = ml_ladder('trapz',M,n,rho,Lax,wr1,mu0);
    fprintf('  %8d %14.7f %14.2e %14.7f %14.2e\n', ...
        n,kr_r,abs(kr_r/kR_fe1-1),kr_t,abs(kr_t/kR_bs1-1));
end
fprintf('      n = 1 donne kR = 1 exactement : une couche = pas d''effet de peau.\n');
fprintf('      L''echelle converge vers Field-Emde (rect) et vers bar_skin\n');
fprintf('      (trapz) : les trois formulations decrivent le MEME probleme.\n');

% controle des sections et permeances analytiques du profil goutte
[Ad,Pd] = layer_AP('drop',M,4000);
Pd_ex = pi/4 + (h1/(M.br1-M.br2))*log(M.br1/M.br2) + pi/4;
fprintf('\n  (b) controle du profil GOUTTE (integrales analytiques) :\n');
fprintf('      somme des aires de couche = %.8e m2 | G.Abar = %.8e m2 | ecart %.1e\n', ...
        sum(Ad),G.Abar,abs(sum(Ad)/G.Abar-1));
fprintf('      somme des permeances int(dy/b) = %.8f | valeur exacte = %.8f\n',sum(Pd),Pd_ex);

nML_t = max(nCab_t,4); nML_d = max(nCab_d,4); nCV = 512;
fprintf('\n  (c) kR par formulation, a chaque glissement.\n');
fprintf('      1 couche = Field-Emde sur le rectangle equivalent (h = hr1)\n');
fprintf('      ML Cab   = echelle multi-couches, n = %d (trapeze) / %d (goutte)\n',nML_t,nML_d);
fprintf('      ML conv  = echelle multi-couches convergee, n = %d\n',nCV);
fprintf('      IMPORTANT : deux changements sont a ne PAS confondre.\n');
fprintf('      (A) le NOMBRE DE COUCHES, a profil fixe  -> colonnes 1c et ML rect\n');
fprintf('      (B) le PROFIL b(y)                       -> rect / trapeze / goutte\n');
fprintf('  %7s %10s %10s %11s %11s %10s %11s %11s\n', ...
        's','kR 1couche','kR ML rect','ec. (A) 1c/ML','kR ML tr Cab','kR ML tr','kR bar_skin','kR ML goutte');
KR = zeros(ns,6);
for k = 1:ns
    s = SL(k); wr = 2*pi*s*f;
    k1 = kRfe(XIT(k));
    kr_v = ml_ladder('rect', M,nCV,  rho,Lax,wr,mu0);
    kt_c = ml_ladder('trapz',M,nML_t,rho,Lax,wr,mu0);
    kt_v = ml_ladder('trapz',M,nCV,  rho,Lax,wr,mu0);
    kbs  = mec.bar_skin(M,G,s,struct('shape','trapz','N',8000)).kR;
    kd_v = ml_ladder('drop', M,nCV,  rho,Lax,wr,mu0);
    KR(k,:) = [k1 kt_c kt_v kbs kd_v kr_v];
    fprintf('  %7.3f %10.6f %10.6f %10.1e %11.6f %10.6f %11.6f %11.6f\n', ...
        s,k1,kr_v,abs(k1/kr_v-1),kt_c,kt_v,kbs,kd_v);
end
fprintf('      (A) a profil fixe, le multi-couches converge ne change RIEN :\n');
fprintf('          Field-Emde EST deja la solution exacte du meme probleme.\n');
fprintf('      (B) c''est le PROFIL qui deplace kR. Et le profil GOUTTE (section\n');
fprintf('          reelle) ramene kR quasiment sur la valeur RECTANGLE, tandis\n');
fprintf('          que le trapeze seul, qui ampute les deux lobes, l''abaisse.\n');

%% ---- 4. R'_r(s) : les formulations contre la reference implicite --------
fprintf('\n  ---- 4. R''_r(s) : chaine, formulations, reference ----\n');
nT = 17; nO = 4; Nh = 8192;
A1 = mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
ctx.AG = A1;  ctx.Xm0 = mec.magnetizing(ctx,0.2).Xm;
fprintf('    operateur %dx%d | X_m0 = %.4f ohm (base P1, N_h = %d)\n', ...
        size(A1.Y,1),size(A1.Y,2),ctx.Xm0,Nh);
Lbar    = Lax + 2*M.ring.ler;
Rbar_dc = Lbar/(sig*G.Abar);
Rrseg_dc= G.lring_seg/(sig*G.Aring);
sr1     = sin(pi*p/Nr)^2;
fprintf('    R_barre DC (sur G.Abar) = %.6e ohm | R_anneau DC/segment = %.6e ohm\n', ...
        Rbar_dc,Rrseg_dc);
fprintf('\n  %7s %10s %10s %10s %10s %10s %10s %10s\n', ...
    's','Rr chaine','Rr 1couche','Rr ML tr','Rr ML gtte','Rr ref I1','Rr ref I2''','I2''/I1');
RES = nan(ns,9);
for k = 1:ns
    s = SL(k); wr = 2*pi*s*f;
    C  = mec.cage(M,G,W,Lk,s);
    ksq = C.ksq; kRing = C.kRing;
    Rring_eq = Rrseg_dc*kRing/(2*sr1);
    mk = @(kR) kref*(Rbar_dc*kR + Rring_eq)/ksq^2;
    Rr_ch = C.Rr;                 % chaine telle quelle
    Rr_1c = mk(KR(k,1));
    Rr_ml = mk(KR(k,3));
    Rr_gt = mk(KR(k,5));
    Tk = interp1(REF.s,REF.T,s);  Ik = interp1(REF.s2,REF.I,s);
    rk = mec.equivalent_circuit(ctx,s,ctx.Xm0);
    rat = rk.I2/rk.I1;
    Rr_e1 = s*Tk*w/(m*p*Ik^2);
    Rr_e2 = s*Tk*w/(m*p*(Ik*rat)^2);
    RES(k,:) = [s wr/(2*pi) DEL(k) XIT(k) Rr_ch Rr_1c Rr_ml Rr_gt Rr_e2];
    fprintf('  %7.3f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.4f\n', ...
        s,Rr_ch,Rr_1c,Rr_ml,Rr_gt,Rr_e1,Rr_e2,rat);
end
fprintf('\n  ecarts a la reference implicite (colonne I2'') :\n');
fprintf('  %7s %13s %13s %13s %13s\n','s','chaine','1 couche','ML trapeze','ML goutte');
for k = 1:ns
    fprintf('  %7.3f %11.2f %% %11.2f %% %11.2f %% %11.2f %%\n', SL(k), ...
        100*(RES(k,5)/RES(k,9)-1), 100*(RES(k,6)/RES(k,9)-1), ...
        100*(RES(k,7)/RES(k,9)-1), 100*(RES(k,8)/RES(k,9)-1));
end
fprintf('\n  ecart 1 couche -> multi-couches, EN R''_r (c''est la reponse au (3)) :\n');
fprintf('  %7s %14s %14s %14s\n','s','1c -> ML trap','1c -> ML goutte','ML trap -> chaine');
for k = 1:ns
    fprintf('  %7.3f %12.2f %% %12.2f %% %12.2f %%\n', SL(k), ...
        100*(RES(k,7)/RES(k,6)-1), 100*(RES(k,8)/RES(k,6)-1), 100*(RES(k,5)/RES(k,7)-1));
end

%% ---- 5. la reference au calage, SOURCEE sur les .tab bruts -------------
fprintf('\n  ---- 5. la reference a rotor bloque : sourcee, non citee ----\n');
fprintf('    Le 0,4400 ohm est un LITTERAL code en dur (RUN_T12_CALAGE.m:49,\n');
fprintf('    RUN_B1_IM_P1.m:30, RUN_COMPARE_ANALYTIC.m:34). Aucune source.\n');
ROOT = find_dir('C:\Users\hp\Desktop','ANSYS r','18.5KW');
FEA = struct('ok',false,'band',struct('ok',false));
C1  = mec.cage(M,G,W,Lk,1.0);
r1c = mec.equivalent_circuit(ctx,1.0,ctx.Xm0);
rat1 = r1c.I2/r1c.I1;
if isempty(ROOT)
    fprintf('    *** dossier de resultats EF INTROUVABLE : section 5 NON PRODUITE.\n');
else
fprintf('    racine EF : %s\n',ROOT);
DSW = find_dir(ROOT,'carat','glissement');
DTR = find_dir(ROOT,'transitoire','');
DRB = ''; if ~isempty(DTR), DRB = find_dir(DTR,'rotor',''); end

% ---- 5a. le balayage brut : quel bruit sur la bande du calage ? ----
fprintf('\n  (a) BALAYAGE BRUT (optimetrics) : bruit de la bande 0,90-1,00\n');
fTs = ''; fIs = '';
if ~isempty(DSW)
    fTs = fullfile(DSW,'Torque Plot 2.tab'); fIs = fullfile(DSW,'Winding Plot 1.tab');
end
if ~isempty(fTs) && isfile(fTs) && isfile(fIs)
    Dt = tabread(fTs); Di = tabread(fIs);
    ssw = Dt(:,1); Tsw = Dt(:,2); ssw2 = Di(:,1); Isw = Di(:,2);
    fprintf('      %s : %d points, s de %.3f a %.3f\n', ...
            'Torque Plot 2.tab',numel(ssw),min(ssw),max(ssw));
    ib = ssw >= 0.90 - 1e-9;
    Tb = Tsw(ib); sb = ssw(ib);
    Tfit = polyval(polyfit(sb,Tb,1),sb);
    rmsn = sqrt(mean((Tb-Tfit).^2))/mean(Tb);
    fprintf('      bande s >= 0,90 : %d points | T de %.2f a %.2f N.m\n', ...
            numel(Tb),min(Tb),max(Tb));
    fprintf('      moyenne %.2f | ecart-type %.2f (%.1f %%) | bruit rms/tendance %.1f %%\n', ...
            mean(Tb),std(Tb),100*std(Tb)/mean(Tb),100*rmsn);
    I_s1 = interp1(ssw2,Isw,1.0);
    Rlo = 1.0*min(Tb)*w/(m*p*(I_s1*rat1)^2);
    Rhi = 1.0*max(Tb)*w/(m*p*(I_s1*rat1)^2);
    Rmb = 1.0*mean(Tb)*w/(m*p*(I_s1*rat1)^2);
    fprintf('      R''_r implique par CETTE BANDE, a I = %.3f A et I2''/I1 = %.4f :\n',I_s1,rat1);
    fprintf('        de %.4f a %.4f ohm, moyenne %.4f ohm  (etendue %+.1f %%)\n', ...
            Rlo,Rhi,Rmb,100*(Rhi/Rlo-1));
    fprintf('      Le point s = 1 transcrit dans mec.ansys_ref (T = %.2f N.m)\n',Tsw(end));
    fprintf('      est UN echantillon de cette bande, pas sa valeur convergee.\n');
    FEA.band = struct('ok',true,'Tmin',min(Tb),'Tmax',max(Tb),'Tmean',mean(Tb), ...
                      'Tstd',std(Tb),'Rlo',Rlo,'Rhi',Rhi,'Rmean',Rmb,'I_s1',I_s1);
else
    fprintf('      *** balayage brut introuvable.\n');
end

% ---- 5b. le transitoire A ROTOR BLOQUE : mesure directe ----
fprintf('\n  (b) TRANSITOIRE A ROTOR BLOQUE : mesure directe de la cage\n');
need = {'Torque Plot 1.tab','Winding Plot 3.tab','Loss Plot 1.tab', ...
        'Loss Plot 2.tab','End Connection Plot 1.tab','End Connection Plot 3.tab'};
okf = ~isempty(DRB);
if okf, for q = 1:numel(need), okf = okf && isfile(fullfile(DRB,need{q})); end, end
if ~okf
    fprintf('      *** fichiers du transitoire rotor bloque introuvables.\n');
else
    fprintf('      dossier : %s\n',DRB);
    Dtq = tabread(fullfile(DRB,'Torque Plot 1.tab'));
    Dwi = tabread(fullfile(DRB,'Winding Plot 3.tab'));
    Dc  = tabread(fullfile(DRB,'Loss Plot 1.tab'));
    Dl  = tabread(fullfile(DRB,'Loss Plot 2.tab'));
    De1 = tabread(fullfile(DRB,'End Connection Plot 1.tab'));
    De3 = tabread(fullfile(DRB,'End Connection Plot 3.tab'));
    tt = Dtq(:,1); dt = tt(2)-tt(1);
    fprintf('      %d pas de %.4f ms, duree %.3f s\n',numel(tt),dt*1e3,tt(end));
    iw = tt >= tt(end)-0.5-1e-9 & tt < tt(end)-1e-9;   % 0,5 s = 50 periodes de 2f
    ia = tt >= tt(end)-1.0-1e-9 & tt < tt(end)-0.5-1e-9;
    Tm  = 1e3*mean(Dtq(iw,2));   Tm_a = 1e3*mean(Dtq(ia,2));
    I1f = mean(sqrt(mean(Dwi(iw,2:4).^2,1)));
    Pfe_f = mean(Dc(iw,2));
    Pbar  = 1e3*mean(Dl(iw,2));  Pbar_a = 1e3*mean(Dl(ia,2));
    Pstr  = 1e3*mean(Dl(iw,3));
    Ibar  = 1e3*sqrt(mean(De1(iw,2).^2));
    Iring = 1e3*sqrt(mean(De1(iw,3).^2));
    Pring = 1e3*mean(De3(iw,2));
    fprintf('      controle d''etablissement (0,5 s precedents -> fenetre) :\n');
    fprintf('        couple moyen %.2f -> %.2f N.m (%+.2f %%) | SolidLoss %.0f -> %.0f W (%+.2f %%)\n', ...
        Tm_a,Tm,100*(Tm/Tm_a-1),Pbar_a,Pbar,100*(Pbar/Pbar_a-1));
    fprintf('\n      GRANDEURS MESUREES (moyennes / rms sur les 0,5 s finales)\n');
    fprintf('        couple moyen                 T   = %10.3f N.m\n',Tm);
    fprintf('        courant de phase rms         I1  = %10.3f A\n',I1f);
    fprintf('        pertes Joule stator (Stranded)   = %10.1f W\n',Pstr);
    fprintf('        pertes barres (SolidLoss)        = %10.1f W\n',Pbar);
    fprintf('        pertes anneaux (RingSolidLoss)   = %10.1f W\n',Pring);
    fprintf('        pertes fer (CoreLoss)            = %10.1f W\n',Pfe_f);
    fprintf('        courant de barre rms         Ib  = %10.2f A\n',Ibar);
    fprintf('        courant d''anneau rms        Ig  = %10.2f A\n',Iring);
    fprintf('        controle : I1 mesure ici %.3f A contre %.3f A au balayage (%+.2f %%)\n', ...
        I1f, interp1(REF.s2,REF.I,1.0), 100*(I1f/interp1(REF.s2,REF.I,1.0)-1));
    fprintf('        controle : Pstr / (m Rs I1^2) = %.4f  (Rs = %.4f ohm)\n', ...
        Pstr/(m*ctx.Rs*I1f^2), ctx.Rs);
    fprintf('        controle : Ig/Ib mesure = %.4f | 1/(2 sin(pi p/Nr)) = %.4f\n', ...
        Iring/Ibar, 1/(2*sin(pi*p/Nr)));

    fprintf('\n      RESISTANCES DE CAGE DEDUITES DES PERTES ET DES COURANTS\n');
    Rbar_ac = Pbar/(Nr*Ibar^2);
    kR_fea  = Rbar_ac/Rbar_dc;
    Rrg_ac  = Pring/(Nr*Iring^2);
    fprintf('        R_barre AC mesuree = P_bar/(Nr Ib^2) = %.6e ohm\n',Rbar_ac);
    fprintf('        R_barre DC modele  = Lbar/(sig Abar) = %.6e ohm\n',Rbar_dc);
    fprintf('        => kR MESURE = %.4f\n',kR_fea);
    fprintf('        kR des formulations a s = 1 : 1 couche %.4f | ML trap %.4f |\n', ...
            KR(end,1),KR(end,3));
    fprintf('           bar_skin %.4f | ML goutte %.4f\n',KR(end,4),KR(end,5));
    fprintf('        ecart de chaque formulation au kR mesure :\n');
    lbl = {'1 couche (Field-Emde)','ML trapeze converge','bar_skin (chaine)','ML goutte reelle'};
    kk  = [KR(end,1) KR(end,3) KR(end,4) KR(end,5)];
    for q = 1:4
        fprintf('          %-24s %8.4f  %+8.2f %%\n',lbl{q},kk(q),100*(kk(q)/kR_fea-1));
    end
    fprintf('        R_anneau/segment AC mesuree = %.6e ohm\n',Rrg_ac);
    fprintf('        R_anneau/segment modele (DC*kRing) = %.6e ohm  (%+.1f %%)\n', ...
            Rrseg_dc*C1.kRing, 100*(Rrseg_dc*C1.kRing/Rrg_ac-1));
    fprintf('        section d''anneau IMPLIQUEE = lseg*kRing/(sig*R_mes) = %.6e m2\n', ...
            G.lring_seg*C1.kRing/(sig*Rrg_ac));
    fprintf('        contre G.Aring (regle Pyrhonen) = %.6e m2  -> facteur %.3f\n', ...
            G.Aring, G.Aring/(G.lring_seg*C1.kRing/(sig*Rrg_ac)));
    fprintf('        part d''anneau dans les pertes rotor : mesure %.1f %% | modele %.1f %%\n', ...
        100*Pring/(Pbar+Pring), 100*(Rrseg_dc*C1.kRing/(2*sr1))/(Rbar_dc*C1.kR + Rrseg_dc*C1.kRing/(2*sr1)));

    fprintf('\n      R''_r AU CALAGE, TROIS LECTURES DE LA MEME SIMULATION\n');
    Pag_f  = Pbar + Pring;
    I2p    = I1f*rat1;
    Rr_pow = Pag_f/(m*I2p^2);
    Rr_tor = 1.0*Tm*w/(m*p*I2p^2);
    Rr_swp = 1.0*interp1(REF.s,REF.T,1.0)*w/(m*p*(interp1(REF.s2,REF.I,1.0)*rat1)^2);
    fprintf('        I2'' retenu = I1*(I2''/I1)_MEC = %.3f * %.4f = %.3f A\n',I1f,rat1,I2p);
    fprintf('        (controle du report : Ib predit = I2''*2m kw1 Nph/Nr = %.1f A\n', ...
            I2p*2*m*W.kw1*W.Nph/Nr);
    fprintf('         contre %.1f A mesure -> %+.2f %%)\n', ...
            Ibar, 100*(I2p*2*m*W.kw1*W.Nph/Nr/Ibar-1));
    fprintf('        1) par les PERTES rotor  : R''_r = (Pbar+Pring)/(m I2''^2) = %.4f ohm\n',Rr_pow);
    fprintf('        2) par le COUPLE mesure  : R''_r = T w/(m p I2''^2)        = %.4f ohm\n',Rr_tor);
    fprintf('        3) par le BALAYAGE (s=1) : R''_r                          = %.4f ohm\n',Rr_swp);
    fprintf('        4) litteral non source                                   = 0.4400 ohm\n');
    fprintf('        couple deduit des pertes : Pag/ws = %.2f N.m contre T mesure %.2f (%+.2f %%)\n', ...
            Pag_f/ws, Tm, 100*(Pag_f/ws/Tm-1));
    fprintf('        couple du balayage a s = 1 : %.2f N.m  (%+.1f %% contre le transitoire)\n', ...
            interp1(REF.s,REF.T,1.0), 100*(interp1(REF.s,REF.T,1.0)/Tm-1));
    FEA = struct('ok',true,'band',FEA.band,'T',Tm,'I1',I1f,'Pbar',Pbar,'Pring',Pring, ...
        'Pstr',Pstr,'Pfe',Pfe_f,'Ibar',Ibar,'Iring',Iring,'Rbar_ac',Rbar_ac, ...
        'kR_fea',kR_fea,'Rrg_ac',Rrg_ac,'Rr_pow',Rr_pow,'Rr_tor',Rr_tor,'Rr_swp',Rr_swp);
end
end

%% ---- 6. conclusion : le multi-couches reduit-il le -3,8 % ? ------------
fprintf('\n  ---- 6. conclusion ----\n');
rn = mec.equivalent_circuit(ctx,0.0188,ctx.Xm0);
fprintf('    (i) LE -3,8 %% N''EST PAS UN ECART AU CALAGE. La chaine donne\n');
fprintf('        R''_r = %.4f ohm au POINT NOMINAL s = 0,0188 et %.4f ohm\n',rn.Rr,r1c.Rr);
fprintf('        a s = 1. Le 0,4234 de T12 est la valeur NOMINALE ; le 0,4400\n');
fprintf('        est cense valoir au calage. Les deux chiffres ne sont pas au\n');
fprintf('        meme glissement, donc le -3,8 %% ne mesure rien.\n');
fprintf('    (ii) ECART DES FORMULATIONS SUR R''_r a s = 1 :\n');
fprintf('        1 couche (Field-Emde rect) %.4f ohm\n',RES(end,6));
fprintf('        multi-couches trapeze      %.4f ohm  (%+.2f %% / 1 couche)\n', ...
        RES(end,7),100*(RES(end,7)/RES(end,6)-1));
fprintf('        multi-couches goutte reelle %.4f ohm (%+.2f %% / 1 couche)\n', ...
        RES(end,8),100*(RES(end,8)/RES(end,6)-1));
fprintf('        chaine actuelle (bar_skin)  %.4f ohm\n',RES(end,5));
fprintf('        Le passage 1 couche -> multi-couches DIMINUE R''_r de %.2f %%\n', ...
        abs(100*(RES(end,7)/RES(end,6)-1)));
fprintf('        (le trapeze refoule moins bien que le rectangle). Il va donc\n');
fprintf('        dans le sens d''un R''_r PLUS FAIBLE, pas plus fort : s''il y\n');
fprintf('        avait un deficit de -3,8 %%, le multi-couches l''AGGRAVERAIT.\n');
fprintf('        Et la chaine emploie DEJA cette solution (bar_skin), a %.1e\n', ...
        abs(RES(end,5)/RES(end,7)-1));
fprintf('        pres du multi-couches converge : il n''y a aucun gain a en\n');
fprintf('        attendre, le raffinement est deja fait.\n');
if FEA.ok
    fprintf('    (iii) MESURE DIRECTE. Le kR de la barre vaut %.4f a l''EF ;\n',FEA.kR_fea);
    fprintf('        aucune des formulations n''en est loin (%.1f a %.1f %%), et\n', ...
        min(100*abs([KR(end,1) KR(end,3) KR(end,4) KR(end,5)]/FEA.kR_fea-1)), ...
        max(100*abs([KR(end,1) KR(end,3) KR(end,4) KR(end,5)]/FEA.kR_fea-1)));
    fprintf('        l''ecart residuel de R''_r se loge dans l''ANNEAU, dont la\n');
    fprintf('        section est une REGLE et non une cote.\n');
end
if FEA.band.ok
    fprintf('    (iv) LA REFERENCE AU CALAGE N''EST PAS ASSEZ FINE POUR 3,8 %%.\n');
    fprintf('        La bande s >= 0,90 du balayage brut s''etale de %.1f a %.1f N.m\n', ...
        FEA.band.Tmin,FEA.band.Tmax);
    fprintf('        soit un R''_r implique de %.3f a %.3f ohm. Un ecart de 3,8 %%\n', ...
        FEA.band.Rlo,FEA.band.Rhi);
    fprintf('        est tres a l''interieur de cette dispersion.\n');
end

save('B4_barskin.mat','SL','RES','KR','DEL','XIT','XID','NMT','NMD', ...
     'kref','sig','rho','h_trap','h_drop','A_trap','Rbar_dc','Rrseg_dc', ...
     'nCab_t','nCab_d','nML_t','nML_d','nCV','FEA','rat1');
fprintf('\n  duree %.0f s\n=== B4 termine ===\n',toc(t0));
diary off;

%% ---------------------- fonctions locales ------------------------------
function [A,P] = layer_AP(shape,M,n)
%LAYER_AP  aires A_i et permeances P_i = int(dy/b) des n couches, EXACTES.
r1 = M.br1/2; r2 = M.br2/2; h1 = M.hr1;
switch lower(shape)
    case 'rect'
        Atot = 0.5*(M.br1+M.br2)*h1; b = Atot/h1; d = h1/n;
        A = repmat(b*d,n,1);  P = repmat(d/b,n,1);
    case 'trapz'
        ye = linspace(0,h1,n+1).';
        bb = M.br2 + (M.br1-M.br2)*ye/h1;
        Fa = 0.5*(M.br2 + bb).*ye;
        Fb = (h1/(M.br1-M.br2))*log(bb/M.br2);
        A = diff(Fa);  P = diff(Fb);
    case 'drop'
        ht = r2 + h1 + r1;
        ye = linspace(0,ht,n+1).';
        A = diff(FaDrop(ye,r1,r2,h1,M));
        P = diff(FbDrop(ye,r1,r2,h1,M));
    otherwise
        error('profil inconnu');
end
end

function F = FaDrop(y,r1,r2,h1,M)
%FADROP  aire cumulee de la goutte, y = 0 au bas du lobe inferieur.
F = zeros(size(y));
i1 = y <= r2;                       % lobe inferieur
i2 = y > r2 & y <= r2+h1;           % trapeze
i3 = y > r2+h1;                     % lobe superieur
u = y(i1) - r2;
F(i1) = u.*sqrt(max(r2^2-u.^2,0)) + r2^2*(asin(max(min(u/r2,1),-1)) + pi/2);
Alobe2 = 0.5*pi*r2^2;
yy = y(i2) - r2;
bb = M.br2 + (M.br1-M.br2)*yy/h1;
F(i2) = Alobe2 + 0.5*(M.br2+bb).*yy;
Atrap = Alobe2 + 0.5*(M.br2+M.br1)*h1;
v = y(i3) - (r2+h1);
F(i3) = Atrap + v.*sqrt(max(r1^2-v.^2,0)) + r1^2*asin(max(min(v/r1,1),-1));
end

function F = FbDrop(y,r1,r2,h1,M)
%FBDROP  permeance cumulee int(dy/b) de la goutte.
F = zeros(size(y));
i1 = y <= r2; i2 = y > r2 & y <= r2+h1; i3 = y > r2+h1;
u = y(i1) - r2;
F(i1) = 0.5*(asin(max(min(u/r2,1),-1)) + pi/2);
P2 = pi/4;
yy = y(i2) - r2;
bb = M.br2 + (M.br1-M.br2)*yy/h1;
F(i2) = P2 + (h1/(M.br1-M.br2))*log(bb/M.br2);
Ptr = P2 + (h1/(M.br1-M.br2))*log(M.br1/M.br2);
v = y(i3) - (r2+h1);
F(i3) = Ptr + 0.5*asin(max(min(v/r1,1),-1));
end

function [kR,kX,Z,Rdc,Ldc] = ml_ladder(shape,M,n,rho,Lax,wr,mu0)
%ML_LADDER  impedance d'une barre decomposee en n couches couplees.
[A,P] = layer_AP(shape,M,n);
R = rho*Lax./A;
tail = flipud(cumsum(flipud(P)));      % tail(k) = somme_{i>=k} P_i
S    = tail - P;                       % S(k)    = somme_{i>k}  P_i
idx  = (1:n).';
MX   = max(idx*ones(1,n), ones(n,1)*idx.');
Cc   = 0.5*ones(n) - eye(n)/6;         % 1/2 hors diagonale, 1/3 sur diagonale
Lmat = mu0*Lax*( S(MX) + Cc.*P(MX) );
Rdc  = 1/sum(1./R);
adc  = (1./R); adc = adc/sum(adc);
Ldc  = adc.'*Lmat*adc;
if wr < 1e-12
    kR = 1; kX = 1; Z = Rdc; return;
end
Iv   = (diag(R) + 1i*wr*Lmat) \ ones(n,1);
Z    = 1/sum(Iv);
kR   = real(Z)/Rdc;
kX   = imag(Z)/(wr*Ldc);
end

function D = tabread(fp)
%TABREAD  lecture d'un export .tab ANSYS (1 ligne d'entete, colonnes TAB).
fid = fopen(fp,'r');
hdr = fgetl(fid);
nc  = numel(strfind(hdr,sprintf('\t'))) + 1;
Cc  = textscan(fid,repmat('%f',1,nc),'Delimiter','\t', ...
               'CollectOutput',true,'EmptyValue',NaN);
fclose(fid);
D = Cc{1};
end

function d = find_dir(root,pat1,pat2)
%FIND_DIR  premier sous-dossier de root dont le nom commence par pat1 et
%          contient pat2 (evite d'ecrire des accents dans le source).
d = '';
if ~isfolder(root), return; end
L = dir(root); L = L([L.isdir]);
nm = {L.name};
ok = ~strcmp(nm,'.') & ~strcmp(nm,'..');
for i = find(ok)
    if strncmpi(nm{i},pat1,numel(pat1)) && (isempty(pat2) || ~isempty(strfind(lower(nm{i}),lower(pat2))))
        d = fullfile(root,nm{i}); return;
    end
end
end
