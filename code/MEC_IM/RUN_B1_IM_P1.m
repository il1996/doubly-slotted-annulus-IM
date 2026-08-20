%% RUN_B1_IM_P1 - chaine MAS complete en base P1
%
%  B1. TOUTES les grandeurs derivees de X_m changent avec la base. On
%  regenere donc en UNE execution, base P1, pavage et N_h converges (B2) :
%    - schema equivalent : Rs, Rr', Xsigma_s, Xsigma_r, Xm, Rfe
%    - point nominal s = 0,0188 : couple, I1, rendement, cos phi, pertes fer
%    - les trois essais : en charge, a vide, rotor bloque
%    - B_g1 et B_t a mi-entrefer, a vide et en charge
%    - T(s) et I(s) sur 30 glissements, decrochage
%
%  CONFIGURATION DECLAREE : pavage nT=17 nO=4, N_h=8192, base P1.
%  B2 y etablit k_C = 1,3320 stable a 0,62 % et X_m = 61,015 ohm.
clear; clc; t0=tic;
diary('B1_im_p1_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G;
nT=17; nO=4; Nh=8192; s_ch=0.0188;
A1=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
ctx.AG=A1; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;

fprintf('=== B1 : chaine MAS complete en base P1 ===\n');
fprintf('  pavage nT=%d nO=%d | N_h = %d | base P1\n',nT,nO,Nh);
fprintf('  operateur %dx%d | X_m0 = %.3f ohm\n\n',size(A1.Y,1),size(A1.Y,2),ctx.Xm0);

%% ---- 1. schema equivalent -------------------------------------------
r=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
fprintf('  ---- 1. SCHEMA EQUIVALENT (point nominal s = %.4f) ----\n',s_ch);
fprintf('  %-22s %12s %12s\n','parametre','P1','EF');
pr=@(l,a,f)fprintf('  %-22s %12.4f %12.4f\n',l,a,f);
pr('Rs [ohm]',ctx.Rs,0.4450);
pr('Rr'' [ohm]',r.Rr,0.4400);
%  Les reactances de fuite ne portent pas le meme nom selon la source ; on
%  les cherche dans l'ordre de priorite et on DECLARE celle qui est trouvee.
gf=@(s,c) deal_first(s,c);
[Xss,nss]=gf(struct('r',r,'Lk',ctx.Lk),{'r.Xs','r.Xsig_s','Lk.Xs','Lk.Xsig','Lk.Xs_slot'});
[Xrr,nrr]=gf(struct('r',r,'Lk',ctx.Lk),{'r.Xr','r.Xsig_r','Lk.Xr','Lk.Xr_slot'});
fprintf('  %-22s %12.4f %12s   (source : %s)\n','Xsigma_s [ohm]',Xss,'--',nss);
fprintf('  %-22s %12.4f %12s   (source : %s)\n','Xsigma_r [ohm]',Xrr,'--',nrr);
fprintf('  (les EF ne separent pas les deux fuites : pas de colonne EF)\n');
pr('Xm sature [ohm]',r.Xm,46.0);
pr('Rfe [ohm]',r.Rfe,1740);
fprintf('  X_m NON sature : %.3f ohm\n',ctx.Xm0);

%% ---- 2. point nominal ------------------------------------------------
fprintf('\n  ---- 2. POINT NOMINAL ----\n');
fprintf('  %-22s %12s %12s %10s\n','grandeur','P1','EF','ecart');
q=@(l,a,f)fprintf('  %-22s %12.4f %12.4f %9.1f %%\n',l,a,f,100*(a-f)/f);
q('couple [N.m]',r.Tem,121.53);
q('courant I1 [A]',r.I1,19.70);
q('rendement',r.eta,0.9160);
q('cos(phi)',r.cosphi,0.8640);
q('pertes fer [W]',r.Pfe,232.6);

%% ---- 3. les trois essais ---------------------------------------------
fprintf('\n  ---- 3. TROIS CONDITIONS ----\n');
r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
rb=mec.equivalent_circuit(ctx,1.0,ctx.Xm0);
fprintf('  EN CHARGE (s = %.4f)\n',s_ch);
q('  couple [N.m]',r.Tem,121.63);
q('  I1 [A]',r.I1,19.73);
q('  I barre [A]',r.I2*(2*M.m*ctx.W.kw1*ctx.W.Nph)/M.Nr,324.7);
q('  pertes fer [W]',r.Pfe,232.6);
fprintf('  A VIDE\n');
q('  I magnetisant [A]',r0.I1,8.49);
q('  f.e.m. entrefer [V]',r0.E1,382.1);
q('  pertes fer [W]',r0.Pfe,249.3);
fprintf('  ROTOR BLOQUE (s = 1)\n');
q('  couple [N.m]',rb.Tem,104.31);
q('  I1 [A]',rb.I1,108.21);
q('  I barre [A]',rb.I2*(2*M.m*ctx.W.kw1*ctx.W.Nph)/M.Nr,1929);

%% ---- 4. champs a mi-entrefer -----------------------------------------
fprintf('\n  ---- 4. CHAMPS A MI-ENTREFER ----\n');
Rm=0.5*(G.Rs+G.Rr); thq=linspace(0,2*pi,2001); thq(end)=[];
%  a vide : branche magnetisante seule
Rnl=mec.magnetizing(ctx,r0.Im);
[Br0,Bt0]=A1.field(Rnl.S.Usurf,Rm,thq); Br0=Br0(:).'; Bt0=Bt0(:).';
%  en charge : courant statorique TOTAL + onde de FMM de barres
p1=angle(r.I1c);
i3=sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
Fs=ctx.W.slotMMF(i3);
p2=angle(r.I2c); tb=2*pi*(0:M.Nr-1)/M.Nr;
Fu=cumsum(cos(M.p*tb-p2).'); Fu=Fu-mean(Fu);
c1=(2/M.Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
Fr=-Fu*((3/2)*(4/pi)*(ctx.W.kw1*ctx.W.Nph/(2*M.p))*sqrt(2)*abs(r.I2c)/abs(c1));
SL=mec.solve_network(ctx.net,G,ctx.BH,A1,Fs,Fr,M.opt);
[Br1,Bt1]=A1.field(SL.Usurf,Rm,thq); Br1=Br1(:).'; Bt1=Bt1(:).';
ff=@(y)abs((2/numel(y))*sum(y.*exp(-1j*M.p*thq)));
q('Bg1 a vide [T]',ff(Br0),0.942);
q('Bt rms a vide [T]',sqrt(mean(Bt0.^2)),0.128);
q('Bg1 en charge [T]',ff(Br1),0.920);
q('Bt rms en charge [T]',sqrt(mean(Bt1.^2)),0.131);

%% ---- 5. caracteristique sur 30 glissements ---------------------------
fprintf('\n  ---- 5. CARACTERISTIQUE (30 glissements) ----\n');
sl=[0.005:0.005:0.12,0.15,0.2,0.3,0.5,0.7,1.0];
T=nan(size(sl)); I=nan(size(sl)); Xp=ctx.Xm0; b6=0;
for k=1:numel(sl)
    rk=mec.equivalent_circuit(ctx,sl(k),Xp); Xp=rk.Xm;
    T(k)=rk.Tem; I(k)=rk.I1;
    c=mec.power_balance(rk,M); b6=max(b6,c.err_global);
end
[Tmx,im]=max(T);
fprintf('  bilan de puissance B6 : %.2e sur %d glissements\n',b6,numel(sl));
fprintf('  decrochage : %.1f N.m a s = %.3f  (EF 324.9 a s = 0.105)\n',Tmx,sl(im));
fprintf('  %6s %10s %10s\n','s','T [N.m]','I1 [A]');
for k=[1 4 10 20 24 numel(sl)]
    fprintf('  %6.3f %10.2f %10.2f\n',sl(k),T(k),I(k));
end
save('B1_im_p1.mat','sl','T','I');
fprintf('\n  duree %.0f s\n=== B1 termine ===\n',toc(t0));
diary off;

% ======================================================================
function [v,nm]=deal_first(S,cands)
%  Renvoie la premiere valeur trouvee parmi les chemins 'champ.sous_champ'
%  et le NOM du chemin retenu, pour que la source soit declaree et non
%  supposee. Renvoie NaN et 'introuvable' si aucun ne repond.
    v=NaN; nm='introuvable';
    for k=1:numel(cands)
        p=strsplit(cands{k},'.');
        if ~isfield(S,p{1}), continue; end
        x=S.(p{1});
        if numel(p)==2
            if ~isstruct(x)||~isfield(x,p{2}), continue; end
            x=x.(p{2});
        end
        if isnumeric(x)&&isscalar(x)&&isfinite(x)
            v=x; nm=cands{k}; return;
        end
    end
end
