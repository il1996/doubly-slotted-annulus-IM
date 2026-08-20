%% RUN_B10_B1_SKEWOFF - B1 complet, vrillage harmonique neutralise
%
%  POURQUOI. B8 a mesure l'exposition au vrillage mais n'a collecte qu'un
%  sous-ensemble des grandeurs de B1. Trois cellules de la section 5 de
%  l'Article II restent donc a la configuration de production alors que le
%  reste de la section est a la configuration LOYALE :
%
%    - R_fe au point nominal        (tab:im_ec, note b)
%    - pertes fer A VIDE            (tab:im_tests, ligne absente)
%    - la caracteristique T(s), I(s) et le decrochage
%
%  Ce script est RUN_B1_IM_P1 a l'identique, avec M.opt.skew_harm = 0.
%  Sa sortie remplace les trois cellules signalees et permet de retirer la
%  note (b) de tab:im_ec.
%
%  RAPPEL DU MOTIF. Le projet EF declare UseSkewModel=true, SkewAngle
%  7,5 deg, mais NumberOfSlices = 1 : une seule tranche ne produit aucune
%  rotation relative, donc aucun moyennage axial, donc une section droite
%  NON vrillee. La chaine MEC, elle, applique ksq_nu aux branches
%  harmoniques par defaut (cage.m:59-61). Neutraliser est la seule facon
%  de comparer deux modeles dans le meme etat.
%
%  RESIDU IRREDUCTIBLE : au fondamental cage.m:61 donne asq = asq1 dans
%  les deux cas, donc ksq = 0,997147 subsiste et R'_r, X'_r restent
%  divises par ksq^2 -- +0,57 % que la reference n'a pas. A DECLARER.
%
%  REFERENCES EF : famille TRANSITOIRE (decision 1 du 6 aout 2026).
%    en charge 121,63 N.m / 19,73 A | a vide 8,49 A | calage 104,31 N.m
%  Le point s = 1 de la caracteristique (98,98 N.m) n'est PAS employe :
%  sa dispersion locale atteint 16 % pour s >= 0,90.
%
%  Sortie : B10_b1_skewoff_out.txt
clear; clc; t0=tic;
diary('B10_b1_skewoff_out.txt'); diary on;

M = mec.machine_18_5kW();
M.opt.skew_harm = 0;                       % <<< la seule difference avec B1
ctx = mec.build_context(M); G = ctx.G;
ctx.M.opt.skew_harm = 0;                   % securite : c'est ctx.M que cage lit
nT=17; nO=4; Nh=8192; s_ch=0.0188;
A1 = mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
ctx.AG = A1; ctx.Xm0 = mec.magnetizing(ctx,0.2).Xm;

fprintf('=== B10 : chaine MAS complete, vrillage harmonique NEUTRALISE ===\n');
fprintf('  pavage nT=%d nO=%d | N_h = %d | base P1 | skew_harm = 0\n',nT,nO,Nh);
fprintf('  operateur %dx%d | X_m0 = %.3f ohm\n',size(A1.Y,1),size(A1.Y,2),ctx.Xm0);
fprintf('  references EF : famille TRANSITOIRE\n\n');

%% ---- 1. schema equivalent -------------------------------------------
r = mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
fprintf('  ---- 1. SCHEMA EQUIVALENT (point nominal s = %.4f) ----\n',s_ch);
fprintf('  %-24s %12s %12s\n','parametre','P1 skew=0','EF');
pr=@(l,a,f)fprintf('  %-24s %12.4f %12.4f\n',l,a,f);
pr('Rs [ohm]',ctx.Rs,0.4450);
pr('Rr'' [ohm]',r.Rr,NaN);
fprintf('  %-24s %12.4f %12s   (source : Lk.Xs_slot)\n','Xsigma_s [ohm]',ctx.Lk.Xs_slot,'--');
fprintf('  %-24s %12.4f %12s   (source : r.Xr)\n','Xsigma_r [ohm]',r.Xr,'--');
pr('Xm sature [ohm]',r.Xm,46.0);
pr('Rfe [ohm]',r.Rfe,1740);
fprintf('  X_m NON sature : %.3f ohm\n',ctx.Xm0);
fprintf('  (Rr'' et Xr sont identiques a la configuration de production :\n');
fprintf('   cage.m:61 donne asq = asq1 au fondamental dans les deux cas)\n');

%% ---- 2. point nominal ------------------------------------------------
fprintf('\n  ---- 2. POINT NOMINAL ----\n');
fprintf('  %-24s %12s %12s %10s\n','grandeur','P1 skew=0','EF','ecart');
q=@(l,a,f)fprintf('  %-24s %12.4f %12.4f %9.2f %%\n',l,a,f,100*(a-f)/f);
q('couple [N.m]',r.Tem,121.63);
q('courant I1 [A]',r.I1,19.73);
q('cos(phi)',r.cosphi,0.8640);
q('pertes fer [W]',r.Pfe,232.6);
fprintf('  %-24s %12.4f %12.4f %8.2f pp\n','rendement',r.eta,0.9160,100*(r.eta-0.9160));
fprintf('  pertes Joule rotor [W]   %12.4f   (dont harmoniques %.4f)\n',r.Pcu_r,r.Pcu_r_h);

%% ---- 3. a vide et calage ---------------------------------------------
r0 = mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
rb = mec.equivalent_circuit(ctx,1.0 ,ctx.Xm0);
kbar = (2*M.m*ctx.W.kw1*ctx.W.Nph)/M.Nr;
fprintf('\n  ---- 3. A VIDE ----\n');
q('  I magnetisant [A]',r0.I1,8.49);
q('  f.e.m. entrefer [V]',r0.E1,382.1);
q('  pertes fer [W]',r0.Pfe,249.3);            %  <<< la cellule manquante
fprintf('\n  ---- 4. ROTOR BLOQUE (diagnostic, NON publie en tableau) ----\n');
q('  couple [N.m]',rb.Tem,104.31);
q('  I1 [A]',rb.I1,108.21);
q('  I barre [A]',rb.I2*kbar,1929);
fprintf('  pertes Joule rotor [W]   %12.4f   (dont harmoniques %.4f)\n',rb.Pcu_r,rb.Pcu_r_h);
fprintf('  RAPPEL : a s = 1 tous les glissements harmoniques valent 1 et\n');
fprintf('  la cage harmonique gouverne le resultat. Ce point est traite en\n');
fprintf('  section 5.6 de l''Article II, hors du tableau de comparaison.\n');

%% ---- 5. champs a mi-entrefer -----------------------------------------
fprintf('\n  ---- 5. CHAMPS A MI-ENTREFER ----\n');
Rm=0.5*(G.Rs+G.Rr); thq=linspace(0,2*pi,2001); thq(end)=[];
Rnl=mec.magnetizing(ctx,r0.Im);
[Br0,Bt0]=A1.field(Rnl.S.Usurf,Rm,thq); Br0=Br0(:).'; Bt0=Bt0(:).';
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

%% ---- 6. caracteristique ----------------------------------------------
fprintf('\n  ---- 6. CARACTERISTIQUE (30 glissements) ----\n');
sl=[0.005:0.005:0.12,0.15,0.2,0.3,0.5,0.7,1.0];
T=nan(size(sl)); I=nan(size(sl)); Xp=ctx.Xm0; b6=0;
for k=1:numel(sl)
    rk=mec.equivalent_circuit(ctx,sl(k),Xp); Xp=rk.Xm;
    T(k)=rk.Tem; I(k)=rk.I1;
    c=mec.power_balance(rk,M); b6=max(b6,c.err_global);
end
[Tmx,im]=max(T);
fprintf('  bilan de puissance : %.2e sur %d glissements\n',b6,numel(sl));
fprintf('  decrochage : %.1f N.m a s = %.3f\n',Tmx,sl(im));
fprintf('  (reference : 324.9 N.m a s = 0.105, bande LISSE du balayage --\n');
fprintf('   rugosite locale 0.26 %% pour s <= 0.05, 0.64 %% pour s <= 0.15)\n');
fprintf('  %6s %10s %10s\n','s','T [N.m]','I1 [A]');
for k=[1 4 10 20 24 numel(sl)]
    fprintf('  %6.3f %10.2f %10.2f\n',sl(k),T(k),I(k));
end
save('B10_b1_skewoff.mat','sl','T','I','r','r0','rb');
fprintf('\n  duree %.0f s\n=== B10 termine ===\n',toc(t0));
diary off;
