%% RUN_ARTICLE  -  FICHIER MAITRE : comparaison ANALYTIQUE / MEC / FEM (ANSYS Maxwell 2D)
%
%  Moteur asynchrone a cage 18,5 kW / 690 V / 4 poles / 48 encoches / 44 barres,
%  toles M800-50A, geometrie d'encoches CONFORME AUX FIGURES ANSYS (stator tc6 ;
%  rotor a epaulement, cercle superieur tangent a l'ouverture).
%
%  Ce script UNIQUE regenere l'ensemble des resultats, figures et tableaux de
%  l'etude (niveau presentation scientifique), en s'appuyant sur le paquet +mec.
%  L'etude comparative confronte TROIS approches sur la MEME machine :
%    * ANALYTIQUE : methode de conception Pyrhonen-Jokinen-Hrabovcova
%      (analytical_18_5kW, version fonction de main_motor_IM.m) ;
%    * MEC        : reseau de reluctances non lineaire (paquet +mec) ;
%    * FEM        : reference ANSYS Maxwell 2D (mec.ansys_ref + essais).
%
%    I.    Machine et geometrie (cotes derivees, controles)
%    II.   Regime etabli : parametres du schema, T(s), I(s), rendement,
%          cosphi, pertes  ->  ANALYTIQUE vs MEC vs FEM
%    III.  Validation multi-essais : en charge / a vide / rotor bloque
%    IV.   Champs : profils Br/Bt au mi-entrefer, spectres spatiaux,
%          distribution du flux, cartes 2D |B|, sondes locales vs FEM
%    V.    TRANSITOIRES MONO-TRANCHE (single slice) : la carte d'ondulation
%          de denture Delta_T(theta) du champ MST a UNE tranche (machine non
%          vrillee, homologue du modele EF) est superposee a l'equation
%          mecanique dq -> vitesse, couple AVEC ondulation, courant
%          statorique, zoom de regime etabli et caracteristique
%          couple-vitesse, a VIDE et EN CHARGE (TL(t) releve d'ANSYS)
%    VI.   Bilan des pertes et pertes supplementaires (381,5 W)
%    VII.  Synthese finale
%
%  Les transitoires LISSES (dq seul) restent dans RUN_TRANSIENT ; l'etude
%  multi-tranches (vrillage) est consignee dans ETUDE_CRITIQUE_ANSYS_MEC.md.
%
%  POURQUOI les autres RUN_* ne sont PAS re-executes ici : ce fichier est
%  l'orchestrateur de PRESENTATION (validation + resultats). Les bancs
%  d'IDENTIFICATION (RUN_FEM_IDENT, RUN_SLOTLEAK, RUN_P4_XSIGMA) exigent la
%  PDE Toolbox et ecrivent des CONSTANTES gelees dans +mec (tracees en §0) :
%  les re-executer a chaque presentation serait circulaire (l'identification
%  doit rester separee de la validation). Les DIAGNOSTICS (RUN_LEAKAGE,
%  RUN_HARMONICS, RUN_BARSKIN, RUN_P7_NU3, RUN_C1, RUN_STEPPING, RUN_MESH)
%  ont repondu a des questions aujourd'hui closes et documentees ; leurs
%  conclusions sont INTEGREES au modele (+mec) que ce script exerce.
%
%  Modele : reseau de reluctances Newton exact + entrefer harmonique de
%  Laplace (P1) + cage a effet de peau exact barre trapezoidale et anneau
%  (P2) + courants de cage harmoniques (P3) + fuites identifiees FEM +
%  extremite 3D (P4) + vrillage (C5). Bilan de puissance B6 a chaque point.
%  Duree totale ~5 min (dont ~2,5 min pour les cartes d'ondulation de la
%  section V et ~1,5 min pour les integrations temporelles a pas fin).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
OUT=fileparts(mfilename('fullpath'));
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
t_all=tic;

fprintf('======================================================================\n');
fprintf('   ANALYTIQUE / MEC / FEM — moteur asynchrone 18,5 kW (RUN_ARTICLE)\n');
fprintf('======================================================================\n');

M=mec.machine_18_5kW(); ctx=mec.build_context(M); ref=mec.ansys_ref();
G=ctx.G; W=ctx.W; BH=ctx.BH; p=M.p; Nr=M.Nr; Ns=M.Ns; mu0=4*pi*1e-7;
% --- Methode ANALYTIQUE (Pyrhonen-Jokinen-Hrabovcova) : 3e voie de comparaison ---
%     version fonction de main_motor_IM.m (requiert Function_magn/Function_magnd)
addpath('C:\Users\hp\Desktop\Matlab program\conception\New Folder');
AN=analytical_18_5kW();
%  §5.5 — PLUS DE DISTINCTION M / Mf. L'entrefer harmonique n'est plus une
%  option a activer : mec.mesh_refined le prend par defaut, et le reseau de
%  performance passe par l'operateur condense (mec.build_context). M suffit.
%  ctxC est le MODELE DE COMPARAISON, et ne sert qu'au tableau de cloture
%  §VII : meme geometrie, meme materiau, meme bobinage, meme solveur — seule
%  la fermeture d'entrefer change (coefficient de Carter + normalisation des
%  permeances). C'est l'etat anterieur du modele, conserve pour etre mesure.
%  La comparaison se fait au niveau du RESEAU, ou les deux fermetures sont
%  strictement substituables ; il n'y a donc pas de variante de maillage.
ctxC=ctx; ctxC.AG=ctx.AGcarter;
thbar=2*pi*(0:Nr-1)'/Nr; Rmid=0.5*(G.Rs+G.Rr);
thq=linspace(0,2*pi,2001); thq(end)=[];
rd=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
rmsw=@(A,c,t0)sqrt(mean(A(A(:,1)>=t0,c).^2,'omitnan'));
avgw=@(A,c,t0)mean(A(A(:,1)>=t0,c),'omitnan');
fuU=@(y,th,k)(2/numel(y))*sum(y(:).*exp(-1i*k*th(:)));

%% ============ I. MACHINE ET GEOMETRIE ============
fprintf('\n===== I. GEOMETRIE (encoches conformes aux figures ANSYS) =====\n');
%  kC / g_eff appartiennent au modele de COMPARAISON : l'operateur DtN ne
%  porte aucun de ces deux parametres. Ils sont affiches ici a titre de
%  repere geometrique, et mesures a nouveau en §7.
fprintf('entrefer g=%.3f mm  [repere Carter : kC=%.3f  g_eff=%.3f mm]\n', ...
    M.g*1e3,ctx.AGcarter.kC,ctx.AGcarter.g_eff*1e3);
fprintf('stator tc6 : bs0/bs1/bs2=%.2f/%.2f/%.2f  hs=%.3f mm  (dent parallele, err %.1e)\n',...
    M.bs0*1e3,M.bs1*1e3,M.bs2*1e3,G.hs*1e3,G.bts_parallel_err);
fprintf('rotor (epaulement, cercle sup. tangent) : br0/br1/br2=%.2f/%.2f/%.2f  hr=%.3f mm\n',...
    M.br0*1e3,M.br1*1e3,M.br2*1e3,G.hr*1e3);
fprintf('bts=%.3f  btr=%.3f mm   Abar=%.2f mm2   hys/hyr=%.2f/%.2f mm\n',...
    G.bts*1e3,G.btr*1e3,G.Abar*1e6,G.hys*1e3,G.hyr*1e3);
fprintf('bobinage : kw1=%.4f  Nph=%d  q=%d  (5/6)   Rs=%.4f ohm   Xs_end3D=%.3f ohm\n',...
    W.kw1,round(W.Nph),W.q,ctx.Rs,ctx.Lk.Xs_end3D);
% --- Tracabilite des constantes IDENTIFIEES (bancs executes UNE fois) ---
fprintf('--- Constantes identifiees du modele (provenance) ---\n');
fprintf('sigma0 = 0,53*tau        <- RUN_FEM_IDENT   (banc 2 dents, PDE Toolbox)\n');
fprintf('lam_slot = 2,722 / 1,543 <- RUN_SLOTLEAK    (bancs encoche, valides 0,0-0,2 %%)\n');
fprintf('Xs_end3D = %.3f ohm      <- RUN_P4_XSIGMA   (Xsigma implique EF, 3 glissements)\n',ctx.Lk.Xs_end3D);
fprintf('kd/ky, kfw, kcu          <- RUN_VALIDATION  (recalage pertes 16/07/2026)\n');
fprintf('P_stray = %.1f W (nom.)  <- pertes suppl. en charge, allocation IEC 60034-2-1 (loi I2^2)\n',M.stray.Pn_ref);
fprintf('J = %.2f kg.m2           <- run-up EF (RUN_TRANSIENT)\n',M.mech.J);
fprintf('(bancs et diagnostics NON re-executes ici : identification separee de la\n');
fprintf(' validation ; constantes gelees dans +mec avec commentaires de provenance)\n');

%% ============ II. REGIME ETABLI : ANALYTIQUE vs MEC vs FEM ============
fprintf('\n===== II. REGIME ETABLI : ANALYTIQUE vs MEC vs FEM =====\n');
% --- Parametres du schema equivalent : ANALYTIQUE / MEC / EF ---
rC0=mec.equivalent_circuit(ctx,0.0188,ctx.Xm0);
fprintf('--- Parametres du schema equivalent (ANALYTIQUE / MEC / EF) ---\n');
fprintf('%-16s %11s %11s %11s\n','parametre','ANALYTIQUE','MEC','EF');
fprintf('%-16s %11.4f %11.4f %11.4f\n','Rs [ohm]',AN.Rs,ctx.Rs,0.445);
fprintf('%-16s %11.4f %11.4f %11.4f\n','Rr'' [ohm]',AN.Rr,rC0.Rr,0.44);
fprintf('%-16s %11.2f %11.2f %11.1f\n','Xm [ohm]',AN.Xm,rC0.Xm,46);
fprintf('%-16s %11.3f %11.3f %11s\n','Xsigma_s [ohm]',AN.Xssigma,ctx.Lk.Xs_leak,'-');
fprintf('%-16s %11.3f %11.3f %11s\n','Xsigma_r [ohm]',AN.Xrsigma,rC0.Xr,'-');
fprintf('%-16s %11.1f %11.1f %11.0f\n','Rfe [ohm]',AN.Rfe,rC0.Rfe,1740);
fprintf('%-16s %11.2f %11.2f %11s\n','Lm [mH]',AN.Lm*1e3,rC0.Xm/M.w*1e3,'-');
s_list=[0.005:0.005:0.12,0.15,0.2,0.3,0.5,0.7,1.0]; N=numel(s_list);
Res=cell(N,1); Xp=ctx.Xm0; maxB6=0;
for k=1:N
    Res{k}=mec.equivalent_circuit(ctx,s_list(k),Xp); Xp=Res{k}.Xm;
    c=mec.power_balance(Res{k},M); maxB6=max(maxB6,c.err_global);
end
gv=@(f)cellfun(@(r)r.(f),Res);
T=gv('Tem'); I1=gv('I1'); eta=gv('eta'); PF=gv('cosphi');
Pcus=gv('Pcu_s'); Pfe=gv('Pfe'); Pcur=gv('Pcu_r'); Pfw=gv('Pfw'); Padd=gv('Padd'); nrpm=gv('n_rpm');
fprintf('Bilan de puissance B6 : erreur max = %.2e sur %d glissements\n',maxB6,N);
fprintf('%6s | %8s %8s %7s | %8s %8s %7s\n','s','T_MEC','T_FEM','err%','I_MEC','I_FEM','err%');
for s=[0.005 0.02 0.05 0.10 0.20 0.50 1.0]
    r=mec.equivalent_circuit(ctx,s,ctx.Xm0);
    Tf=interp1(ref.s,ref.T,s); If=interp1(ref.s2,ref.I,s);
    fprintf('%6.3f | %8.1f %8.1f %6.1f | %8.1f %8.1f %6.1f\n',...
        s,r.Tem,Tf,100*(r.Tem-Tf)/Tf,r.I1,If,100*(r.I1-If)/If);
end
[Tmax,imax]=max(T);
fprintf('Decrochage : MEC %.1f N.m (s=%.3f) ; Analytique %.1f (s=%.3f) ; FEM %.1f (s=%.3f)\n',...
    Tmax,s_list(imax),AN.Tb,AN.sb,ref.anchor.T_break,ref.anchor.s_break);
% --- Point nominal (s=0.0188, glissement EF en charge) : ANALYTIQUE / MEC / EF ---
sc=0.0188; rC=mec.equivalent_circuit(ctx,sc,ctx.Xm0);
TA=interp1(AN.ss,AN.T,sc); IA=interp1(AN.ss,AN.i1,sc);
etaA=interp1(AN.ss,AN.eff,sc); cosA=interp1(AN.ss,AN.cosfi,sc);
Tf=interp1(ref.s,ref.T,sc); If=interp1(ref.s2,ref.I,sc);
fprintf('\n--- Point nominal (s=%.4f) : ANALYTIQUE / MEC / EF ---\n',sc);
fprintf('%-18s %11s %11s %11s\n','grandeur','ANALYTIQUE','MEC','EF');
fprintf('%-18s %11.2f %11.2f %11.2f\n','Couple [N.m]',TA,rC.Tem,Tf);
fprintf('%-18s %11.2f %11.2f %11.2f\n','Courant I1 [A]',IA,rC.I1,If);
fprintf('%-18s %11.4f %11.4f %11.4f\n','Rendement',etaA,rC.eta,0.916);
fprintf('%-18s %11.4f %11.4f %11.4f\n','cos(phi)',cosA,rC.cosphi,0.864);
fprintf('%-18s %11.1f %11.1f %11.1f\n','Pertes fer [W]',AN.Pfe,rC.Pfe,232.6);
fprintf('%-18s %11.1f %11.1f %11.1f\n','Pertes suppl. [W]',AN.PLL,rC.Padd,381.5);
fprintf('(les 3 methodes utilisent la MEME allocation de pertes suppl. 381,5 W, IEC 60034-2-1 ;\n');
fprintf(' rendement MEC %.1f %% > plage exp. 90,5-91,3 %% : ses pertes identifiees sont ~215 W sous l''EF)\n',100*rC.eta);

figure('Name','II. Caracteristiques (Analytique/MEC/FEM)','Position',[40 40 1150 700]);
%  BLOC A8 / decision 1 du 6 aout 2026. Le panneau (a) tracait les 30
%  points TRANSCRITS de mec.ansys_ref. Il trace desormais les 129 points
%  BRUTS du balayage : les deux bandes bruitees deviennent visibles --
%  0,135-0,21 (13,2 % rms, etendue 247,6 a 336,6 N.m) et s >= 0,90
%  (9,9 % rms, ecart-type 16,2 % sur 21 points) -- ce qui justifie devant
%  le lecteur que l'usage quantitatif soit borne a s <= 0,13.
tcA=rd(fullfile(ROOT,'caratéristique en fonction glissement','Torque Plot 2.tab'));
subplot(2,2,1);
plot(tcA(:,1),tcA(:,2),'.','Color',[.88 .60 .60],'MarkerSize',7); hold on;
mq=tcA(:,1)<=0.13;
plot(tcA(mq,1),tcA(mq,2),'r-','LineWidth',1.2);
plot(s_list,T,'b-o','LineWidth',1.3,'MarkerSize',3);
plot(AN.ss,AN.T,'g-','LineWidth',1.4); grid on;
xline(0.13,'k:','LineWidth',1);
xlabel('slip s'); ylabel('T_{em} [N·m]');
legend('FEA, 129 raw points','FEA, quantitative band','MEC','Analytical', ...
    'Location','best','FontSize',7);
title('(a) Torque-slip: FEA raw sweep, quantitative for s \leq 0.13');
subplot(2,2,2); plot(s_list,I1,'b-o','LineWidth',1.3); hold on;
plot(ref.s2,ref.I,'r--s','LineWidth',1.0);
plot(AN.ss,AN.i1,'g-','LineWidth',1.4); grid on;
xlabel('slip s'); ylabel('I_1 [A]');
legend('MEC','FEA','Analytical','Location','best'); title('(b) Current-slip');
subplot(2,2,3); plot(s_list,eta,'b-o',s_list,PF,'m-^','LineWidth',1.2); hold on;
plot(AN.ss,AN.eff,'g-',AN.ss,AN.cosfi,'g--','LineWidth',1.2);
grid on; ylim([0 1]); xlim([0 0.2]); xlabel('slip s');
legend('\eta MEC','cos\phi MEC','\eta Analyt.','cos\phi Analyt.','Location','best');
title('(c) Efficiency and power factor');
subplot(2,2,4); area(s_list,[Pcus(:),Pfe(:),Pcur(:),Pfw(:),Padd(:)]); grid on;
xlabel('slip s'); ylabel('losses [W]');
legend('Stator Cu','iron','Rotor Cu','mech.','stray','Location','northwest');
title('(d) Loss breakdown (MEC)');
savefig2(gcf,OUT,'article_II_caracteristiques');

%% ============ III. VALIDATION MULTI-ESSAIS ============
fprintf('\n===== III. VALIDATION MULTI-ESSAIS (moyennes de regime etabli) =====\n');
% --- en charge ---
D=fullfile(ROOT,'transitoire','en charge'); t0=1.0;
spd=rd(fullfile(D,'Speed Plot 1.tab'));  n_ans=avgw(spd,2,t0);
s_ch=(M.ns*60-n_ans)/(M.ns*60);
trq=rd(fullfile(D,'Plot 1.tab'));        pw=rd(fullfile(D,'Output Variables Plot 3.tab'));
cl=rd(fullfile(D,'Loss Plot 1.tab'));    sl=rd(fullfile(D,'Loss Plot 2.tab'));
rng_=rd(fullfile(D,'End Connection Plot 3.tab')); ec=rd(fullfile(D,'End Connection Plot 1.tab'));
cur=rd(fullfile(D,'Winding Plot 4.tab'));
r=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
fprintf('--- EN CHARGE (n=%.1f tr/min -> s=%.4f) ---\n',n_ans,s_ch);
fprintf('%-26s %10s %10s %8s\n','grandeur','MEC','FEM','ecart%');
pr2('Couple [N.m]',r.Tem,avgw(trq,3,t0));
pr2('I1 [A]',r.I1,rmsw(cur,2,t0));
pr2('I_barre [A]',r.Ibar,1e3*rmsw(ec,2,t0));
pr2('I_anneau [A]',r.Iring,1e3*rmsw(ec,3,t0));
pr2('P_in [W]',r.Pin,avgw(pw,4,t0)*1e3);
pr2('Pertes fer [W]',r.Pfe,avgw(cl,2,t0));
pr2('Pertes Cu stator [W]',r.Pcu_s,avgw(sl,3,t0)*1e3);
pr2('Pertes rotor [W]',r.Pcu_r,(avgw(sl,2,t0)+avgw(rng_,2,t0))*1e3);
% --- a vide ---
D0=fullfile(ROOT,'transitoire','a vide');
r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
cur0=rd(fullfile(D0,'Winding Plot 4.tab')); vi0=rd(fullfile(D0,'Winding Plot 2.tab'));
cl0=rd(fullfile(D0,'Loss Plot 1.tab'));
fprintf('--- A VIDE ---\n');
pr2('I0 [A]',r0.I1,rmsw(cur0,2,t0));
pr2('E1 [V]',r0.E1,rmsw(vi0,2,t0));
pr2('Pertes fer [W]',r0.Pfe,avgw(cl0,2,t0));
% --- rotor bloque ---
Db=fullfile(ROOT,'transitoire','rotor bloqué');
curb=rd(fullfile(Db,'Winding Plot 3.tab')); trqb=rd(fullfile(Db,'Torque Plot 1.tab'));
ecb=rd(fullfile(Db,'End Connection Plot 1.tab')); t0b=0.5*max(curb(:,1));
rb=mec.equivalent_circuit(ctx,1.0,ctx.Xm0);
fprintf('--- ROTOR BLOQUE ---\n');
pr2('Couple [N.m]',rb.Tem,avgw(trqb,2,t0b)*1e3);
pr2('I1 [A]',rb.I1,rmsw(curb,2,t0b));
pr2('I_barre [A]',rb.Ibar,1e3*rmsw(ecb,2,t0b));
pr2('I_anneau [A]',rb.Iring,1e3*rmsw(ecb,3,t0b));

%% ============ IV. CHAMPS : PROFILS, SPECTRES, CARTES, SONDES ============
fprintf('\n===== IV. CHAMPS vs FEM (entrefer harmonique de Laplace, P1) =====\n');
P0=read_profile(fullfile(ROOT,'transitoire','a vide','Calculator Expressions Plot 1.tab'));
P1p=read_profile(fullfile(ROOT,'transitoire','en charge','Calculator Expressions Plot 1.tab'));
n1=1469.772776; s1=(M.ns*60-n1)/(M.ns*60);
r1f=mec.equivalent_circuit(ctx,s1,ctx.Xm0);
[i30,ib0]=inst_currents(M,W,r0);
[i31,ib1]=inst_currents(M,W,r1f);
me0=mec.mesh_refined(M,G,6,i30,0,3,ib0); Se0=mec.solve_mesh(me0,BH,M.opt);
me1=mec.mesh_refined(M,G,6,i31,0,3,ib1); Se1=mec.solve_mesh(me1,BH,M.opt);
[Us0,Ur0]=surfU(me0,Se0); [Us1,Ur1]=surfU(me1,Se1);
[Br0,Bt0]=me0.gapF.field(Us0,Ur0,Rmid,thq);
[Br1,Bt1]=me1.gapF.field(Us1,Ur1,Rmid,thq);
fprintf('Bg1 mi-entrefer : vide %.3f (FEM %.3f) ; charge %.3f (FEM %.3f) T\n',...
    abs(fuU(Br0,thq,p)),abs(fuU(P0.Br,P0.th,p)),abs(fuU(Br1,thq,p)),abs(fuU(P1p.Br,P1p.th,p)));
fprintf('Bt RMS          : vide %.3f (FEM %.3f) ; charge %.3f (FEM %.3f) T\n',...
    rms(Bt0),rms(P0.Bt),rms(Bt1),rms(P1p.Bt));
fprintf('Flux/pole       : vide %.2f (FEM %.2f) mWb\n',...
    (max(cumsum(Br0)*mean(diff(thq))*Rmid)-min(cumsum(Br0)*mean(diff(thq))*Rmid))*M.L*1e3,...
    (max(P0.A)-min(P0.A))*M.L*1e3);

figure('Name','IV. Champs entrefer','Position',[50 50 1150 700]);
subplot(2,2,1);
plot(P0.th*180/pi,P0.Br,'r-','LineWidth',0.8); hold on;
plot(al(thq,Br0,P0,p)*180/pi,Br0,'b.','MarkerSize',3);
grid on; xlim([0 90]); xlabel('\theta [deg]'); ylabel('B_r [T]');
legend('FEA','MEC'); title('(a) B_r at mid-gap, no load');
subplot(2,2,2);
plot(P0.th*180/pi,P0.Bt,'r-','LineWidth',0.8); hold on;
plot(al(thq,Bt0,P0,p)*180/pi,Bt0,'b.','MarkerSize',3);
grid on; xlim([0 90]); xlabel('\theta [deg]'); ylabel('B_t [T]');
legend('FEA','MEC'); title('(b) B_t at mid-gap, no load (P1)');
nus=[1 5 7 11 13 17 19 23 25];
Am=arrayfun(@(n)abs(fuU(Br0,thq,n*p)),nus);
Aa=arrayfun(@(n)abs(fuU(P0.Br,P0.th,n*p)),nus);
subplot(2,2,3);
bar(nus,[Aa(:) Am(:)],1.0); set(gca,'YScale','log'); ylim([1e-4 2]); grid on;
xlabel('harmonic order \nu (electrical)'); ylabel('|B_{r,\nu}| [T]');
legend('FEA','MEC'); title('(c) Spatial spectrum of B_r (no load)'); xticks(nus);
subplot(2,2,4);
plot(P1p.th*180/pi,(P1p.A-mean(P1p.A))*1e3,'r-','LineWidth',1.1); hold on;
Ac=cumsum(Br1)*mean(diff(thq))*Rmid;
plot(al(thq,Br1,P1p,p)*180/pi,(Ac-mean(Ac))*1e3,'b.','MarkerSize',3);
grid on; xlabel('\theta [deg]'); ylabel('A [mWb/m]');
legend('FEA','MEC'); title('(d) Flux distribution, on load');
savefig2(gcf,OUT,'article_IV_champs');

% sondes locales (cartes .docx ANSYS) vs maxima regionaux MEC
[Bts0,Bys0,Btr0,Byr0]=regional_max(me0,Se0);
[Bts1,Bys1,Btr1,Byr1]=regional_max(me1,Se1);
fprintf('--- Sondes locales FEM vs maxima regionaux MEC [T] ---\n');
pr2('culasse stator (vide)',Bys0,1.837); pr2('dent stator (vide)',Bts0,1.633);
pr2('dent rotor (vide)',Btr0,1.929);     pr2('culasse rotor (vide)',Byr0,1.674);
pr2('culasse stator (chg)',Bys1,1.896);  pr2('dent stator (chg)',Bts1,1.684);
pr2('dent rotor (chg)',Btr1,1.870);      pr2('culasse rotor (chg)',Byr1,1.577);
% cartes 2D STRUCTUREES (geometrie reelle des encoches, lignes de flux =
% contours de la fonction de flux psi integree des flux de branches, et
% carte de densite de courant J) — rendu comparable aux cartes ANSYS,
% dans les DEUX conditions : a VIDE et EN CHARGE.
%  LANGUE DES FIGURES : ANGLAIS. Titres, libelles d'echelle et noms de
%  fenetre sont passes ICI, et non figes dans mec.draw_cross_section : la
%  langue d'affichage est une decision de la figure, pas de la fonction de
%  trace. Les chaines sont volontairement en ASCII pur -- ces figures
%  partent dans l'article, et aucun caractere ne depend alors de
%  l'encodage du fichier .m ni de la police du moteur de rendu.
%  Les NOMS DE FICHIER restent inchanges (article_IV_carte_*) : ils sont
%  cites par le .tex, et libelle affiche n'est pas nom de fichier.
figEN=struct( ...
    'titleB','Flux density |B| and flux lines (actual slot geometry)', ...
    'titleJ','Instantaneous current density (bars + two stator layers)', ...
    'labB','|B| [T]', ...
    'labJ','J [A/m^2]', ...
    'nameB','Flux density map |B|', ...
    'nameJ','Current density map J');

o0=figEN; o0.prefix='No load: '; o0.Jmax=3e6;
mec.draw_cross_section(me0,Se0,M,G,W,i30,ib0,o0);
figs=get(groot,'Children');
%  raster = true : cartes de champ DENSES, seule exception a l'export
%  vectoriel (bloc A8). PNG 600 dpi.
savefig2(figs(2),OUT,'article_IV_carte_B_avide',true);
savefig2(figs(1),OUT,'article_IV_carte_J_avide',true);

o1=figEN; o1.prefix='On load: ';
mec.draw_cross_section(me1,Se1,M,G,W,i31,ib1,o1);
figs=get(groot,'Children');
savefig2(figs(2),OUT,'article_IV_carte_B_charge',true);
savefig2(figs(1),OUT,'article_IV_carte_J_charge',true);

%% ============ V. TRANSITOIRES MONO-TRANCHE (ondulation incluse) ============
fprintf('\n===== V. TRANSITOIRES MONO-TRANCHE (ondulation de denture) =====\n');
Npos=19; Pth=2*(2*pi/Ns); th0=linspace(0,Pth,Npos);
ptsT={{'a vide',r0},{'en charge',r1f}};
DTmap=cell(2,1);
tic;
for j=1:2
    rj=ptsT{j}{2};
    psi1=angle(rj.I1c); psi2=angle(rj.I2c);
    Ibj=rj.I2*(2*M.m*W.kw1*W.Nph)/Nr;
    Tk=zeros(Npos,1);
    for k=1:Npos
        ph=p*th0(k);
        i3k=sqrt(2)*rj.I1*[cos(psi1+ph);cos(psi1+ph-2*pi/3);cos(psi1+ph+2*pi/3)];
        ibk=sqrt(2)*Ibj*cos(p*thbar-psi2);
        me=mec.mesh_refined(M,G,6,i3k,th0(k),3,ibk);
        Se=mec.solve_mesh(me,BH,M.opt);
        [Us,Ur]=surfU(me,Se);
        Tk(k)=me.gapF.torque(Us,Ur);
    end
    dTk=Tk-mean(Tk);
    DTmap{j}=@(th) interp1(th0,dTk,mod(th,Pth),'linear');
    fprintf('carte %-9s : ondulation mono-tranche = %.1f N.m crete-a-crete\n',...
        ptsT{j}{1},max(Tk)-min(Tk));
end
fprintf('(2 x %d resolutions de champ en %.0f s)\n',Npos,toc);
o0=mec.dq_startup(ctx,struct('tend',0.4,'TL',0,'Trip',DTmap{1}));
oc=mec.dq_startup(ctx,struct('tend',2.0,'TL',[trq(:,1),trq(:,2)],'Trip',DTmap{2}));
spd0=rd(fullfile(D0,'la vitesse en fonction du temps.tab'));
trq0=rd(fullfile(D0,'Torque Plot 1.tab'));
ppv=@(x)max(x)-min(x);
w0=o0.t>=0.35; wc=oc.t>=1.90;
f0=trq0(:,1)>=1.90; fc=trq(:,1)>=1.90;
s0v=spd0(:,1)>=1.90; scv=spd(:,1)>=1.90;
fprintf('--- Regime etabli : ondulation crete-a-crete (MEC 1 tranche vs FEM) ---\n');
fprintf('%12s | %11s %11s | %13s %13s\n','essai','dT MEC','dT FEM','dn MEC','dn FEM*');
fprintf('%12s | %8.1f Nm %8.1f Nm | %9.3f rpm %9.2f rpm\n','a vide',...
    ppv(o0.Tem(w0)),ppv(trq0(f0,3)),ppv(o0.n_rpm(w0)),ppv(spd0(s0v,2)));
fprintf('%12s | %8.1f Nm %8.1f Nm | %9.3f rpm %9.2f rpm\n','en charge',...
    ppv(oc.Tem(wc)),ppv(trq(fc,3)),ppv(oc.n_rpm(wc)),ppv(spd(scv,2)));
fprintf('(* dn FEM exagere par le repli du pas de 1 ms ; a vide, dT MEC ~x8 :\n');
fprintf('   l''amortissement de denture par la cage est absent du banc de champ)\n');
fprintf('EN CHARGE : t95 = %.3f s (FEM %.3f) ; moyennes n = %.1f (FEM %.1f) tr/min,\n',...
    oc.t(find(oc.n_rpm>=0.95*1500,1)),spd(find(spd(:,2)>=0.95*1500,1),1),...
    mean(oc.n_rpm(wc)),mean(spd(scv,2)));
fprintf('T = %.1f (FEM %.1f) N.m\n',mean(oc.Tem(wc)),mean(trq(fc,3)));

%  lab   : sert de NOM DE FICHIER (reference par l'article) -> inchange
%  labEN : sert de LIBELLE AFFICHE -> anglais, comme le reste des figures
casesT={{'a vide',o0,spd0,trq0,cur0,0.4,0.35,'no load'}, ...
        {'en charge',oc,spd,trq,cur,2.0,1.90,'on load'}};
for j=1:2
    lab=casesT{j}{1}; o=casesT{j}{2}; sp=casesT{j}{3}; tq=casesT{j}{4};
    cu=casesT{j}{5}; tendj=casesT{j}{6}; tz=casesT{j}{7}; labEN=casesT{j}{8};
    figure('Name',['V. single-slice transient - ',labEN],'Position',[40 30 1250 720]);
    subplot(2,3,1);
    plot(o.t,o.n_rpm,'b-','LineWidth',1.2); hold on;
    plot(sp(:,1),sp(:,2),'r--','LineWidth',0.9);
    grid on; xlim([0 tendj]); xlabel('t [s]'); ylabel('n [rpm]');
    legend('MEC single slice','FEA','Location','southeast'); title(['(a) Speed - ',labEN]);
    subplot(2,3,2);
    plot(o.t,o.Tem,'b-','LineWidth',0.6); hold on;
    plot(tq(:,1),tq(:,3),'r--','LineWidth',0.6);
    grid on; xlim([0 tendj]); xlabel('t [s]'); ylabel('T_{em} [N·m]');
    legend('MEC single slice','FEA'); title('(b) Torque (with ripple)');
    subplot(2,3,3);
    plot(o.t,o.ia,'b-','LineWidth',0.6); hold on;
    plot(cu(:,1),cu(:,2),'r--','LineWidth',0.6);
    grid on; xlim([0 tendj]); xlabel('t [s]'); ylabel('i_A [A]');
    legend('MEC','FEA'); title('(c) Stator current');
    subplot(2,3,4);
    mM=o.t>=tz; mF=tq(:,1)>=tz & tq(:,1)<=tz+0.05;
    plot((tq(mF,1)-tz)*1e3,tq(mF,3),'r-','LineWidth',1.0); hold on;
    plot((o.t(mM)-tz)*1e3,o.Tem(mM),'b-','LineWidth',1.0);
    plot((o.t(mM)-tz)*1e3,o.Tem1(mM),'k:','LineWidth',1.1);
    grid on; xlim([0 30]); xlabel('t [ms]'); ylabel('T_{em} [N·m]');
    legend('FEA (1 kHz aliasing)','MEC single slice','MEC dq only','Location','best');
    title('(d) Zoom: slot-harmonic ripple');
    subplot(2,3,[5 6]);
    plot(o.n_rpm,o.Tem,'b-','LineWidth',0.5); hold on;
    nF=interp1(sp(:,1),sp(:,2),tq(:,1),'linear','extrap');
    plot(nF,tq(:,3),'r--','LineWidth',0.5);
    grid on; xlabel('n [rpm]'); ylabel('T_{em} [N·m]');
    legend('MEC single slice','FEA','Location','best');
    title('(e) Torque-speed characteristic (dynamic trajectory)');
    savefig2(gcf,OUT,sprintf('article_V_1tranche_%s',strrep(lab,' ','_')));
end

%% ============ VI. BILAN DES PERTES ET PERTES SUPPLEMENTAIRES ============
fprintf('\n===== VI. BILAN DES PERTES ET PERTES SUPPLEMENTAIRES (en charge, s=%.4f) =====\n',s_ch);
sumF=avgw(cl,2,t0)+avgw(sl,3,t0)*1e3+avgw(sl,2,t0)*1e3+avgw(rng_,2,t0)*1e3+98.2;
PinF=avgw(pw,4,t0)*1e3; PoutF=avgw(pw,3,t0)*1e3; strayF=PinF-PoutF-sumF;
fprintf('FEM : Pin-Pout = %.0f W ; pertes listees = %.0f W ; NON identifiees = %.1f W (%.1f %% Pn, IEC 60034-2-1)\n',...
    PinF-PoutF,sumF,strayF,100*strayF/M.Pn);
fprintf('MEC : bilan ferme a %.1e (B6) ; Cu_s %.0f, fer %.0f, rotor %.0f (barres %.0f + anneaux %.0f), meca %.0f W\n',...
    maxB6,r.Pcu_s,r.Pfe,r.Pcu_r,r.Pbars,r.Pring,r.Pfw);
fprintf('      pertes SUPPLEMENTAIRES en charge = %.1f W (%.1f %% Pn ; loi I2^2 plafonnee, mec.stray_losses)\n',...
    r.Padd,100*r.Padd/M.Pn);
fprintf('ANA : Cu_s %.0f, fer %.0f, rotor %.0f, meca %.0f, suppl. %.0f W (PLL, IEC 60034-2-1)\n',...
    AN.Pscu,AN.Pfe,AN.Prcu,AN.Pro,AN.PLL);
fprintf('      -> rendement : ANALYTIQUE %.2f %% / MEC %.2f %% / EF %.2f %% (plage experimentale 90,5-91,3 %%)\n',...
    100*AN.eta,100*r.eta,100*PoutF/PinF);
fprintf('      Note : pertes suppl. = allocation IEC 60034-2-1 (381,5 W) COMMUNE aux 3 methodes ; la physique\n');
fprintf('      de denture (RUN_STEPPING) n''en explique que ~17,5 W (surface/inter-barres/3D 2D-inaccessibles).\n');

%% ====== VII. CLOTURE D'ENTREFER : CARTER vs OPERATEUR DtN vs EF ======
%  Tableau de cloture. Les deux colonnes MEC ne different QUE par la façon
%  de fermer l'entrefer ; geometrie, materiau, bobinage et solveur sont
%  identiques. La colonne Carter est l'etat anterieur du modele.
%
%  A LIRE AVEC LA MISE EN GARDE SUIVANTE. Dans le modele de Carter, les
%  permeances sont normalisees pour que leur somme par dent vaille
%  exactement mu0*L*tau_s/g_eff : Xm y est juste PAR CONSTRUCTION, et son
%  accord avec l'EF ne valide rien. L'operateur ne contient aucune
%  normalisation de ce type — ses ecarts sont des erreurs de PREDICTION.
%  Un chiffre qui se degrade en passant de gauche a droite n'est donc pas
%  une regression du modele : c'est le prix du retrait d'un ajustement.
fprintf('\n===== VII. CLOTURE D''ENTREFER : CARTER vs DtN vs EF =====\n');
ctxC.Xm0=mec.magnetizing(ctxC,0.2).Xm;
cA0=mec.equivalent_circuit(ctxC,1e-4,ctxC.Xm0);
cAL=mec.equivalent_circuit(ctxC,s_ch,ctxC.Xm0);
[bgA,~ ]=bg1load(ctxC,cAL);
[bgB,BtB]=bg1load(ctx ,r   );
I0f=rmsw(cur0,2,t0); Tf7=avgw(trq,3,t0);
bg1f=abs(fuU(P1p.Br,P1p.th,p)); btf=rms(P1p.Bt);
fprintf('%-24s %11s %11s %11s   %s\n','grandeur','CARTER','DtN','EF','ecart DtN');
p7=@(lab,a,b,f,u)fprintf('%-24s %11.3f %11.3f %11.3f   %+7.1f %%  %s\n', ...
    lab,a,b,f,100*(b-f)/f,u);
p7('Xm non sature [ohm]',ctxC.Xm0,ctx.Xm0,NaN,'');
p7('Xm charge [ohm]',cAL.Xm,r.Xm,46,'');
p7('I0 a vide [A]',cA0.I1,r0.I1,I0f,'');
p7('E1 a vide [V]',cA0.E1,r0.E1,382.1,'');
p7('Couple charge [N.m]',cAL.Tem,r.Tem,Tf7,'');
p7('Bg1 charge [T]',bgA,bgB,bg1f,'(flux/pas)');
fprintf('%-24s %11s %11.3f %11.3f   %+7.1f %%  %s\n','Bt rms charge [T]', ...
    'AUCUN',rms(Bt1),btf,100*(rms(Bt1)-btf)/btf,'(maillage)');
fprintf('%-24s %11.3f %11.3f %11s\n','Carter kC (implicite)', ...
    ctx.AGcarter.kC,kc_implicite(ctx,M,G),'--');
fprintf(['  Bt : le modele de Carter n''en produit AUCUN — ses permeances\n' ...
         '  sont purement radiales. La ligne n''existait pas avant.\n']);
fprintf(['  kC : 1.266 impose par formule a gauche ; a droite il EMERGE de\n' ...
         '  la resolution de la couronne, aucun coefficient n''est entre.\n']);

%% ============ VIII. SYNTHESE ============
fprintf('\n===== VIII. SYNTHESE =====\n');
fprintf('Comparaison 3 voies (sec. II) ANALYTIQUE / MEC / EF : schema equivalent, T(s), I(s),\n');
fprintf('  rendement (ANA %.1f / MEC %.1f / EF 91,6 %%) et pertes ; courbe article_II.\n',100*AN.eta,100*r.eta);
fprintf('Geometrie conforme figures ANSYS (rotor : epaulement, hr=%.2f mm, Abar=%.2f mm2).\n',G.hr*1e3,G.Abar*1e6);
fprintf('Zone utile s<=0,05 : couple < 9 %%, courant < 7 %% ; fort glissement < 7 %% ;\n');
fprintf('I0 a vide %.1f %% ; culasses locales a +/-0,1 %% (vide) ; B6 <= %.1e partout.\n',...
    100*(r0.I1-rmsw(cur0,2,t0))/rmsw(cur0,2,t0),maxB6);
fprintf('(maxima de dents MEC > sondes EF de +10..19 %% : le maillage conforme\n');
fprintf(' resout les BECS satures, les sondes EF sont au corps de dent)\n');
fprintf('Transitoires mono-tranche : ondulation en charge %.1f N.m cc (FEM %.1f) ;\n',...
    ppv(oc.Tem(wc)),ppv(trq(fc,3)));
fprintf('moyennes et courant reproduits (RUN_TRANSIENT pour la version dq lisse).\n');
fprintf('\nFigures : article_II_caracteristiques / article_IV_champs /\n');
fprintf('          article_IV_carte_B_avide / article_IV_carte_J_avide /\n');
fprintf('          article_IV_carte_B_charge / article_IV_carte_J_charge /\n');
fprintf('          article_V_1tranche_a_vide / article_V_1tranche_en_charge  (.png)\n');
fprintf('\nDuree totale : %.0f s\n=== RUN_ARTICLE termine ===\n',toc(t_all));

%% ================= fonctions locales =================
function [bg1,btr]=bg1load(cx,rr)
%  Champ d'entrefer au point de CHARGE, chemin dentaire : courant
%  statorique total + onde de FMM de barres (meme construction que
%  mec.torque_vw). La branche magnetisante seule ne porterait pas la
%  reaction d'induit, donc pas de composante tangentielle.
    Mx=cx.M; Gx=cx.G; Wx=cx.W;
    q1=angle(rr.I1c);
    i3=sqrt(2)*abs(rr.I1c)*[cos(q1);cos(q1-2*pi/3);cos(q1+2*pi/3)];
    Fs=Wx.slotMMF(i3);
    q2=angle(rr.I2c); tb=2*pi*(0:Mx.Nr-1)/Mx.Nr;
    Fu=cumsum(cos(Mx.p*tb-q2).'); Fu=Fu-mean(Fu);
    c1=(2/Mx.Nr)*sum(Fu.'.*exp(-1j*Mx.p*tb));
    Fr=-Fu*((3/2)*(4/pi)*(Wx.kw1*Wx.Nph/(2*Mx.p))*sqrt(2)*abs(rr.I2c)/abs(c1));
    S=mec.solve_network(cx.net,Gx,cx.BH,cx.AG,Fs,Fr,Mx.opt);
    th=2*pi*(0:Mx.Ns-1)/Mx.Ns; Bg=S.Bgap_avg_i(:).';
    bg1=abs((2/Mx.Ns)*sum(Bg.*exp(-1j*Mx.p*th)));
    btr=NaN;
    if isfield(S,'Usurf') && isfield(cx.AG,'expand')
        tq=linspace(0,2*pi,2001); tq(end)=[];
        [~,Bt]=cx.AG.field(S.Usurf,0.5*(Gx.Rs+Gx.Rr),tq);
        btr=sqrt(mean(Bt(:).^2));
    end
end
function kc=kc_implicite(cx,Mx,Gx)
%  Coefficient de Carter EMERGENT : rapport de la reactance a surfaces
%  lisses (arcs pavant 100 % de l'alesage) a celle de l'operateur dente.
%  Aucun coefficient n'est injecte — il est MESURE sur la solution.
    ths=2*pi*(0:Mx.Ns-1)/Mx.Ns; thr=2*pi*(0:Mx.Nr-1)/Mx.Nr;
    AL=mec.airgap_fourier(ths,repmat(2*pi/Mx.Ns,1,Mx.Ns), ...
                          thr,repmat(2*pi/Mx.Nr,1,Mx.Nr), ...
                          Gx.Rs,Gx.Rr,Mx.L,cx.AG.Nh);
    cl=cx; cl.AG=AL;
    kc=mec.magnetizing(cl,0.2).Xm/cx.Xm0;
end
function pr2(lab,vm,vf)
    fprintf('%-26s %10.2f %10.2f %+7.1f\n',lab,vm,vf,100*(vm-vf)/vf);
end
function savefig2(f,OUT,name,raster)
%SAVEFIG2  Export d'une figure d'article : PDF VECTORIEL par defaut.
%
%   BLOC A8. Les quatre figures de l'Article II etaient des gabarits dans
%   article/figures/ -- des PDF de ~20 ko portant le texte « placeholder,
%   regenerate as vector PDF from MATLAB ». Elles sont desormais produites
%   par la chaine qui les calcule.
%
%   raster = true reserve l'exception aux CARTES DE CHAMP DENSES
%   (article_IV_carte_*), qui comptent des dizaines de milliers de patches
%   et depassent la dizaine de Mo en vectoriel : PNG 600 dpi, extension a
%   changer dans le .tex. Toutes les autres sortent en .pdf vectoriel.
    if nargin<4||isempty(raster), raster=false; end
    if raster
        exportgraphics(f,fullfile(OUT,[name '.png']), ...
            'Resolution',600,'BackgroundColor','none');
    else
        exportgraphics(f,fullfile(OUT,[name '.pdf']), ...
            'ContentType','vector','BackgroundColor','none');
        exportgraphics(f,fullfile(OUT,[name '.png']),'Resolution',130);
    end
end
function [i3,ib]=inst_currents(M,W,r)
    psi1=angle(r.I1c); psi2=angle(r.I2c);
    i3=sqrt(2)*r.I1*[cos(psi1);cos(psi1-2*pi/3);cos(psi1+2*pi/3)];
    Ibar=r.I2*(2*M.m*W.kw1*W.Nph)/M.Nr;
    ib=sqrt(2)*Ibar*cos(M.p*(2*pi*(0:M.Nr-1)'/M.Nr)-psi2);
end
function [Us,Ur]=surfU(me,Se)
    Us=Se.U(me.gapF.ids(1:me.Ms)); Ur=Se.U(me.gapF.ids(me.Ms+1:end));
end
function th2=al(thq,Bmec,Pans,p)
    ca=(2/numel(Pans.Br))*sum(Pans.Br(:).*exp(-1i*p*Pans.th(:)));
    cm=(2/numel(Bmec))*sum(Bmec(:).*exp(-1i*p*thq(:)));
    th2=mod(thq(:)+(angle(cm)-angle(ca))/p,2*pi);
end
function P=read_profile(f)
    fid=fopen(f); hdr=fgetl(fid); fclose(fid);
    tk=regexp(hdr,'"([^"]+)"','tokens'); names=cellfun(@(c)c{1},tk,'UniformOutput',false);
    D=readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    gc=@(pat)D(:,find(contains(names,pat),1));
    dist=gc('Distance'); A=gc('Flux_Lines'); Br=gc('Br'); Bt=gc('Bt');
    C=dist(end)*1e-3; n=numel(dist)-1;
    P.th=2*pi*dist(1:n)*1e-3/C; P.A=A(1:n); P.Br=Br(1:n); P.Bt=Bt(1:n);
end
function [Bts,Bys,Btr,Byr]=regional_max(me,Se)
    nb=me.gapfirst-1; isFe=logical(me.iron(1:nb));
    a=me.a(1:nb); B=abs(Se.B(1:nb));
    Ms=me.Ms; Ls=me.Ls; Mr=me.Mr; nst=me.nr; nrt=me.nr;
    laS=ceil(a/Ms); laS(a>Ms*Ls)=0;
    aR=a-Ms*Ls; laR=ceil(aR/Mr); laR(a<=Ms*Ls)=0;
    Bts=max(B(isFe&laS>=1&laS<=nst)); Bys=max(B(isFe&laS>nst));
    Btr=max(B(isFe&laR>=1&laR<=nrt)); Byr=max(B(isFe&laR>nrt));
end
function map_field(me,Se,M,cmax,ttl)
    Ms=me.Ms; Mr=me.Mr; Ls=me.Ls; Lr=me.Lr;
    rs=me.rs_edges; rr=me.rr_edges;
    nb=me.gapfirst-1; a=me.a(1:nb); b=me.b(1:nb);
    Bsig=zeros(nb,1); ok=me.A(1:nb)>0; Bsig(ok)=Se.Phi(ok)./me.A(ok);
    rs_c=0.5*(rs(1:end-1)+rs(2:end)); rr_c=0.5*(rr(1:end-1)+rr(2:end));
    x=zeros(me.Nnodes,1); y=x;
    for la=1:Ls, id=(la-1)*Ms+(1:Ms); x(id)=rs_c(la)*cos(me.th_s); y(id)=rs_c(la)*sin(me.th_s); end
    for la=1:Lr, id=Ms*Ls+(la-1)*Mr+(1:Mr); x(id)=rr_c(la)*cos(me.th_r); y(id)=rr_c(la)*sin(me.th_r); end
    dx=x(b)-x(a); dy=y(b)-y(a); dl=hypot(dx,dy); dl(dl==0)=1;
    ux=dx./dl; uy=dy./dl;
    rn=hypot(x(a),y(a)); rn(rn==0)=1;
    isRad=abs(ux.*x(a)./rn+uy.*y(a)./rn)>0.7;
    N=me.Nnodes;
    acc=@(m)deal(accumarray([a(m);b(m)],[Bsig(m).*ux(m);Bsig(m).*ux(m)],[N 1]),...
                 accumarray([a(m);b(m)],[Bsig(m).*uy(m);Bsig(m).*uy(m)],[N 1]),...
                 accumarray([a(m);b(m)],1,[N 1]));
    [vxR,vyR,cR]=acc(isRad); [vxT,vyT,cT]=acc(~isRad);
    cR(cR==0)=1; cT(cT==0)=1;
    Bn=hypot(vxR./cR+vxT./cT,vyR./cR+vyT./cT);
    [XX,YY,CC]=deal([]);
    for la=1:Ls
        [X,Y]=cellpoly(me.th_s,me.dth_s,rs(la),rs(la+1));
        XX=[XX X]; YY=[YY Y]; CC=[CC Bn((la-1)*Ms+(1:Ms))']; %#ok<AGROW>
    end
    for la=1:Lr
        [X,Y]=cellpoly(me.th_r,me.dth_r,rr(la),rr(la+1));
        XX=[XX X]; YY=[YY Y]; CC=[CC Bn(Ms*Ls+(la-1)*Mr+(1:Mr))']; %#ok<AGROW>
    end
    figure('Name',ttl,'Position',[90 90 700 640]);
    patch(XX,YY,CC,'EdgeColor','none'); axis equal off;
    colormap(jet(11)); clim([0 cmax]); cb=colorbar; cb.Label.String='|B| [T]';
    title(ttl);
end
function [X,Y]=cellpoly(th,dth,r1,r2)
    th=th(:).'; dth=dth(:).';
    t1=th-dth/2; t2=th+dth/2;
    X=[r1*cos(t1);r1*cos(t2);r2*cos(t2);r2*cos(t1)];
    Y=[r1*sin(t1);r1*sin(t2);r2*sin(t2);r2*sin(t1)];
end
