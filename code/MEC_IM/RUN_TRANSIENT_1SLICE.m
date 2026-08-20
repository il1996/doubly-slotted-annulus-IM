%% RUN_TRANSIENT_1SLICE  -  Caracteristiques TRANSITOIRES en modele MONO-TRANCHE
%
%  Etude transitoire de la machine asynchrone avec le modele de champ a UNE
%  SEULE TRANCHE (single slice = machine NON VRILLEE) : l'ondulation de
%  couple de denture Delta_T(theta), extraite de la carte MST mono-tranche
%  (entrefer harmonique P1) au point de fonctionnement, est SUPERPOSEE a
%  l'equation mecanique du modele dq (mec.dq_startup, option Trip ;
%  echelle en (psi/psiN)^2 pendant l'etablissement du flux).
%
%  CONTENU (dans les DEUX conditions : a VIDE et EN CHARGE, TL(t) releve
%  d'ANSYS) :
%    * evolution temporelle de la VITESSE, du COUPLE electromagnetique
%      (fondamental + ondulation) et du COURANT statorique ;
%    * ZOOM de regime etabli mettant en evidence l'ondulation de couple du
%      modele mono-tranche, comparee a celle du transitoire EF (le modele
%      EF de reference est lui-meme non vrille -> comparaison directe) ;
%    * caracteristique COUPLE-VITESSE (trajectoire dynamique complete,
%      bande d'ondulation au point d'equilibre) ;
%    * quantification : ondulation crete-a-crete du couple et fluctuation
%      de vitesse induite, MEC mono-tranche vs FEM.
%
%  NB : une seule tranche = PAS de moyenne de vrillage : l'ondulation est
%  celle de la machine non vrillee (~5x la machine reelle vrillee).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
OUT=fileparts(mfilename('fullpath'));
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';

fprintf('=== Transitoires en modele MONO-TRANCHE (ondulation incluse) ===\n\n');
M=mec.machine_18_5kW(); ctx=mec.build_context(M);
G=ctx.G; W=ctx.W; BH=ctx.BH; p=M.p; Nr=M.Nr; Ns=M.Ns; mu0=4*pi*1e-7;
Mf=M; Mf.opt.gap_fourier=1;
thbar=2*pi*(0:Nr-1)'/Nr; Rmid=0.5*(G.Rs+G.Rr);
rd=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');

%% 1. Cartes d'ondulation MONO-TRANCHE aux deux points de fonctionnement
Npos=19; Pth=2*(2*pi/Ns); th0=linspace(0,Pth,Npos);
n1=1469.772776; s1=(M.ns*60-n1)/(M.ns*60);
r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
r1=mec.equivalent_circuit(ctx,s1,ctx.Xm0);
pts={{'a vide',r0},{'en charge',r1}};
DT=cell(2,1);
tic;
for j=1:2
    rj=pts{j}{2};
    psi1=angle(rj.I1c); psi2=angle(rj.I2c);
    Ibj=rj.I2*(2*M.m*W.kw1*W.Nph)/Nr;
    Tk=zeros(Npos,1);
    for k=1:Npos
        ph=p*th0(k);
        i3k=sqrt(2)*rj.I1*[cos(psi1+ph);cos(psi1+ph-2*pi/3);cos(psi1+ph+2*pi/3)];
        ibk=sqrt(2)*Ibj*cos(p*thbar-psi2);
        me=mec.mesh_refined(Mf,G,6,i3k,th0(k),3,ibk);
        Se=mec.solve_mesh(me,BH,Mf.opt);
        Us=Se.U(me.gapF.ids(1:me.Ms)); Ur=Se.U(me.gapF.ids(me.Ms+1:end));
        Tk(k)=me.gapF.torque(Us,Ur);
    end
    dTk=Tk-mean(Tk);
    DT{j}=@(th) interp1(th0,dTk,mod(th,Pth),'linear');
    fprintf('carte %s : ondulation mono-tranche = %.1f N.m crete-a-crete\n',...
        pts{j}{1},max(Tk)-min(Tk));
end
fprintf('(2 x %d resolutions de champ en %.0f s)\n\n',Npos,toc);

%% 2. Transitoires dq AVEC ondulation superposee
o0=mec.dq_startup(ctx,struct('tend',0.4,'TL',0,'Trip',DT{1}));
D=fullfile(ROOT,'transitoire','en charge');
trq=rd(fullfile(D,'Plot 1.tab')); spd=rd(fullfile(D,'Speed Plot 1.tab'));
cur=rd(fullfile(D,'Winding Plot 4.tab'));
oc=mec.dq_startup(ctx,struct('tend',2.0,'TL',[trq(:,1),trq(:,2)],'Trip',DT{2}));
D0=fullfile(ROOT,'transitoire','a vide');
spd0=rd(fullfile(D0,'la vitesse en fonction du temps.tab'));
trq0=rd(fullfile(D0,'Torque Plot 1.tab'));
cur0=rd(fullfile(D0,'Winding Plot 4.tab'));

%% 3. Quantification en regime etabli
pp=@(x)max(x)-min(x);
w0=o0.t>=0.35;  wc=oc.t>=1.90;               % fenetres MEC
f0=trq0(:,1)>=1.90; fc=trq(:,1)>=1.90;       % fenetres FEM
s0f=spd0(:,1)>=1.90; scf=spd(:,1)>=1.90;
fprintf('--- Regime etabli : ondulation (cc) du couple et de la vitesse ---\n');
fprintf('%12s | %14s %14s | %14s %14s\n','essai',...
    'dT MEC 1tr','dT FEM','dn MEC 1tr','dn FEM');
fprintf('%12s | %11.1f N.m %11.1f N.m | %10.3f tr/min %8.2f tr/min\n','a vide',...
    pp(o0.Tem(w0)),pp(trq0(f0,3)),pp(o0.n_rpm(w0)),pp(spd0(s0f,2)));
fprintf('%12s | %11.1f N.m %11.1f N.m | %10.3f tr/min %8.2f tr/min\n','en charge',...
    pp(oc.Tem(wc)),pp(trq(fc,3)),pp(oc.n_rpm(wc)),pp(spd(scf,2)));
fprintf(['(FEM non vrille = comparaison directe ; la fluctuation de vitesse FEM\n' ...
    ' est exageree par le repli du pas de 1 ms — cf. etude critique §8)\n\n']);
fprintf('Moyennes en charge : n = %.1f (FEM %.1f) tr/min ; T = %.1f (FEM %.1f) N.m\n',...
    mean(oc.n_rpm(wc)),mean(spd(scf,2)),mean(oc.Tem(wc)),mean(trq(fc,3)));

%% 4. Figures (une par condition) : n(t), T(t), zoom, ia(t), T(n)
cases={{'a vide',o0,spd0,trq0,cur0,0.4,0.35},{'en charge',oc,spd,trq,cur,2.0,1.90}};
for j=1:2
    lab=cases{j}{1}; o=cases{j}{2}; sp=cases{j}{3}; tq=cases{j}{4};
    cu=cases{j}{5}; tend=cases{j}{6}; tz=cases{j}{7};
    figure('Name',['Transitoire mono-tranche — ',lab],'Position',[50 30 1150 760]);
    subplot(2,2,1);
    plot(o.t,o.n_rpm,'b-','LineWidth',1.2); hold on;
    plot(sp(:,1),sp(:,2),'r--','LineWidth',0.9);
    grid on; xlim([0 tend]); xlabel('t [s]'); ylabel('n [tr/min]');
    legend('MEC 1 tranche','FEM','Location','southeast');
    title(['(a) Vitesse — ',lab]);
    subplot(2,2,2);
    plot(o.t,o.Tem,'b-','LineWidth',0.7); hold on;
    plot(tq(:,1),tq(:,3),'r--','LineWidth',0.7);
    grid on; xlim([0 tend]); xlabel('t [s]'); ylabel('T_{em} [N·m]');
    legend('MEC 1 tranche','FEM'); title(['(b) Couple — ',lab]);
    subplot(2,2,3);
    mM=o.t>=tz; mF=tq(:,1)>=tz & tq(:,1)<=tz+0.05;
    plot((tq(mF,1)-tz)*1e3,tq(mF,3),'r-','LineWidth',1.0); hold on;
    plot((o.t(mM)-tz)*1e3,o.Tem(mM),'b-','LineWidth',1.0);
    plot((o.t(mM)-tz)*1e3,o.Tem1(mM),'k:','LineWidth',1.0);
    grid on; xlim([0 30]); xlabel('t [ms]'); ylabel('T_{em} [N·m]');
    legend('FEM','MEC 1 tranche','MEC dq seul (sans ondulation)','Location','best');
    title('(c) Zoom regime etabli : ONDULATION de denture');
    subplot(2,2,4);
    plot(o.n_rpm,o.Tem,'b-','LineWidth',0.6); hold on;
    nF=interp1(sp(:,1),sp(:,2),tq(:,1),'linear','extrap');
    plot(nF,tq(:,3),'r--','LineWidth',0.6);
    grid on; xlabel('n [tr/min]'); ylabel('T_{em} [N·m]');
    legend('MEC 1 tranche','FEM','Location','best');
    title('(d) Caracteristique couple–vitesse (trajectoire)');
    fn=sprintf('transitoire_1tranche_%s.png',strrep(lab,' ','_'));
    saveas(gcf,fullfile(OUT,fn));
    fprintf('Figure enregistree : %s\n',fn);
end
fprintf('\n=== Termine ===\n');
