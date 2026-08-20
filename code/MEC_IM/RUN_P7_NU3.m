%% RUN_P7_NU3  -  P7 : harmonique de saturation nu=3 — attribution PAR LES DONNEES
%
%  QUESTION. Les barres ANSYS portent 27,3 A rms a VIDE. Trois mecanismes
%  candidats, aux SIGNATURES FREQUENTIELLES distinctes dans le repere rotor :
%    * denture stator nu=23 (directe inverse) : f_rot = |1-23(1-s)|*f
%         a vide ~22f = 1100 Hz ;  en charge 21,57f = 1078 Hz
%    * denture stator nu=25 (inverse)         : f_rot = (1+25(1-s))*f
%         a vide 26f = 1300 Hz ;   en charge 26,53f = 1327 Hz
%    * harmonique de SATURATION nu=3 (co-tournant au synchronisme) :
%         f_rot = 3*s*f  ->  ~0 Hz a vide, 2,8 Hz en charge.
%  Les .tab ANSYS sont echantillonnes a fs = 1 kHz (Nyquist 500 Hz) : les
%  raies de denture se REPLIENT a |f-1000| = 100/300 Hz (vide) et
%  78/327 Hz (charge) — bandes ou rien d'autre n'existe : discriminant.
%
%  Ce script : (1) spectre des courants de barre/anneau ANSYS (vide+charge) ;
%  (2) amplitude B3 du profil d'entrefer (MEC Fourier vs ANSYS) — la
%  saturation aplatit bien l'onde ; (3) verdict d'attribution.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
fprintf('=== P7 : attribution des courants de cage a vide (nu=3 ?) ===\n\n');

M=mec.machine_18_5kW(); ctx=mec.build_context(M);
G=ctx.G; W=ctx.W; BH=ctx.BH; p=M.p; Nr=M.Nr;

%% 1. Spectres des courants de cage ANSYS
cases={{'a vide',1e-4},{'en charge',0.0188}};
for cc=1:2
    lab=cases{cc}{1}; s=cases{cc}{2};
    ec=readmatrix(fullfile(ROOT,'transitoire',lab,'End Connection Plot 1.tab'),...
        'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    m=ec(:,1)>=1.0; t=ec(m,1); x=1e3*ec(m,2);          % I_barre [A]
    fs=1/mean(diff(t)); N=numel(x);
    w=0.5*(1-cos(2*pi*(0:N-1)'/(N-1)));               % fenetre de Hann
    X=fft((x-mean(x)).*w); f=(0:N-1)*fs/N;
    A=2*abs(X)/sum(w); n2=floor(N/2);
    Ah=A(1:n2); loc=find(Ah(2:end-1)>Ah(1:end-2) & Ah(2:end-1)>Ah(3:end))+1;
    [pk,o]=sort(Ah(loc),'descend'); ip=loc(o);
    nk=min(6,numel(pk)); pk=pk(1:nk); ip=ip(1:nk);
    fprintf('--- I_barre %s (fs=%.0f Hz, RMS=%.1f A) : raies dominantes ---\n',...
        lab,fs,rms(x));
    % sens 6k±1 : nu=23 (6k-1) INVERSE -> (1+23(1-s))f ; nu=25 (6k+1)
    % DIRECT -> |1-25(1-s)|f : les DEUX tombent a ~24f (1200 Hz), comme les
    % bandes de permeance Ns±p -> TOUTE la denture stator se replie ensemble.
    fA=(1+23*(1-s))*M.f;                     % ~24f
    f22=abs(1-23*(1-s))*M.f; f26=(1+25*(1-s))*M.f;   % raies secondaires
    fprintf(['attendu (replie) : denture nu=23/25+permeance -> %.1f Hz ; ',...
        'secondaires -> %.1f / %.1f Hz ; nu=3 -> %.2f Hz\n'],...
        abs(fA-fs*round(fA/fs)),abs(f22-fs*round(f22/fs)),...
        abs(f26-fs*round(f26/fs)),3*s*M.f);
    for q=1:numel(pk)
        fprintf('   %7.1f Hz : %7.2f A crete\n',f(ip(q)),pk(q));
    end
    fprintf('\n');
end

%% 2. Aplatissement de saturation : B3 du profil d'entrefer
% ANSYS (profil mi-entrefer, snapshot)
fu=@(y,th,k)abs((2/numel(y))*sum(y(:).*exp(-1i*k*th(:))));
Pr=readmatrix(fullfile(ROOT,'transitoire','a vide','Calculator Expressions Plot 1.tab'),...
    'FileType','text','NumHeaderLines',1,'Delimiter','\t');
n=size(Pr,1)-1; thA=2*pi*Pr(1:n,3)/Pr(end,3); BrA=Pr(1:n,6);
% MEC (entrefer Fourier, a vide)
Mf=M; Mf.opt.gap_fourier=1;
r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
psi1=angle(r0.I1c);
i30=sqrt(2)*r0.I1*[cos(psi1);cos(psi1-2*pi/3);cos(psi1+2*pi/3)];
me=mec.mesh_refined(Mf,G,6,i30,0,3,zeros(Nr,1));
Se=mec.solve_mesh(me,BH,Mf.opt);
Us=Se.U(me.gapF.ids(1:me.Ms)); Ur=Se.U(me.gapF.ids(me.Ms+1:end));
thq=linspace(0,2*pi,2001); thq(end)=[];
Br=me.gapF.field(Us,Ur,0.5*(G.Rs+G.Rr),thq);
fprintf('--- 2. Harmonique de saturation B3 (mi-entrefer, a vide) ---\n');
fprintf('%12s | %8s %8s %10s\n','','B1 [T]','B3 [T]','B3/B1');
fprintf('%12s | %8.3f %8.3f %9.1f %%\n','EF (ANSYS)',fu(BrA,thA,p),fu(BrA,thA,3*p),...
    100*fu(BrA,thA,3*p)/fu(BrA,thA,p));
fprintf('%12s | %8.3f %8.3f %9.1f %%\n','MEC Fourier',fu(Br,thq,p),fu(Br,thq,3*p),...
    100*fu(Br,thq,3*p)/fu(Br,thq,p));

%% 3. Verdict (etabli sur les donnees du 18/07/2026)
fprintf('\n--- 3. Verdict (mesure) ---\n');
fprintf('1) A VIDE, la raie dominante est ~200 Hz = repli de 24f (36,6 A crete\n');
fprintf('   ~ la quasi-totalite du RMS) : c''est la DENTURE STATOR complete\n');
fprintf('   (nu=23 inverse + nu=25 direct + bandes de permeance -> toutes a 24f).\n');
fprintf('   nu=3 (%.2f Hz) : rien. => Le modele EF est effectivement NON VRILLE\n',3*1e-4*M.f);
fprintf('   pour la cage : l''ecart I_barre a vide post-C5 (7,6 vs 27,3 A) est un\n');
fprintf('   ecart de REFERENCE (machine vrillee vs EF non vrille), pas une\n');
fprintf('   physique manquante.\n');
fprintf('2) EN CHARGE : fondamental de glissement s*f (~1 Hz, ~308 A rms,\n');
fprintf('   coherent avec le modele), denture repliee ~178 Hz, nu=5/7 ~295 Hz.\n');
fprintf('3) B3 : EF = 0,8 %% seulement — en alimentation en TENSION, le courant\n');
fprintf('   magnetisant pointu maintient le flux quasi sinusoidal ; le 11,6 %%\n');
fprintf('   du MEC est l''artefact d''une excitation sinusoidale en COURANT du\n');
fprintf('   banc de champ. => L''hypothese « nu=3 manquant » est REFUTEE.\n');
fprintf('=== Termine ===\n');
