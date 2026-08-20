%% RUN_VERIFY_OPERATOR - tests d'acceptation de l'operateur DtN (MAS)
%
%  Cahier des charges SPEC_REFACTOR_OPERATEUR_MAS.md, section 6.
%  T4 est passe EN PREMIER : c'est celui qui attrape une erreur de signe,
%  laquelle converge vers une solution fausse sans rien signaler.
clear; clc;
M=mec.machine_18_5kW(); ctx=mec.build_context(M);
G=ctx.G; AT=ctx.AG; AF=AT.AF; g=ctx.gap; mu0=4*pi*1e-7;
Rs=G.Rs; Rr=G.Rr; L=M.L; Rm=0.5*(Rs+Rr); gg=Rs-Rr;
thsT=2*pi*(0:M.Ns-1)/M.Ns; thrT=2*pi*(0:M.Nr-1)/M.Nr;
ok=@(b)char(9989*b+10060*(~b));
fprintf('=== Tests d''acceptation de l''operateur DtN ===\n');
fprintf('  Rs=%.4f m  Rr=%.4f m  g=%.4f mm  X=ln(Rs/Rr)=%.3e\n', ...
    Rs,Rr,gg*1e3,AT.X);
fprintf('  grille fine %d+%d colonnes -> operateur condense %dx%d, Nh=%d\n\n', ...
    AT.Msf,AT.Mrf,size(AT.Y,1),size(AT.Y,2),AT.Nh);

%% ---- T4 : limites analytiques des noyaux ------------------------------
n=AF.n; Kn=mu0*pi*L*n;
lim1=mu0*pi*Rm*L/(2*gg);                 % un seul cote alimente
e1=abs((Kn.*AF.Cn/2)/lim1-1);
lim2=mu0*pi*(n.^2)*gg*L/(2*Rm);          % memes potentiels des 2 cotes
e2=abs((Kn.*(AF.Cn-AF.Dn))./lim2-1);
fprintf('T4 limites analytiques (erreur relative, harmoniques 1..10)\n');
fprintf('   un cote alimente   -> mu0*pi*R*L/(2g)   : max %.2e\n',max(e1(1:10)));
fprintf('   deux cotes egaux   -> mu0*pi*n^2*g*L/2R : max %.2e   %s\n\n', ...
    max(e2(1:10)),ok(max(e2(1:10))<2e-3));

%% ---- T2 : symetrie ----------------------------------------------------
%  La condensation de Schur doit PRESERVER la symetrie : Y_c = P'(Ytt -
%  Yto*Yoo^-1*Yot)P est symetrique des que Y l'est. On teste les deux.
for c={'fin (grille de surface)',AF.Y; 'condense (noeuds dentaires)',AT.Y}.'
    Y=c{2}; s2=norm(Y-Y.','fro')/norm(Y,'fro');
    ev=sort(eig(full((Y+Y.')/2)));
    fprintf('T2 %-28s ||Y-Y''||/||Y|| = %.2e   %s\n',c{1},s2,ok(s2<1e-12));
    fprintf('   vp min %10.3e  max %10.3e   (semi-def. pos. : le mode nul\n',ev(1),ev(end));
    fprintf('   est le mode homopolaire, U = cste ne debite aucun flux)\n');
end
fprintf('\n');

%% ---- T3 : troncature ---------------------------------------------------
%  Nh doit resoudre la PLUS FINE colonne de la grille, pas seulement le
%  couplage 1/X. Le critere automatique doit se voir converger.
fprintf('T3 troncature : Xm selon Nh (defaut automatique = %d)\n',AT.Nh);
fprintf('   %6s %8s %12s %12s\n','Nh','rang','Xm [ohm]','Bg1 [T]');
for Nh=[192 384 768 1536 AT.Nh]
    A=mec.airgap_dtn_tooth(M,G,0,g.nT,g.nO,Nh);
    c2=ctx; c2.AG=A; R=mec.magnetizing(c2,0.2);
    fprintf('   %6d %8d %12.2f %12.5f\n',Nh,rank(A.Y),R.Xm,R.Bg1);
end

%% ---- CONTROLE PHYSIQUE : entrefer LISSE -------------------------------
%  Si l'operateur et son assemblage sont justes, une couronne a surfaces
%  LISSES (arcs pavant exactement l'alesage) doit redonner la permeance
%  d'entrefer lisse analytique, et Xm doit alors DEPASSER la valeur
%  encochee (pas d'effet Carter). Le sens de cet ecart est le test.
fprintf('\nCONTROLE entrefer lisse (arcs pavant 100 %% de l''alesage)\n');
taus=2*pi/M.Ns; taur=2*pi/M.Nr;
AFs=mec.airgap_fourier(thsT,repmat(taus,1,M.Ns),thrT,repmat(taur,1,M.Nr), ...
                       Rs,Rr,L,AT.Nh);
cs=ctx; cs.AG=AFs; Rl=mec.magnetizing(cs,0.2);
kCimp=Rl.Xm/ctx.Xm0;
fprintf('   Xm surfaces LISSES   = %.2f ohm\n',Rl.Xm);
fprintf('   Xm operateur condense= %.2f ohm\n',ctx.Xm0);
fprintf('   => Carter IMPLICITE  = %.3f   [Carter classique %.3f]\n', ...
    kCimp,ctx.AGcarter.kC);
fprintf('   attendu : Xm(lisse) > Xm(encoche)   %s\n',ok(Rl.Xm>ctx.Xm0));

%% ---- T1 : invariance au raffinement de la grille de surface ----------
%  DEUX EFFETS A NE PAS CONFONDRE. Raffiner la grille cree des colonnes
%  plus etroites, qui reclament PLUS d'harmoniques : laisser Nh suivre le
%  critere automatique melange l'invariance au maillage et la troncature,
%  et le plafond Nh=4000 fausse alors les deux grilles les plus fines.
%  On mesure donc l'invariance a Nh FIXE (c'est le test), et on donne la
%  colonne Nh automatique a cote (c'est l'usage).
NhF=AT.Nh;
fprintf('\nT1 invariance au raffinement (nT x nO colonnes par pas dentaire)\n');
fprintf('   %9s %8s %10s %10s | %8s %10s\n', ...
    'nT x nO','colonnes','Xm(NhF)','Bg1(NhF)','Nh auto','Xm(auto)');
Xf=[]; Xa=[];
for cfg=[3 1; 6 2; 9 3; 12 4].'
    A=mec.airgap_dtn_tooth(M,G,0,cfg(1),cfg(2),NhF);
    cx=ctx; cx.AG=A; R=mec.magnetizing(cx,0.2);
    B=mec.airgap_dtn_tooth(M,G,0,cfg(1),cfg(2));
    cy=ctx; cy.AG=B; Rb=mec.magnetizing(cy,0.2);
    Xf(end+1)=R.Xm; Xa(end+1)=Rb.Xm; %#ok<SAGROW>
    fprintf('   %5dx%-3d %8d %10.2f %10.5f | %8d %10.2f\n', ...
        cfg(1),cfg(2),A.Msf+A.Mrf,R.Xm,R.Bg1,B.Nh,Rb.Xm);
end
df=100*(max(Xf)-min(Xf))/mean(Xf);
da=100*(max(Xa)-min(Xa))/mean(Xa);
fprintf('   dispersion a Nh=%d fixe : %.2f %%   %s\n',NhF,df,ok(df<1));
fprintf('   dispersion a Nh automatique : %.2f %% (inclut la troncature)\n',da);

%% ---- T5 : bilan de puissance sur 30 glissements ----------------------
s_list=[0.005:0.005:0.12,0.15,0.2,0.3,0.5,0.7,1.0];
Xp=ctx.Xm0; b6=0;
for s=s_list
    r=mec.equivalent_circuit(ctx,s,Xp); Xp=r.Xm;
    c=mec.power_balance(r,M); b6=max(b6,c.err_global);
end
fprintf('\nT5 bilan de puissance B6, %d glissements : %.2e   %s\n', ...
    numel(s_list),b6,ok(b6<=2e-5));

%% ---- T6 : non-regression a vide (tout ecart EST un resultat) ---------
r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
fprintf('\nT6 non-regression a vide (ancien modele Carter -> operateur)\n');
fprintf('   %-10s %10s %10s %10s\n','','Carter','DtN','EF');
fprintf('   %-10s %10.2f %10.2f %10.2f\n','I0 [A]',8.58,r0.I1,8.49);
fprintf('   %-10s %10.1f %10.1f %10.1f\n','E1 [V]',378.5,r0.E1,382.1);
fprintf('   ecart DtN/EF : I0 %+.1f %% , E1 %+.1f %%\n', ...
    100*(r0.I1-8.49)/8.49,100*(r0.E1-382.1)/382.1);

%% ---- T7 : composante tangentielle en charge --------------------------
%  Ce que Carter ne produit pas du tout, l'operateur le produit.
%  EXCITATION COMPLETE : courant statorique TOTAL + onde de FMM de barres,
%  comme mec.torque_vw. Le champ de la seule branche magnetisante (Im) ne
%  porte quasiment pas de composante tangentielle : c'est la reaction
%  d'induit qui la cree.
sc=0.0188; rc=mec.equivalent_circuit(ctx,sc,ctx.Xm0);
ps1=angle(rc.I1c);
i3=sqrt(2)*abs(rc.I1c)*[cos(ps1);cos(ps1-2*pi/3);cos(ps1+2*pi/3)];
Fs=ctx.W.slotMMF(i3);
ps2=angle(rc.I2c); thb=2*pi*(0:M.Nr-1)/M.Nr;
Fu=cumsum(cos(M.p*thb-ps2).'); Fu=Fu-mean(Fu);
cc1=(2/M.Nr)*sum(Fu.'.*exp(-1j*M.p*thb));
Fr=-Fu*((3/2)*(4/pi)*(ctx.W.kw1*ctx.W.Nph/(2*M.p))*sqrt(2)*abs(rc.I2c)/abs(cc1));
Sl=mec.solve_network(ctx.net,G,ctx.BH,AT,Fs,Fr,M.opt);
thq=linspace(0,2*pi,2001); thq(end)=[];
[Br,Bt]=AT.field(Sl.Usurf,Rm,thq);
Br=Br(:).'; Bt=Bt(:).';
Bg1L=abs((2/numel(thq))*sum(Br.*exp(-1j*M.p*thq)));
fprintf('\nT7 champ au mi-entrefer EN CHARGE (s=%.4f), operateur condense\n',sc);
fprintf('   %-14s %10s %10s %9s\n','','MEC/DtN','EF','ecart');
fprintf('   %-14s %10.4f %10.4f %+8.1f %%\n','Bg1 [T]',Bg1L,0.920,100*(Bg1L-0.920)/0.920);
Btr=sqrt(mean(Bt.^2));
fprintf('   %-14s %10.4f %10.4f %+8.1f %%\n','Bt rms [T]',Btr,0.131,100*(Btr-0.131)/0.131);
fprintf('   rappel : le modele de Carter ne produit AUCUN Bt (permeances\n');
fprintf('   purement radiales) — cette ligne n''existait pas avant.\n');

%% ---- T8 : absence de Carter sur le chemin actif ----------------------
fprintf('\nT8 recherche de Carter hors modele de comparaison\n');
d=dir(fullfile(fileparts(mfilename('fullpath')),'+mec','*.m'));
hit=0;
for k=1:numel(d)
    if any(strcmp(d(k).name,{'airgap_permeance.m','carter.m'})), continue; end
    txt=fileread(fullfile(d(k).folder,d(k).name));
    L=regexp(txt,'[^\n]*','match');
    for j=1:numel(L)
        s=L{j}; if ~isempty(regexp(s,'^\s*%','once')), continue; end
        if isempty(regexp(s,'kC|carter|g_eff','once')), continue; end
        %  Seule tolerance : l'instanciation du modele de COMPARAISON,
        %  qui n'est pas sur le chemin de resolution.
        if ~isempty(regexp(s,'AGcarter\s*=','once'))
            fprintf('   (tolere) %s:%d  %s\n',d(k).name,j,strtrim(s));
        else
            fprintf('   %s:%d  %s\n',d(k).name,j,strtrim(s)); hit=hit+1;
        end
    end
end
fprintf('   %d occurrence(s) hors commentaire   %s\n',hit,ok(hit==0));

%% ---- reperes ----------------------------------------------------------
fprintf('\nReperes\n');
Xm_carter=NaN;
try
    cc=ctx; cc.AG=ctx.AGcarter; Rc=mec.magnetizing(cc,0.2); Xm_carter=Rc.Xm;
catch ME
    fprintf('   (modele Carter indisponible : %s)\n',ME.message);
end
fprintf('   Xm  Carter (ancien defaut) : %8.2f ohm\n',Xm_carter);
fprintf('   Xm  operateur DtN          : %8.2f ohm\n',ctx.Xm0);
fprintf('   Xm  reference EF           : %8.2f ohm\n',46);
fprintf('   kC = %.3f  g_eff = %.4f mm  (modele Carter)\n', ...
    ctx.AGcarter.kC,ctx.AGcarter.g_eff*1e3);
