%% RUN_Z9_BASIS_P1A - sensibilite de base au niveau machine : chapeau symetrique
%  (p1) contre chapeau asymetrique (p1a), fermeture cavites, pavages (33,16)
%  et (33,32), N_h = 8192, configuration loyale de RUN_Z1_CAVITY (vrillage
%  harmonique neutralise, references EF relues des .tab, fenetre t >= 1 s).
%
%  Le paquet de production +mec/airgap_fourier.m ne connait que 'p0' et 'p1'.
%  Il n'est PAS modifie : la base 'p1a' vient de variant_p1a/+mec/airgap_fourier.m,
%  copie de la production augmentee du cas 'p1a' (identique a la production pour
%  'p0'/'p1', identique a dsop.assemble pour 'p1a' : MES_R3_VALIDATE, 17/09/2026),
%  rendue prioritaire en en faisant le DOSSIER COURANT (le dossier courant prime
%  sur le path, et le +mec de production, dans le dossier de ce script, prime sur
%  toute entree de path tant qu'il est le dossier courant -- addpath seul ne
%  suffit donc pas). Le marqueur AF.variant_p1a est imprime pour chaque point.
%  Garde 1 : les points p1 sont recalcules ici avec la variante et compares aux
%  memes points resolus par RUN_Z9_BASIS_P1A_PROD avec la production seule
%  (Z9_basis_p1a_prod.mat, session MATLAB separee) : ils doivent coincider.
%  Garde 2 : (33,16) p1 doit reproduire la ligne 'operateur+cavite (33,16)' de
%  Z1_cavity.mat / Z1_cavity_out.txt.
%  Les matrices de cavite nO = 32 viennent de cavity_nO32.mat
%  (export_cavity_nO32.py, meme construction que export_cavity.py, lc_mouth 0.005).
%  Lancement (deux sessions) :
%     matlab -batch "RUN_Z9_BASIS_P1A_PROD"
%     matlab -batch "RUN_Z9_BASIS_P1A"
clear; clc; t0=tic;
here=fileparts(mfilename('fullpath')); if isempty(here), here=pwd; end
vp=fullfile(here,'variant_p1a'); addpath(here); cd(vp);          % variante = dossier courant, production sur le path
OUTF=fullfile(here,'Z9_basis_p1a_out.txt'); if isfile(OUTF), delete(OUTF); end
diary(OUTF); diary on;
Nh=8192; s_ch=0.0188; Im_ref=0.2;
ROOT=fullfile(here,'..','..','reference','ANSYS_18_5kW');
rd=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
rmsw=@(A,c,tw)sqrt(mean(A(A(:,1)>=tw,c).^2));
avgw=@(A,c,tw)mean(A(A(:,1)>=tw,c));
fprintf('=== Z9 : sensibilite de base p1 / p1a, fermeture cavites, N_h = %d, s = %.4f ===\n',Nh,s_ch);
fprintf('  variante p1a : %s (en tete du path)\n\n',which('mec.airgap_fourier'));
%% ---- references EF, relues (comme Z1) --------------------------------------
D0=fullfile(ROOT,'transitoire','a vide'); DC=fullfile(ROOT,'transitoire','en charge'); tw=1.0;
A=rd(fullfile(D0,'Winding Plot 4.tab')); I0F=mean([rmsw(A,2,tw) rmsw(A,3,tw) rmsw(A,4,tw)]);
A=rd(fullfile(D0,'Winding Plot 2.tab')); E0F=mean([rmsw(A,2,tw) rmsw(A,3,tw) rmsw(A,4,tw)]);
A=rd(fullfile(DC,'Plot 1.tab'));         TF=avgw(A,3,tw);
A=rd(fullfile(DC,'Winding Plot 4.tab')); I1F=mean([rmsw(A,2,tw) rmsw(A,3,tw) rmsw(A,4,tw)]);
A=rd(fullfile(DC,'Loss Plot 1.tab'));    PfeF=avgw(A,2,tw);
A=rd(fullfile(DC,'Loss Plot 2.tab'));    PbarF=1e3*avgw(A,2,tw);
A=rd(fullfile(DC,'End Connection Plot 3.tab')); PringF=1e3*avgw(A,2,tw);
PrF=PbarF+PringF;
fprintf('  references EF (t >= %.1f s) : I0 %.4f A | E(2D) a vide %.2f V | T %.3f N.m | I1 %.4f A | Pfe %.1f W | Pcu_r (barres+anneaux) %.1f W\n\n',tw,I0F,E0F,TF,I1F,PfeF,PrF);
NAMES={'Xm0 (0,2 A) [ohm]','I0 [A]','E1(0) [V]','Xm(s) [ohm]','I1 [A]','T [N.m]','Pfe [W]','Pcu_r [W]','cosphi','eta'};
REF=[NaN, I0F, E0F, NaN, I1F, TF, PfeF, PrF, NaN, NaN];
%% ---- les quatre points ---------------------------------------------------------
CAS={33,16,'p1'; 33,16,'p1a'; 33,32,'p1'; 33,32,'p1a'}; nc=size(CAS,1); V=nan(nc,numel(NAMES)); MK=false(nc,1);
for k=1:nc
    M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
    ctx=mec.build_context(M); G=ctx.G; ctx.M.opt.skew_harm=0;
    cav=load(fullfile(here,sprintf('cavity_nO%d.mat',CAS{k,2})));
    ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,CAS{k,1},CAS{k,2},Nh,CAS{k,3},cav); MK(k)=isfield(ctx.AG.AF,'variant_p1a')&&ctx.AG.AF.variant_p1a;
    Rm=mec.magnetizing(ctx,Im_ref); ctx.Xm0=Rm.Xm;
    r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0); rn=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
    V(k,:)=[ctx.Xm0, r0.I1, r0.E1, rn.Xm, rn.I1, rn.Tem, rn.Pfe, rn.Pcu_r, rn.cosphi, rn.eta];
    fprintf('  (%d,%d) %-3s fait (%.0f s) : Xm0 %.3f | I0 %.3f | T %.2f  [variante %d]\n',CAS{k,1},CAS{k,2},CAS{k,3},toc(t0),V(k,1),V(k,2),V(k,6),MK(k));
end
if ~all(MK), fprintf('  ATTENTION : la variante n''a pas ete employee pour tous les points\n'); end
%% ---- garde 1 : les points p1 contre la production seule (session separee) -----
if isfile(fullfile(here,'Z9_basis_p1a_prod.mat'))
    P=load(fullfile(here,'Z9_basis_p1a_prod.mat'));
    d1=max(abs(V([1 3],:)-P.V)./abs(P.V),[],'all');
    fprintf('  garde 1 : points p1 (33,16) et (33,32), variante contre production seule (Z9_basis_p1a_prod.mat, variante sur le path : %d %d) : ecart relatif max %.1e\n',P.MK(1),P.MK(2),d1);
else
    fprintf('  garde 1 : Z9_basis_p1a_prod.mat absent (lancer RUN_Z9_BASIS_P1A_PROD d''abord)\n');
end
%% ---- garde 2 : (33,16) p1 contre Z1 ---------------------------------------------
if isfile(fullfile(here,'..','..','outputs','MEC_IM','Z1_cavity.mat'))
    Z=load(fullfile(here,'..','..','outputs','MEC_IM','Z1_cavity.mat')); z=Z.V(5,:);   % ligne 'operateur+cavite (33,16)' : [Xm0 Bg1f I0 E1 Bg1 Xm I1 E1s T cosphi eta Pfe Pcu_r Pcu_s]
    d2=max(abs(V(1,1:8)-z([1 3 4 6 7 9 12 13]))./abs(z([1 3 4 6 7 9 12 13])));
    fprintf('  garde 2 : (33,16) p1 contre Z1_cavity.mat, ligne (33,16) : ecart relatif max %.1e\n',d2);
else
    fprintf('  garde 2 : Z1_cavity.mat absent ; Z1_cavity_out.txt donne Xm0 70.029 | I0 9.153 | T 115.53\n');
end
%% ---- tableau ------------------------------------------------------------------------
fprintf('\n  %-20s %12s %12s %12s | %12s %12s %12s | %10s\n','grandeur','(33,16) p1','(33,16) p1a','ecart','(33,32) p1','(33,32) p1a','ecart','ref. EF');
for j=1:numel(NAMES)
    if isnan(REF(j)), e1=100*(V(2,j)-V(1,j))/V(1,j); e2=100*(V(4,j)-V(3,j))/V(3,j); lab='(% p1a/p1 - 1)';
    else, e1=100*(V(2,j)-REF(j))/REF(j)-100*(V(1,j)-REF(j))/REF(j); e2=100*(V(4,j)-REF(j))/REF(j)-100*(V(3,j)-REF(j))/REF(j); lab='(points, ecart a l''EF)'; end
    if isnan(REF(j)), rs='---'; else, rs=sprintf('%.4f',REF(j)); end
    fprintf('  %-20s %12.4f %12.4f %+12.3f | %12.4f %12.4f %+12.3f | %10s  %s\n',NAMES{j},V(1,j),V(2,j),e1,V(3,j),V(4,j),e2,rs,lab);
end
fprintf('\n  ecarts a la reference EF (%%) :\n  %-20s %12s %12s %12s %12s\n','grandeur','(33,16) p1','(33,16) p1a','(33,32) p1','(33,32) p1a');
for j=find(~isnan(REF))
    fprintf('  %-20s %+12.2f %+12.2f %+12.2f %+12.2f\n',NAMES{j},100*(V(:,j)-REF(j))/REF(j));
end
save(fullfile(here,'Z9_basis_p1a.mat'),'V','CAS','NAMES','REF','MK'); cd(here);
fprintf('\n  duree %.0f s\n=== Z9 termine ===\n',toc(t0)); diary off;
