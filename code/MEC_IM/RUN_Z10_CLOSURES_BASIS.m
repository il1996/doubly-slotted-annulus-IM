%% RUN_Z10_CLOSURES_BASIS - les fermetures de RUN_Z1_CAVITY (Carter, Phi_O = 0,
%  cavites) dans les DEUX bases de surface, chapeau symetrique p1 et chapeau
%  asymetrique p1a, aux pavages (17,4), (33,16), (33,32), N_h = 8192,
%  configuration de RUN_Z1_CAVITY (vrillage harmonique neutralise, X_m0 a
%  0,2 A, a vide s = 1e-4, nominal s = 0,0188, references EF relues des .tab,
%  fenetre t >= 1 s). Une seule execution. Le gain 'Phi_O = 0 -> cavites' est
%  imprime dans chaque base separement.
%
%  Comme RUN_Z9_BASIS_P1A : le paquet de production +mec/airgap_fourier.m ne
%  connait que 'p0'/'p1' et n'est PAS modifie ; la base 'p1a' vient de
%  variant_p1a/+mec/airgap_fourier.m (copie de la production + cas 'p1a',
%  MES_R3_VALIDATE), rendue prioritaire comme DOSSIER COURANT (le dossier
%  courant prime sur le path ; addpath seul ne suffit pas). Le marqueur
%  AF.variant_p1a est imprime. La fermeture de Carter (mec.airgap_permeance)
%  ne depend d'aucune base et n'est calculee qu'une fois.
%  Garde 1 : tous les points p1 et Carter contre RUN_Z10_CLOSURES_BASIS_PROD
%  (production seule, session separee). Garde 2 : Carter, Phi_O=0 (17,4),
%  cavites (17,4) et (33,16) contre Z1_cavity.mat (Table 8).
%  Lancement (deux sessions) :
%     matlab -batch "RUN_Z10_CLOSURES_BASIS_PROD"
%     matlab -batch "RUN_Z10_CLOSURES_BASIS"
clear; clc; t0=tic;
here=fileparts(mfilename('fullpath')); if isempty(here), here=pwd; end
vp=fullfile(here,'variant_p1a'); addpath(here); cd(vp);
OUTF=fullfile(here,'Z10_closures_basis_out.txt'); if isfile(OUTF), delete(OUTF); end
diary(OUTF); diary on;
Nh=8192; s_ch=0.0188; Im_ref=0.2;
ROOT=fullfile(here,'..','..','reference','ANSYS_18_5kW');
rd=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
rmsw=@(A,c,tw)sqrt(mean(A(A(:,1)>=tw,c).^2));
avgw=@(A,c,tw)mean(A(A(:,1)>=tw,c));
fprintf('=== Z10 : fermetures de Z1 dans les deux bases p1 / p1a, N_h = %d, s = %.4f ===\n',Nh,s_ch);
fprintf('  airgap_fourier employe : %s\n\n',which('mec.airgap_fourier'));
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
%% ---- les treize points : Carter, puis Phi_O=0 et cavites x 3 pavages x 2 bases --
CAS={'Carter',0,0,'-',false};
for cl={'Phi_O=0','cavites'}
    for pv=[17 4; 33 16; 33 32].'
        for b={'p1','p1a'}, CAS(end+1,:)={cl{1},pv(1),pv(2),b{1},strcmp(cl{1},'cavites')}; end %#ok<SAGROW>
    end
end
nc=size(CAS,1); V=nan(nc,numel(NAMES)); MK=false(nc,1);
for k=1:nc
    M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
    ctx=mec.build_context(M); G=ctx.G; ctx.M.opt.skew_harm=0;
    if strcmp(CAS{k,1},'Carter'), ctx.AG=ctx.AGcarter; MK(k)=true;
    elseif CAS{k,5}, cav=load(fullfile(here,sprintf('cavity_nO%d.mat',CAS{k,3}))); ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,CAS{k,2},CAS{k,3},Nh,CAS{k,4},cav); MK(k)=isfield(ctx.AG.AF,'variant_p1a')&&ctx.AG.AF.variant_p1a;
    else, ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,CAS{k,2},CAS{k,3},Nh,CAS{k,4},[]); MK(k)=isfield(ctx.AG.AF,'variant_p1a')&&ctx.AG.AF.variant_p1a; end
    Rm=mec.magnetizing(ctx,Im_ref); ctx.Xm0=Rm.Xm;
    r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0); rn=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
    V(k,:)=[ctx.Xm0, r0.I1, r0.E1, rn.Xm, rn.I1, rn.Tem, rn.Pfe, rn.Pcu_r, rn.cosphi, rn.eta];
    fprintf('  %-8s (%2d,%2d) %-3s fait (%.0f s) : Xm0 %.3f | I0 %.3f | T %.2f  [variante %d]\n',CAS{k,1},CAS{k,2},CAS{k,3},CAS{k,4},toc(t0),V(k,1),V(k,2),V(k,6),MK(k));
end
if ~all(MK), fprintf('  ATTENTION : la variante n''a pas ete employee pour tous les points\n'); end
%% ---- garde 1 : p1 et Carter contre la production seule ---------------------------
if isfile(fullfile(here,'Z10_closures_basis_prod.mat'))
    P=load(fullfile(here,'Z10_closures_basis_prod.mat')); d1=0;
    for k=1:nc
        if strcmp(CAS{k,4},'p1a'), continue; end
        j=find(strcmp(P.CAS(:,1),CAS{k,1})&cell2mat(P.CAS(:,2))==CAS{k,2}&cell2mat(P.CAS(:,3))==CAS{k,3});
        d1=max(d1,max(abs(V(k,:)-P.V(j,:))./abs(P.V(j,:))));
    end
    fprintf('  garde 1 : Carter et points p1, variante contre production seule (Z10_closures_basis_prod.mat, variante sur le path : %s) : ecart relatif max %.1e\n',mat2str(double(P.MK.')),d1);
else
    fprintf('  garde 1 : Z10_closures_basis_prod.mat absent (lancer RUN_Z10_CLOSURES_BASIS_PROD d''abord)\n');
end
%% ---- garde 2 : contre Z1_cavity.mat (Table 8) ----------------------------------------
zf=fullfile(here,'..','..','outputs','MEC_IM','Z1_cavity.mat');
if isfile(zf)
    Z=load(zf); pairs={'Carter',0,0,1; 'Phi_O=0',17,4,2; 'cavites',17,4,3; 'cavites',33,16,5}; d2=0;
    for q=1:size(pairs,1)
        k=find(strcmp(CAS(:,1),pairs{q,1})&cell2mat(CAS(:,2))==pairs{q,2}&cell2mat(CAS(:,3))==pairs{q,3}&~strcmp(CAS(:,4),'p1a'));
        z=Z.V(pairs{q,4},:); d2=max(d2,max(abs(V(k,1:8)-z([1 3 4 6 7 9 12 13]))./abs(z([1 3 4 6 7 9 12 13]))));
    end
    fprintf('  garde 2 : Carter, Phi_O=0 (17,4), cavites (17,4) et (33,16) en p1 contre Z1_cavity.mat (Table 8) : ecart relatif max %.1e\n',d2);
end
%% ---- tableau complet -------------------------------------------------------------------
fprintf('\n  %-8s %7s %4s | %9s %8s %8s %9s %8s %8s %7s %7s %7s %7s\n','fermeture','pavage','base','Xm0','I0','E1(0)','Xm(s)','I1','T','Pfe','Pcu_r','cosphi','eta');
for k=1:nc
    fprintf('  %-8s (%2d,%2d) %-4s | %9.3f %8.3f %8.2f %9.3f %8.3f %8.2f %7.1f %7.1f %7.4f %7.4f\n',CAS{k,1},CAS{k,2},CAS{k,3},CAS{k,4},V(k,:));
end
fprintf('  %-8s %7s %4s | %9s %8.3f %8.2f %9s %8.3f %8.2f %7.1f %7.1f\n','ref. EF','','','---',I0F,E0F,'---',I1F,TF,PfeF,PrF);
fprintf('\n  ecarts a la reference EF (%%) :\n  %-8s %7s %4s | %8s %8s %8s %8s %8s %8s\n','fermeture','pavage','base','I0','E1(0)','I1','T','Pfe','Pcu_r');
E=100*(V(:,[2 3 5 6 7 8])-REF([2 3 5 6 7 8]))./REF([2 3 5 6 7 8]);
for k=1:nc, fprintf('  %-8s (%2d,%2d) %-4s | %+8.2f %+8.2f %+8.2f %+8.2f %+8.2f %+8.2f\n',CAS{k,1},CAS{k,2},CAS{k,3},CAS{k,4},E(k,:)); end
%% ---- ecart entre bases, par fermeture et pavage (points d'ecart a l'EF ; % pour Xm) ----
fprintf('\n  ecart p1a - p1 (points d''ecart a l''EF ; pour Xm0, Xm(s), cosphi, eta : %% p1a/p1 - 1) :\n');
fprintf('  %-8s %7s | %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n','fermeture','pavage','Xm0','I0','E1(0)','Xm(s)','I1','T','Pfe','Pcu_r','cosphi','eta');
for k=2:2:nc
    ka=k+1; row=nan(1,10);
    for j=1:10, if isnan(REF(j)), row(j)=100*(V(ka,j)-V(k,j))/V(k,j); else, row(j)=100*(V(ka,j)-REF(j))/REF(j)-100*(V(k,j)-REF(j))/REF(j); end, end
    fprintf('  %-8s (%2d,%2d) | %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f\n',CAS{k,1},CAS{k,2},CAS{k,3},row);
end
%% ---- gain Phi_O = 0 -> cavites, par base ---------------------------------------------------
fprintf('\n  gain Phi_O=0 -> cavites (points d''ecart a l''EF : ecart(cavites) - ecart(Phi_O=0)) :\n');
fprintf('  %-28s %4s | %8s %8s %8s %8s %8s %8s\n','paire de pavages','base','I0','E1(0)','I1','T','Pfe','Pcu_r');
idx=@(cl,nT,nO,b) find(strcmp(CAS(:,1),cl)&cell2mat(CAS(:,2))==nT&cell2mat(CAS(:,3))==nO&strcmp(CAS(:,4),b));
PAIRS={[17 4],[33 16],'Phi_O=0 (17,4) -> cav (33,16)'; [17 4],[17 4],'Phi_O=0 (17,4) -> cav (17,4)'; [33 16],[33 16],'Phi_O=0 (33,16) -> cav (33,16)'; [33 32],[33 32],'Phi_O=0 (33,32) -> cav (33,32)'};
for q=1:size(PAIRS,1)
    for b={'p1','p1a'}
        i1=idx('Phi_O=0',PAIRS{q,1}(1),PAIRS{q,1}(2),b{1}); i2=idx('cavites',PAIRS{q,2}(1),PAIRS{q,2}(2),b{1});
        fprintf('  %-28s %-4s | %+8.2f %+8.2f %+8.2f %+8.2f %+8.2f %+8.2f\n',PAIRS{q,3},b{1},E(i2,:)-E(i1,:));
    end
end
fprintf('  (la paire du resume, "+21.5 %% to +7.7 %%", est Phi_O=0 (17,4) -> cavites (33,16) sur I0)\n');
save(fullfile(here,'Z10_closures_basis.mat'),'V','CAS','NAMES','REF','E','MK'); cd(here);
fprintf('\n  duree %.0f s\n=== Z10 termine ===\n',toc(t0)); diary off;
