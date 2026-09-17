%% RUN_Z9_BASIS_P1A_PROD - premiere moitie de Z9 : les points en base p1 avec le
%  paquet de production SEUL (aucune variante sur le path), sauvegardes dans
%  Z9_basis_p1a_prod.mat pour la garde 1 de RUN_Z9_BASIS_P1A. Meme
%  configuration que RUN_Z1_CAVITY (vrillage harmonique neutralise), fermeture
%  cavites, N_h = 8192 ; cavity_nO32.mat de export_cavity_nO32.py.
clear; clc; t0=tic;
Nh=8192; s_ch=0.0188; Im_ref=0.2;
here=fileparts(mfilename('fullpath')); if isempty(here), here=pwd; end
vp=fullfile(here,'variant_p1a'); if contains(path,vp), error('Z9_PROD : la variante p1a est sur le path ; ce script doit tourner avec la production seule'); end
CAS={33,16,'p1'; 33,32,'p1'}; V=nan(2,10); MK=false(2,1);
for k=1:2
    M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
    ctx=mec.build_context(M); G=ctx.G; ctx.M.opt.skew_harm=0;
    cav=load(sprintf('cavity_nO%d.mat',CAS{k,2}));
    ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,CAS{k,1},CAS{k,2},Nh,CAS{k,3},cav); MK(k)=isfield(ctx.AG.AF,'variant_p1a');
    Rm=mec.magnetizing(ctx,Im_ref); ctx.Xm0=Rm.Xm;
    r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0); rn=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
    V(k,:)=[ctx.Xm0, r0.I1, r0.E1, rn.Xm, rn.I1, rn.Tem, rn.Pfe, rn.Pcu_r, rn.cosphi, rn.eta];
    fprintf('  Z9_PROD (%d,%d) %s : Xm0 %.3f | I0 %.3f | T %.2f  [variante %d]  (%.0f s)\n',CAS{k,1},CAS{k,2},CAS{k,3},V(k,1),V(k,2),V(k,6),MK(k),toc(t0));
end
save('Z9_basis_p1a_prod.mat','V','CAS','MK');
