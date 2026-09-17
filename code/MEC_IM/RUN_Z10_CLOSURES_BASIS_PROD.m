%% RUN_Z10_CLOSURES_BASIS_PROD - premiere moitie de Z10 : toutes les fermetures de
%  RUN_Z1_CAVITY (Carter, Phi_O = 0, cavites) en base p1, avec le paquet de
%  production SEUL, sauvegardees dans Z10_closures_basis_prod.mat pour la garde 1
%  de RUN_Z10_CLOSURES_BASIS. Configuration de RUN_Z1_CAVITY (vrillage harmonique
%  neutralise, X_m0 a 0,2 A, s = 1e-4 et 0,0188), N_h = 8192.
clear; clc; t0=tic;
Nh=8192; s_ch=0.0188; Im_ref=0.2;
here=fileparts(mfilename('fullpath')); if isempty(here), here=pwd; end
if contains(path,fullfile(here,'variant_p1a')), error('Z10_PROD : la variante p1a est sur le path'); end
CAS={'Carter',0,0,'p1',false; 'Phi_O=0',17,4,'p1',false; 'Phi_O=0',33,16,'p1',false; 'Phi_O=0',33,32,'p1',false; ...
     'cavites',17,4,'p1',true; 'cavites',33,16,'p1',true; 'cavites',33,32,'p1',true};
nc=size(CAS,1); V=nan(nc,10); MK=false(nc,1);
for k=1:nc
    M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
    ctx=mec.build_context(M); G=ctx.G; ctx.M.opt.skew_harm=0;
    if strcmp(CAS{k,1},'Carter'), ctx.AG=ctx.AGcarter;
    elseif CAS{k,5}, cav=load(fullfile(here,sprintf('cavity_nO%d.mat',CAS{k,3}))); ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,CAS{k,2},CAS{k,3},Nh,CAS{k,4},cav); MK(k)=isfield(ctx.AG.AF,'variant_p1a');
    else, ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,CAS{k,2},CAS{k,3},Nh,CAS{k,4},[]); MK(k)=isfield(ctx.AG.AF,'variant_p1a'); end
    Rm=mec.magnetizing(ctx,Im_ref); ctx.Xm0=Rm.Xm;
    r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0); rn=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
    V(k,:)=[ctx.Xm0, r0.I1, r0.E1, rn.Xm, rn.I1, rn.Tem, rn.Pfe, rn.Pcu_r, rn.cosphi, rn.eta];
    fprintf('  Z10_PROD %-8s (%2d,%2d) %s : Xm0 %.3f | I0 %.3f | T %.2f  [variante %d]  (%.0f s)\n',CAS{k,1},CAS{k,2},CAS{k,3},CAS{k,4},V(k,1),V(k,2),V(k,6),MK(k),toc(t0));
end
save(fullfile(here,'Z10_closures_basis_prod.mat'),'V','CAS','MK');
