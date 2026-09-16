%% RUN_Z3_LEAKAGE - sensibilite a la reactance de fuite d'extremite : la valeur
%  identifiee (0.426 ohm) contre la valeur declaree par le projet EF (6.251 mH
%  = 1.964 ohm en serie, dont on retranche le terme empirique de tetes deja
%  present dans le reseau pour que le total d'extremite soit celui du projet).
clear; clc; t0=tic;
diary('Z3_leakage_out.txt'); diary on;
Nh=8192; s_ch=0.0188; Im_ref=0.2;
M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
ctx=mec.build_context(M); G=ctx.G;
Lk=mec.leakage(M,G,ctx.W);
fprintf('reseau : Xs_slot %.4f  Xs_tip %.4f  Xs_ew %.4f  Xs_end3D %.4f  total %.4f ohm\n',Lk.Xs_slot,Lk.Xs_tip,Lk.Xs_ew,Lk.Xs_end3D,Lk.Xs_leak);
Xfe=2*pi*50*6.251120231328e-3;
fprintf('projet EF : inductance externe 6.2511 mH = %.4f ohm\n',Xfe);
cav=load('cavity_nO16.mat');
for variant=1:2
    for k=1:2
        M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
        if variant==2, M.opt.Xs_end3D=Xfe-Lk.Xs_ew; end
        ctx=mec.build_context(M); G=ctx.G; ctx.M.opt.skew_harm=0;
        if k==1, ctx.AG=ctx.AGcarter; lab='Carter'; else, ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,33,16,Nh,'p1',cav); lab='cavites (33,16)'; end
        Rm=mec.magnetizing(ctx,Im_ref); ctx.Xm0=Rm.Xm;
        r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0); rn=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
        fprintf('variant %d (%s) %-16s : Xs_leak %.4f | I0 %.3f  E1(0) %.2f | I1 %.3f  T %.2f  cosphi %.4f  Pcu_r %.1f  Pfe %.1f  eta %.4f  (%.0f s)\n', ...
            variant, tern(variant==1,'identifie 0.426','EF 1.964 ohm'), lab, ctx.Lk.Xs_leak, r0.I1, r0.E1, rn.I1, rn.Tem, rn.cosphi, rn.Pcu_r, rn.Pfe, rn.eta, toc(t0));
    end
end
diary off;
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
