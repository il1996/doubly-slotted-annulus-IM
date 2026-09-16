%% RUN_Z7_NOLOAD_NET - essai a vide NUMERIQUE, cote reseau : balayage de la
%  tension d'alimentation a vide (s = 1e-4) pour trois fermetures
%  d'entrefer, avec la MEME construction de X_m que celle appliquee aux
%  exports du transitoire EF (voir T2_ansys/noload_postprocess.py) :
%
%      E_2D = U_ph - R_s*I - j*X_ext*I          (phaseurs, fondamental)
%      X_m,2D = |E_2D| / |I|
%
%  ou X_ext est la reactance serie DECLAREE hors du champ 2D de chaque
%  modele : cote EF l'inductance externe de tetes du projet (6,2511 mH =
%  1,9638 ohm) ; cote reseau la somme des termes de tetes du reseau
%  (X_s,ew empirique 0,7603 + residu 3D identifie 0,4260 = 1,1863 ohm).
%  E_2D contient donc, des deux cotes, la chute de fuite 2D (encoche + bec)
%  et le champ des harmoniques d'espace, en plus de la f.e.m. d'entrefer.
%
%  Le reseau fournit aussi E1 (f.e.m. d'entrefer) et X_m = E1/I_m (saturee),
%  ainsi que X_m0 a I_m = 0,2 A (definition de la Table 8).
%  Configuration loyale de Z1 (vrillage harmonique neutralise), N_h = 8192.
%  Sortie : Z7_noload_net_out.txt, Z7_noload_net.mat
clear; clc; t0=tic;
if isfile('Z7_noload_net_out.txt'), delete('Z7_noload_net_out.txt'); end
diary('Z7_noload_net_out.txt'); diary on;
Nh=8192; s0=1e-4; Im_ref=0.2;
VL=[690 550 410 275 140 100 60];
cases={ 'Carter',                 [];
        'operateur (17,4) Phi_O=0', struct('nT',17,'nO',4,'cav',false);
        'operateur+cavite (33,16)', struct('nT',33,'nO',16,'cav',true)};
nc=size(cases,1); nv=numel(VL);
fprintf('=== Z7 : essai a vide numerique, cote reseau (s = %.0e, skew_harm = 0, N_h = %d) ===\n\n',s0,Nh);
OUT=nan(nc,nv,9); XM0=nan(nc,1); XEXT=nan(nc,1);
for k=1:nc
    M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
    ctx=mec.build_context(M); G=ctx.G; ctx.M.opt.skew_harm=0;
    cfg=cases{k,2};
    if isempty(cfg)
        ctx.AG=ctx.AGcarter;
    elseif ~cfg.cav
        ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,cfg.nT,cfg.nO,Nh,'p1',[]);
    else
        cav=load(sprintf('cavity_nO%d.mat',cfg.nO));
        ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,cfg.nT,cfg.nO,Nh,'p1',cav);
    end
    Rm=mec.magnetizing(ctx,Im_ref); ctx.Xm0=Rm.Xm; XM0(k)=Rm.Xm;
    Xext=ctx.Lk.Xs_ew+ctx.Lk.Xs_end3D; XEXT(k)=Xext;
    fprintf('  %-28s X_m0(0.2 A) = %.3f ohm | R_s %.4f | X_ext(reseau) = %.4f ohm (ew %.4f + 3D %.4f) | X_slot+tip = %.4f (%.0f s)\n', ...
        cases{k,1},XM0(k),ctx.Rs,Xext,ctx.Lk.Xs_ew,ctx.Lk.Xs_end3D,ctx.Lk.Xs_slot+ctx.Lk.Xs_tip,toc(t0));
    fprintf('  %8s %9s %9s %9s %9s %9s %9s %10s %10s %9s\n','V_line','U_ph','I0 (A)','phi_I(deg)','E1 (V)','Im (A)','Xm=E1/Im','|E_2D| (V)','X_m,2D','Bg1 (T)');
    for v=1:nv
        ctx.M.Uph=VL(v)/sqrt(3);
        r0=mec.equivalent_circuit(ctx,s0,ctx.Xm0);
        I=r0.I1c; Uph=ctx.M.Uph;
        E2D=Uph - ctx.Rs*I - 1i*Xext*I;
        Xm2D=abs(E2D)/abs(I);
        OUT(k,v,:)=[Uph abs(I) angle(I)*180/pi r0.E1 r0.Im r0.Xm abs(E2D) Xm2D r0.Bg1];
        fprintf('  %8.1f %9.3f %9.4f %9.3f %9.3f %9.4f %9.3f %10.3f %10.3f %9.4f\n',VL(v),OUT(k,v,:));
    end
    fprintf('\n');
end
%  courbe de magnetisation pure du reseau : Xm(Im) = E1/Im, aucun circuit
fprintf('  ---- courbe X_m(I_m) du reseau seul (mec.magnetizing), trois fermetures ----\n');
IMS=[0.2 0.5 1 2 3 4 5 6 7 8 9 10 12];
XMI=nan(nc,numel(IMS));
for k=1:nc
    M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
    ctx=mec.build_context(M); G=ctx.G;
    cfg=cases{k,2};
    if isempty(cfg), ctx.AG=ctx.AGcarter;
    elseif ~cfg.cav, ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,cfg.nT,cfg.nO,Nh,'p1',[]);
    else, cav=load(sprintf('cavity_nO%d.mat',cfg.nO)); ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,cfg.nT,cfg.nO,Nh,'p1',cav); end
    for v=1:numel(IMS), XMI(k,v)=mec.magnetizing(ctx,IMS(v)).Xm; end
    fprintf('  %-28s',cases{k,1}); fprintf(' %8.3f',XMI(k,:)); fprintf('\n');
end
fprintf('  %-28s','I_m (A)'); fprintf(' %8.2f',IMS); fprintf('\n');
save('Z7_noload_net.mat','VL','cases','OUT','XM0','XEXT','IMS','XMI');
fprintf('\n  duree %.0f s\n=== Z7 termine ===\n',toc(t0));
diary off;
