%% RUN_Z1_CAVITY - fermetures d'entrefer : Carter, operateur (Phi_O = 0),
%  operateur + cavites d'encoche.  Configuration loyale (vrillage harmonique
%  neutralise), pavage et troncature declares, references EF relues des .tab.
%
%  Les matrices de cavite cav.Qs / cav.Qr sont calculees par elements finis
%  (Python, cavity.py) dans le repere local d'une encoche et chargees depuis
%  cavity_nO<n>.mat.  Voir mec.airgap_dtn_tooth_cav.
clear; clc; t0=tic;
diary('Z1_cavity_out.txt'); diary on;
Nh=8192; s_ch=0.0188; Im_ref=0.2;
ROOT=fullfile('..','..','reference','ANSYS_18_5kW');
rd=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
rmsw=@(A,c,tw)sqrt(mean(A(A(:,1)>=tw,c).^2));
avgw=@(A,c,tw)mean(A(A(:,1)>=tw,c));

fprintf('=== Z1 : fermetures d''entrefer, avec et sans cavites d''encoche ===\n\n');
%% ---- references EF, relues (fenetre t >= 1,0 s) -----------------------
D0=fullfile(ROOT,'transitoire','a vide'); DC=fullfile(ROOT,'transitoire','en charge');
tw=1.0;
A=rd(fullfile(D0,'Winding Plot 4.tab')); I0F=mean([rmsw(A,2,tw) rmsw(A,3,tw) rmsw(A,4,tw)]);
A=rd(fullfile(D0,'Winding Plot 2.tab')); E0F=mean([rmsw(A,2,tw) rmsw(A,3,tw) rmsw(A,4,tw)]);
A=rd(fullfile(D0,'Loss Plot 1.tab'));    Pfe0F=avgw(A,2,tw);
A=rd(fullfile(DC,'Plot 1.tab'));         TF=avgw(A,3,tw);
A=rd(fullfile(DC,'Winding Plot 4.tab')); I1F=mean([rmsw(A,2,tw) rmsw(A,3,tw) rmsw(A,4,tw)]);
A=rd(fullfile(DC,'Winding Plot 2.tab')); E1F=mean([rmsw(A,2,tw) rmsw(A,3,tw) rmsw(A,4,tw)]);
A=rd(fullfile(DC,'Loss Plot 1.tab'));    PfeF=avgw(A,2,tw);
A=rd(fullfile(DC,'Loss Plot 2.tab'));    PbarF=1e3*avgw(A,2,tw); PcusF=1e3*avgw(A,3,tw);
A=rd(fullfile(DC,'End Connection Plot 3.tab')); PringF=1e3*avgw(A,2,tw);
A=rd(fullfile(DC,'Speed Plot 1.tab'));   nF=avgw(A,2,tw); sF=1-nF/1500;
fprintf('  references EF (t >= %.1f s) : I0 %.4f A | E(2D) a vide %.2f V | Pfe0 %.1f W\n',tw,I0F,E0F,Pfe0F);
fprintf('     charge : s %.5f | T %.3f N.m | I1 %.4f A | E(2D) %.2f V | Pfe %.1f W | Pbar %.1f W | Pring %.1f W | Pcu_s %.1f W\n',sF,TF,I1F,E1F,PfeF,PbarF,PringF,PcusF);

%% ---- fermetures --------------------------------------------------------
cases={ 'Carter',                 [];
        'operateur (17,4) Phi_O=0', struct('nT',17,'nO',4,'cav',false);
        'operateur+cavite (17,4)',  struct('nT',17,'nO',4,'cav',true);
        'operateur+cavite (33,8)',  struct('nT',33,'nO',8,'cav',true);
        'operateur+cavite (33,16)', struct('nT',33,'nO',16,'cav',true)};
nc=size(cases,1);
V=nan(nc,14);
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
    Rm=mec.magnetizing(ctx,Im_ref); ctx.Xm0=Rm.Xm;
    r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
    rn=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
    V(k,:)=[ctx.Xm0, Rm.Bg1_field, r0.I1, r0.E1, r0.Bg1, rn.Xm, rn.I1, rn.E1, rn.Tem, rn.cosphi, rn.eta, rn.Pfe, rn.Pcu_r, rn.Pcu_s];
    fprintf('  %-28s fait (%.0f s) : Xm0 %.3f | I0 %.3f | T %.2f\n',cases{k,1},toc(t0),ctx.Xm0,r0.I1,rn.Tem);
end

fprintf('\n  %-28s %9s %8s %8s %8s %9s %8s %8s %8s %7s %7s %7s %7s %7s\n','fermeture','Xm0','I0','E1(0)','Bg1(0)','Xm(s)','I1','E1(s)','T','cosphi','eta','Pfe','Pcu_r','Pcu_s');
for k=1:nc
    fprintf('  %-28s %9.3f %8.3f %8.2f %8.4f %9.3f %8.3f %8.2f %8.2f %7.4f %7.4f %7.1f %7.1f %7.1f\n',cases{k,1},V(k,[1 3 4 5 6 7 8 9 10 11 12 13 14]));
end
fprintf('  %-28s %9s %8.3f %8.2f %8s %9s %8.3f %8.2f %8.2f %7s %7s %7.1f %7.1f %7.1f\n','reference EF','---',I0F,E0F,'---','---',I1F,E1F,TF,'---','---',PfeF,PbarF+PringF,PcusF);
fprintf('\n  ecarts a la reference (%%) :\n');
fprintf('  %-28s %8s %8s %8s %8s %8s %8s\n','fermeture','I0','E(0)','I1','T','Pfe','Pcu_r');
for k=1:nc
    fprintf('  %-28s %+8.2f %+8.2f %+8.2f %+8.2f %+8.2f %+8.2f\n',cases{k,1}, ...
        100*(V(k,3)-I0F)/I0F, 100*(V(k,4)-E0F)/E0F, 100*(V(k,7)-I1F)/I1F, 100*(V(k,9)-TF)/TF, ...
        100*(V(k,12)-PfeF)/PfeF, 100*(V(k,13)-(PbarF+PringF))/(PbarF+PringF));
end
save('Z1_cavity.mat','V','cases','I0F','E0F','I1F','E1F','TF','PfeF','PbarF','PringF','PcusF','sF');
fprintf('\n  duree %.0f s\n=== Z1 termine ===\n',toc(t0));
diary off;
