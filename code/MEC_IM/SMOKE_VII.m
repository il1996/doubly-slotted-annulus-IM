%% SMOKE_VII - essai isole du tableau de cloture (section VII de RUN_ARTICLE)
%
%  Verifie que le tableau s'assemble et que les deux fermetures d'entrefer
%  sont bien comparables terme a terme, SANS relancer les 740 s du script
%  complet ni relire les fichiers .tab d'ANSYS : les references EF sont
%  reprises en dur depuis RUN_ARTICLE (elles y sont mesurees).
clear; clc;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G; p=M.p;
Mc=M; Mc.opt.gap_fourier=0;
ctxC=ctx; ctxC.AG=ctx.AGcarter;

n1=1469.772776; s_ch=(M.ns*60-n1)/(M.ns*60);
r  =mec.equivalent_circuit(ctx ,s_ch,ctx.Xm0);
r0 =mec.equivalent_circuit(ctx ,1e-4,ctx.Xm0);
ctxC.Xm0=mec.magnetizing(ctxC,0.2).Xm;
cA0=mec.equivalent_circuit(ctxC,1e-4,ctxC.Xm0);
cAL=mec.equivalent_circuit(ctxC,s_ch,ctxC.Xm0);

[bgA,~  ]=bg1load(ctxC,cAL);
[bgB,BtB]=bg1load(ctx ,r  );

%  References EF (mesurees dans RUN_ARTICLE sections III et IV)
I0f=8.49; E1f=382.1; Tf7=120.2; bg1f=0.920; btf=0.131;

fprintf('\n===== VII. CLOTURE D''ENTREFER : CARTER vs DtN vs EF =====\n');
fprintf('point de charge : s = %.4f  (n = %.1f tr/min)\n\n',s_ch,n1);
fprintf('%-24s %11s %11s %11s   %s\n','grandeur','CARTER','DtN','EF','ecart DtN');
p7=@(lab,a,b,f,u)fprintf('%-24s %11.3f %11.3f %11.3f   %+7.1f %%  %s\n', ...
    lab,a,b,f,100*(b-f)/f,u);
p7('Xm non sature [ohm]',ctxC.Xm0,ctx.Xm0,NaN,'(pas de ref. EF)');
p7('Xm charge [ohm]',cAL.Xm,r.Xm,46,'');
p7('I0 a vide [A]',cA0.I1,r0.I1,I0f,'');
p7('E1 a vide [V]',cA0.E1,r0.E1,E1f,'');
p7('Couple charge [N.m]',cAL.Tem,r.Tem,Tf7,'');
p7('Bg1 charge [T]',bgA,bgB,bg1f,'(flux/pas)');
fprintf('%-24s %11s %11.3f %11.3f   %+7.1f %%  %s\n','Bt rms charge [T]', ...
    'AUCUN',BtB,btf,100*(BtB-btf)/btf,'(chemin dentaire)');
fprintf('%-24s %11.3f %11.3f %11s\n','Carter kC (implicite)', ...
    ctx.AGcarter.kC,kc_implicite(ctx,M,G),'--');
fprintf('\n=== SMOKE_VII termine ===\n');

%% ---- copies locales des fonctions de RUN_ARTICLE ----
function [bg1,btr]=bg1load(cx,rr)
    Mx=cx.M; Gx=cx.G; Wx=cx.W;
    q1=angle(rr.I1c);
    i3=sqrt(2)*abs(rr.I1c)*[cos(q1);cos(q1-2*pi/3);cos(q1+2*pi/3)];
    Fs=Wx.slotMMF(i3);
    q2=angle(rr.I2c); tb=2*pi*(0:Mx.Nr-1)/Mx.Nr;
    Fu=cumsum(cos(Mx.p*tb-q2).'); Fu=Fu-mean(Fu);
    c1=(2/Mx.Nr)*sum(Fu.'.*exp(-1j*Mx.p*tb));
    Fr=-Fu*((3/2)*(4/pi)*(Wx.kw1*Wx.Nph/(2*Mx.p))*sqrt(2)*abs(rr.I2c)/abs(c1));
    S=mec.solve_network(cx.net,Gx,cx.BH,cx.AG,Fs,Fr,Mx.opt);
    th=2*pi*(0:Mx.Ns-1)/Mx.Ns; Bg=S.Bgap_avg_i(:).';
    bg1=abs((2/Mx.Ns)*sum(Bg.*exp(-1j*Mx.p*th)));
    btr=NaN;
    if isfield(S,'Usurf') && isfield(cx.AG,'expand')
        tq=linspace(0,2*pi,2001); tq(end)=[];
        [~,Bt]=cx.AG.field(S.Usurf,0.5*(Gx.Rs+Gx.Rr),tq);
        btr=sqrt(mean(Bt(:).^2));
    end
end
function kc=kc_implicite(cx,Mx,Gx)
    ths=2*pi*(0:Mx.Ns-1)/Mx.Ns; thr=2*pi*(0:Mx.Nr-1)/Mx.Nr;
    AL=mec.airgap_fourier(ths,repmat(2*pi/Mx.Ns,1,Mx.Ns), ...
                          thr,repmat(2*pi/Mx.Nr,1,Mx.Nr), ...
                          Gx.Rs,Gx.Rr,Mx.L,cx.AG.Nh);
    cl=cx; cl.AG=AL;
    kc=mec.magnetizing(cl,0.2).Xm/cx.Xm0;
end
