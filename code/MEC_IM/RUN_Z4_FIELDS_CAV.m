%% RUN_Z4_FIELDS_CAV - formes d'onde de champ au mi-entrefer avec la
%  fermeture A CAVITES (33,16), a vide et en charge.  Copie de RUN_Z2_FIELDS
%  ou seule la construction de l'operateur change :
%      A1 = mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1')
%  devient
%      cav = load('cavity_nO16.mat');
%      A1  = mec.airgap_dtn_tooth_cav(M,G,0,33,16,Nh,'p1',cav);
%  Le reste (recalage sur la seule phase du fondamental, aucun autre
%  decalage, memes grilles, meme reference EF) est inchange.
%  Sorties : Z4_fields_cav_avide.txt, Z4_fields_cav_charge.txt (meme format
%  que Z2 : theta, Br_EF, Bt_EF, Br_MEC, Bt_MEC), transcript Z4_fields_cav_out.txt.
clear; clc; t0=tic;
if isfile('Z4_fields_cav_out.txt'), delete('Z4_fields_cav_out.txt'); end
diary('Z4_fields_cav_out.txt'); diary on;
ROOT=fullfile('..','..','reference','ANSYS_18_5kW');
M=mec.machine_18_5kW(); M.opt.verbose=false; ctx=mec.build_context(M); G=ctx.G; W=ctx.W;
p=M.p; Rm=0.5*(G.Rs+G.Rr);
nT=33; nO=16; Nh=8192; s_ch=0.0188;
fprintf('=== Z4 : champs au mi-entrefer, fermeture a cavites (%d,%d), N_h=%d, base p1 ===\n',nT,nO,Nh);
cav=load(sprintf('cavity_nO%d.mat',nO));
A1=mec.airgap_dtn_tooth_cav(M,G,0,nT,nO,Nh,'p1',cav);
fprintf('operateur %dx%d condense depuis %d colonnes (%.0f s)\n',size(A1.Y,1),size(A1.Y,2),A1.Msf+A1.Mrf,toc(t0));
ctx.AG=A1; RmS=mec.magnetizing(ctx,0.2); ctx.Xm0=RmS.Xm;
fprintf('X_m0 = %.3f ohm\n',ctx.Xm0);
r =mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
fprintf('a vide : I0 %.4f A, E1 %.2f V | charge : I1 %.4f A, T %.3f N.m\n',r0.I1,r0.E1,r.I1,r.Tem);
Rnl=mec.magnetizing(ctx,r0.Im); U0=Rnl.S.Usurf;
p1=angle(r.I1c);
i3=sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
Fs=W.slotMMF(i3);
p2=angle(r.I2c); tb=2*pi*(0:M.Nr-1)/M.Nr;
Fu=cumsum(cos(M.p*tb-p2).'); Fu=Fu-mean(Fu);
c1=(2/M.Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
Fr=-Fu*((3/2)*(4/pi)*(W.kw1*W.Nph/(2*M.p))*sqrt(2)*abs(r.I2c)/abs(c1));
SL=mec.solve_network(ctx.net,G,ctx.BH,A1,Fs,Fr,M.opt); U1=SL.Usurf;
P0=read_profile(fullfile(ROOT,'transitoire','a vide','Calculator Expressions Plot 1.tab'));
P1=read_profile(fullfile(ROOT,'transitoire','en charge','Calculator Expressions Plot 1.tab'));
cf=@(y,th,k)(2/numel(y))*sum(y(:).*exp(-1j*k*th(:)));
CAS={'avide',U0,P0; 'charge',U1,P1};
for c=1:2
    P=CAS{c,3}; U=CAS{c,2};
    ca=cf(P.Br,P.th,p); phA=atan2(-imag(ca),real(ca));
    thq=P.th(:).';
    [Br,Bt]=A1.field(U,Rm,thq); Br=Br(:).'; Bt=Bt(:).';
    cm=cf(Br,thq,p); phM=atan2(-imag(cm),real(cm));
    dphi=(phM-phA)/p;                     % decalage mecanique a appliquer
    [Br2,Bt2]=A1.field(U,Rm,thq+dphi); Br2=Br2(:).'; Bt2=Bt2(:).';
    cm2=cf(Br2,thq,p); phM2=atan2(-imag(cm2),real(cm2));
    fprintf('%s : phase EF %.4f rad, MEC %.4f -> %.4f apres recalage (dphi=%.4f deg mec)\n',CAS{c,1},phA,phM,phM2,dphi*180/pi);
    rmse=sqrt(mean((Br2-P.Br(:).').^2)); fprintf('   RMSE Br = %.4f T (%.1f %% du rms EF %.4f)\n',rmse,100*rmse/sqrt(mean(P.Br.^2)),sqrt(mean(P.Br.^2)));
    rmst=sqrt(mean((Bt2-P.Bt(:).').^2)); fprintf('   RMSE Bt = %.4f T (%.1f %% du rms EF %.4f)\n',rmst,100*rmst/sqrt(mean(P.Bt.^2)),sqrt(mean(P.Bt.^2)));
    fprintf('   Bg1 MEC %.4f T (EF %.4f) | Bt rms MEC %.4f T (EF %.4f)\n',abs(cm2),abs(ca),sqrt(mean(Bt2.^2)),sqrt(mean(P.Bt.^2)));
    out=[thq(:), P.Br(:), P.Bt(:), Br2(:), Bt2(:)];
    dlmwrite(sprintf('Z4_fields_cav_%s.txt',CAS{c,1}),out,'delimiter','\t','precision','%.6e');
end
fprintf('duree %.0f s\n=== Z4 termine ===\n',toc(t0));
diary off;

function P=read_profile(f)
    fid=fopen(f); hdr=fgetl(fid); fclose(fid);
    tk=regexp(hdr,'"([^"]+)"','tokens'); names=cellfun(@(c)c{1},tk,'UniformOutput',false);
    D=readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    gc=@(pat)D(:,find(contains(names,pat),1));
    dist=gc('Distance'); A=gc('Flux_Lines'); Br=gc('Br'); Bt=gc('Bt');
    C=dist(end)*1e-3; n=numel(dist)-1;
    P.th=2*pi*dist(1:n)*1e-3/C; P.A=A(1:n); P.Br=Br(1:n); P.Bt=Bt(1:n);
end
