%% RUN_Z2_FIELDS - export des formes d'onde de champ au mi-entrefer, reseau
%  (operateur condense, chaine de RUN_B1_IM_P1 / RUN_M11) et reference EF,
%  a vide et en charge, pour la figure de l'article.  Recalage sur la phase
%  du fondamental (comme M11), aucun autre decalage.
clear; clc; t0=tic;
ROOT=fullfile('..','..','reference','ANSYS_18_5kW');
M=mec.machine_18_5kW(); M.opt.verbose=false; ctx=mec.build_context(M); G=ctx.G; W=ctx.W;
p=M.p; Rm=0.5*(G.Rs+G.Rr);
nT=17; nO=4; Nh=8192; s_ch=0.0188;
A1=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
ctx.AG=A1; RmS=mec.magnetizing(ctx,0.2); ctx.Xm0=RmS.Xm;
r =mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
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
    Uf=A1.expand(U);
    A=A1.AF.amp(Uf(1:A1.Msf),Uf(A1.Msf+1:end));
    %  phase exacte du fondamental de B_r du reseau (rang p) : B_r,n ~ -(a_n cos + b_n sin) ...
    thq=P.th(:).';
    [Br,Bt]=A1.field(U,Rm,thq); Br=Br(:).'; Bt=Bt(:).';
    cm=cf(Br,thq,p); phM=atan2(-imag(cm),real(cm));
    dphi=(phM-phA)/p;                     % decalage mecanique a appliquer
    [Br2,Bt2]=A1.field(U,Rm,thq+dphi); Br2=Br2(:).'; Bt2=Bt2(:).';
    cm2=cf(Br2,thq,p); phM2=atan2(-imag(cm2),real(cm2));
    fprintf('%s : phase EF %.4f rad, MEC %.4f -> %.4f apres recalage (dphi=%.4f deg mec)\n',CAS{c,1},phA,phM,phM2,dphi*180/pi);
    rmse=sqrt(mean((Br2-P.Br(:).').^2)); fprintf('   RMSE Br = %.4f T (%.1f %% du rms EF %.4f)\n',rmse,100*rmse/sqrt(mean(P.Br.^2)),sqrt(mean(P.Br.^2)));
    out=[thq(:), P.Br(:), P.Bt(:), Br2(:), Bt2(:)];
    dlmwrite(sprintf('Z2_fields_%s.txt',CAS{c,1}),out,'delimiter','\t','precision','%.6e');
end
fprintf('duree %.0f s\n',toc(t0));

function P=read_profile(f)
    fid=fopen(f); hdr=fgetl(fid); fclose(fid);
    tk=regexp(hdr,'"([^"]+)"','tokens'); names=cellfun(@(c)c{1},tk,'UniformOutput',false);
    D=readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    gc=@(pat)D(:,find(contains(names,pat),1));
    dist=gc('Distance'); A=gc('Flux_Lines'); Br=gc('Br'); Bt=gc('Bt');
    C=dist(end)*1e-3; n=numel(dist)-1;
    P.th=2*pi*dist(1:n)*1e-3/C; P.A=A(1:n); P.Br=Br(1:n); P.Bt=Bt(1:n);
end
