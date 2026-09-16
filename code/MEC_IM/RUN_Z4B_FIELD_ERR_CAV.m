%% RUN_Z4B_FIELD_ERR_CAV - indicateurs d'erreur des formes d'onde de champ
%  avec la fermeture A CAVITES (33,16) : reunit, pour cette fermeture, ce que
%  RUN_M11_FIELD_ERR (indicateurs a phi = 0, energie d'erreur par famille),
%  RUN_M11B_ROTORPOS (balayage de la position rotor sur un pas dentaire,
%  dispersion, position appariee) et RUN_M11C_MATCHED (indicateurs a la
%  position appariee) produisent pour la fermeture Phi_O = 0 (17,4).
%  Meme chaine (copie conforme de RUN_B1_IM_P1 sec. 4 / RUN_Z2), meme
%  recalage sur la phase exacte du fondamental, memes gardes.
%  Sortie : Z4B_field_err_cav_out.txt, Z4B_field_err_cav.mat
clear; clc; t0=tic;
if isfile('Z4B_field_err_cav_out.txt'), delete('Z4B_field_err_cav_out.txt'); end
diary('Z4B_field_err_cav_out.txt'); diary on;
ROOT=fullfile('..','..','reference','ANSYS_18_5kW');
M=mec.machine_18_5kW(); M.opt.verbose=false; ctx0=mec.build_context(M); G=ctx0.G; W=ctx0.W;
p=M.p; Rm=0.5*(G.Rs+G.Rr); Nr=M.Nr; Ns=M.Ns;
nT=33; nO=16; Nh=8192; s_ch=0.0188;
thq=linspace(0,2*pi,2001); thq(end)=[];
cav=load(sprintf('cavity_nO%d.mat',nO));
NPHI=23; phis=linspace(0,2*pi/Nr,NPHI);
fprintf('=== Z4B : indicateurs d''erreur des formes d''onde, fermeture a CAVITES ===\n\n');
fprintf('  operateur : mec.airgap_dtn_tooth_cav(nT=%d, nO=%d, N_h=%d, base p1, cavity_nO%d.mat)\n',nT,nO,Nh,nO);
fprintf('  champ     : AT.field(Usurf, R_moy, theta) ; reference : %s\n',ROOT);
fprintf('  balayage  : %d positions rotor sur un pas dentaire (%.4f deg mec)\n\n',NPHI,360/Nr);

P0=read_profile(fullfile(ROOT,'transitoire','a vide','Calculator Expressions Plot 1.tab'));
P1=read_profile(fullfile(ROOT,'transitoire','en charge','Calculator Expressions Plot 1.tab'));
cf=@(y,th,k)(2/numel(y))*sum(y(:).*exp(-1j*k*th(:)));
NG=numel(P0.th); nmax=floor(NG/2); wP=0.5*ones(nmax,1); if mod(NG,2)==0, wP(end)=0.25; end
nn=(1:nmax).'; fS=false(nmax,1); fR=false(nmax,1);
for k=1:ceil(nmax/Ns), fS=fS|(nn==k*Ns-p)|(nn==k*Ns+p); end
for k=1:ceil(nmax/Nr), fR=fR|(nn==k*Nr-p)|(nn==k*Nr+p); end
ff0=(nn==p);
FAM={'fondamental',ff0; 'denture stator seule',fS&~fR&~ff0; 'denture rotor seule',fR&~fS&~ff0; ...
     'commune aux deux',fS&fR&~ff0; 'reste du spectre',~fS&~fR&~ff0};
NOM={'a vide','en charge'}; REF=[0.942 0.128; 0.920 0.131];

%% ---- balayage de la position rotor (indicateurs complets a chaque phi) ----
R=nan(NPHI,2); RT=nan(NPHI,2); FLD=nan(NPHI,4); DL=nan(NPHI,2); EMAX=nan(NPHI,2);
ENE=nan(NPHI,2,5); DPH=nan(NPHI,2); XM0=nan(NPHI,1); PTS=nan(NPHI,4);
fprintf('  %4s %9s | %10s %10s | %10s %10s | %9s %9s\n','k','phi(deg)','RMSE vide','RMSE chg','Bg1 vide','Bg1 chg','Bt vide','Bt chg');
for q=1:NPHI
    phi=phis(q);
    A=mec.airgap_dtn_tooth_cav(M,G,phi,nT,nO,Nh,'p1',cav);
    ctx=ctx0; ctx.AG=A; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm; XM0(q)=ctx.Xm0;
    r =mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
    r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
    PTS(q,:)=[r0.I1 r.I1 r.Tem r.Pcu_r];
    Rnl=mec.magnetizing(ctx,r0.Im); U0=Rnl.S.Usurf;
    p1=angle(r.I1c);
    i3=sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
    Fs=W.slotMMF(i3);
    p2=angle(r.I2c); tb=2*pi*(0:Nr-1)/Nr + phi;
    Fu=cumsum(cos(M.p*tb-p2).'); Fu=Fu-mean(Fu);
    c1=(2/Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
    Fr=-Fu*((3/2)*(4/pi)*(W.kw1*W.Nph/(2*M.p))*sqrt(2)*abs(r.I2c)/abs(c1));
    SL=mec.solve_network(ctx.net,G,ctx.BH,A,Fs,Fr,M.opt); U1=SL.Usurf;
    [b0,t0f]=A.field(U0,Rm,thq); [b1,t1f]=A.field(U1,Rm,thq);
    ffn=@(y)abs((2/numel(y))*sum(y(:).'.*exp(-1j*p*thq)));
    FLD(q,:)=[ffn(b0) sqrt(mean(t0f(:).^2)) ffn(b1) sqrt(mean(t1f(:).^2))];
    for c=1:2
        if c==1, U=U0; P=P0; else, U=U1; P=P1; end
        ca=cf(P.Br,P.th,p); phA=atan2(-imag(ca),real(ca));
        Uf=A.expand(U);
        [~,~,H]=A.AF.field(Uf(1:A.Msf),Uf(A.Msf+1:end),Rm,0);
        dlt=(atan2(H.Brs(p),H.Brc(p))-phA)/p;
        [Bm,Btm]=A.field(U,Rm,P.th(:)+dlt); Bm=Bm(:); Btm=Btm(:);
        cmv=cf(Bm,P.th,p); DPH(q,c)=abs(mod(atan2(-imag(cmv),real(cmv))-phA+pi,2*pi)-pi);
        e=Bm-P.Br(:); et=Btm-P.Bt(:);
        R(q,c)=sqrt(mean(e.^2)); RT(q,c)=sqrt(mean(et.^2)); EMAX(q,c)=max(abs(e)); DL(q,c)=dlt*180/pi;
        d=arrayfun(@(n)cf(e,P.th,n),(1:nmax).'); Emse=mean(e.^2);
        for z=1:5, m=FAM{z,2}; ENE(q,c,z)=100*sum(wP(m).*abs(d(m)).^2)/Emse; end
        if q==1
            Erec=sum(wP.*abs(d).^2)+mean(e)^2; PARS(c)=abs(Erec-Emse)/Emse; %#ok<SAGROW>
            ERR0{c}=e; A0=A; %#ok<SAGROW>
        end
    end
    fprintf('  %4d %9.4f | %10.5f %10.5f | %10.5f %10.5f | %9.5f %9.5f   (%.0f s)\n', ...
        q,phi*180/pi,R(q,1),R(q,2),FLD(q,1),FLD(q,3),FLD(q,2),FLD(q,4),toc(t0));
end

%% ---- indicateurs a phi = 0 (position du reseau), comme M-11 ----
fprintf('\n  ---- INDICATEURS A phi = 0 (position du reseau) ----\n');
rr0=[sqrt(mean(P0.Br.^2)) sqrt(mean(P1.Br.^2))]; rt0=[sqrt(mean(P0.Bt.^2)) sqrt(mean(P1.Bt.^2))];
fprintf('  %-10s %10s %10s %10s %10s %10s %10s\n','cas','Bg1','Bt rms','RMSE Br','RMSE Br %','RMSE Bt %','err max %');
for c=1:2
    fprintf('  %-10s %10.4f %10.4f %10.4f %9.1f %% %9.1f %% %9.1f %%\n',NOM{c},FLD(1,2*c-1),FLD(1,2*c), ...
        R(1,c),100*R(1,c)/rr0(c),100*RT(1,c)/rt0(c),100*EMAX(1,c)/max(abs(tern(c==1,P0.Br,P1.Br))));
    fprintf('             Bg1 vs EF %.4f : %+.1f %% | Bt rms vs EF %.4f : %+.1f %%\n',REF(c,1),100*(FLD(1,2*c-1)-REF(c,1))/REF(c,1),REF(c,2),100*(FLD(1,2*c)-REF(c,2))/REF(c,2));
    fprintf('             recalage %+8.4f deg mec ; residu des deux estimateurs de phase %.1e rad ; Parseval %.1e\n',DL(1,c),DPH(1,c),PARS(c));
    fprintf('             energie d''erreur : ');
    for z=1:5, fprintf('%s %.1f %% | ',FAM{z,1},ENE(1,c,z)); end
    fprintf('\n');
end
%  localisation ouvertures / faces (masques du pavage) a phi = 0
mS=arc_mask(A0.ths,A0.dths,A0.isFace(1:A0.Msf),P0.th(:));
mR=arc_mask(A0.thr,A0.dthr,A0.isFace(A0.Msf+1:end),P0.th(:));
fprintf('  localisation de l''erreur (phi = 0) :\n  %-11s %-22s %10s %10s %6s\n','cas','zone','rms err','part MSE','pts');
for c=1:2
    e=ERR0{c}; tot=sum(e.^2);
    Z={'ouverture stator',~mS; 'face de dent stator',mS; 'ouverture rotor',~mR; 'face de dent rotor',mR};
    for z=1:size(Z,1), m=Z{z,2}; fprintf('  %-11s %-22s %10.4f %9.1f %% %6d\n',NOM{c},Z{z,1},sqrt(mean(e(m).^2)),100*sum(e(m).^2)/tot,sum(m)); end
end

%% ---- gardes du balayage ----
dper=max(abs(R(end,:)-R(1,:))); dfld=max(abs(FLD(end,:)-FLD(1,:)));
fprintf('\n  GARDE periodicite : RMSE(2pi/Nr)-RMSE(0) = %.2e T (seuil 1e-4) ; champs %.2e T  -> %s\n',dper,dfld,tern(dper<1e-4,'PASSEE','ECHOUEE'));
fprintf('  GARDE recalage    : residu max des deux estimateurs de phase %.1e rad (seuil 1e-3) -> %s\n',max(DPH(:)),tern(max(DPH(:))<1e-3,'PASSEE','ECHOUEE'));

%% ---- position appariee ----
fprintf('\n  ---- POSITION ROTOR APPARIEE (minimum de la RMSE de B_r sur le pas) ----\n');
POS=nan(1,2); RMIN=nan(1,2); SEP=nan(1,2);
for c=1:2
    y=R(1:end-1,c); x=phis(1:end-1).'; [~,im]=min(y);
    ia=mod(im-2,numel(y))+1; ib=mod(im,numel(y))+1; ya=y(ia); yb=y(im); yc=y(ib); h=x(2)-x(1);
    den=(ya-2*yb+yc);
    if abs(den)>eps, dd=0.5*(ya-yc)/den*h; vmin=yb-(ya-yc)^2/(8*den); else, dd=0; vmin=yb; end
    POS(c)=(x(im)+dd)*180/pi; RMIN(c)=vmin;
    loc=find(y<circshift(y,1) & y<circshift(y,-1)); loc=loc(loc~=im);
    if isempty(loc), SEP(c)=Inf; else, SEP(c)=100*(min(y(loc))-yb)/yb; end
    fprintf('    %-10s minimum en phi = %7.4f deg mec ; RMSE %.5f T (phi=0) -> %.5f T (%+.1f %%) ; second minimum local +%.1f %%\n', ...
        NOM{c},POS(c),R(1,c),RMIN(c),100*(RMIN(c)-R(1,c))/R(1,c),SEP(c));
end
fprintf('  dispersion sur le pas dentaire rotor (etendue/moyenne) :\n');
NF={'Bg1 a vide','Bt rms vide','Bg1 charge','Bt rms charge'};
SPR=nan(1,4);
for k=1:4, SPR(k)=100*(max(FLD(:,k))-min(FLD(:,k)))/mean(FLD(:,k)); fprintf('    %-16s %.5f a phi=0 | dispersion %.2f %%\n',NF{k},FLD(1,k),SPR(k)); end
fprintf('    X_m0 : %.3f a phi=0 | dispersion %.2f %% ; I0 %.4f | T %.3f | I1 %.4f a phi=0, dispersions %.2f / %.2f / %.2f %%\n', ...
    XM0(1),100*(max(XM0)-min(XM0))/mean(XM0),PTS(1,1),PTS(1,3),PTS(1,2), ...
    100*(max(PTS(:,1))-min(PTS(:,1)))/mean(PTS(:,1)),100*(max(PTS(:,3))-min(PTS(:,3)))/mean(PTS(:,3)),100*(max(PTS(:,2))-min(PTS(:,2)))/mean(PTS(:,2)));

%% ---- indicateurs a la position appariee (rejoue la chaine a POS) ----
fprintf('\n  ---- INDICATEURS A LA POSITION APPARIEE ----\n');
RESM=nan(2,4); ENEM=nan(2,5); FLDM=nan(2,2);
for c=1:2
    phi=POS(c)*pi/180;
    A=mec.airgap_dtn_tooth_cav(M,G,phi,nT,nO,Nh,'p1',cav);
    ctx=ctx0; ctx.AG=A; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
    r =mec.equivalent_circuit(ctx,s_ch,ctx.Xm0); r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
    if c==1
        Rnl=mec.magnetizing(ctx,r0.Im); U=Rnl.S.Usurf; P=P0;
    else
        p1=angle(r.I1c); i3=sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)]; Fs=W.slotMMF(i3);
        p2=angle(r.I2c); tb=2*pi*(0:Nr-1)/Nr + phi; Fu=cumsum(cos(M.p*tb-p2).'); Fu=Fu-mean(Fu);
        c1=(2/Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
        Fr=-Fu*((3/2)*(4/pi)*(W.kw1*W.Nph/(2*M.p))*sqrt(2)*abs(r.I2c)/abs(c1));
        SL=mec.solve_network(ctx.net,G,ctx.BH,A,Fs,Fr,M.opt); U=SL.Usurf; P=P1;
    end
    [bq,tq]=A.field(U,Rm,thq); FLDM(c,:)=[abs((2/numel(bq))*sum(bq(:).'.*exp(-1j*p*thq))) sqrt(mean(tq(:).^2))];
    ca=cf(P.Br,P.th,p); phA=atan2(-imag(ca),real(ca));
    Uf=A.expand(U); [~,~,H]=A.AF.field(Uf(1:A.Msf),Uf(A.Msf+1:end),Rm,0);
    dlt=(atan2(H.Brs(p),H.Brc(p))-phA)/p;
    [Bm,Btm]=A.field(U,Rm,P.th(:)+dlt); Bm=Bm(:); Btm=Btm(:);
    e=Bm-P.Br(:); et=Btm-P.Bt(:); rmse=sqrt(mean(e.^2)); rref=sqrt(mean(P.Br.^2));
    RESM(c,:)=[rmse 100*rmse/rref sqrt(mean(et.^2)) 100*sqrt(mean(et.^2))/sqrt(mean(P.Bt.^2))];
    d=arrayfun(@(n)cf(e,P.th,n),(1:nmax).'); Emse=mean(e.^2);
    for z=1:5, m=FAM{z,2}; ENEM(c,z)=100*sum(wP(m).*abs(d(m)).^2)/Emse; end
    fprintf('    %-10s phi = %.4f deg : RMSE Br %.5f T (%.1f %%) [estimation parabolique %.5f] | RMSE Bt %.5f T (%.1f %%) | Bg1 %.4f | Bt rms %.4f\n', ...
        NOM{c},POS(c),rmse,RESM(c,2),RMIN(c),RESM(c,3),RESM(c,4),FLDM(c,1),FLDM(c,2));
    fprintf('               energie d''erreur : '); for z=1:5, fprintf('%s %.1f %% | ',FAM{z,1},ENEM(c,z)); end; fprintf('\n');
end

%% ---- synthese pour la colonne "cavites" de la Table 11 ----
fprintf('\n  ---- SYNTHESE (colonne cavites de la Table 11) ----\n');
fprintf('  %-52s %-26s %-26s\n','grandeur','a vide','en charge');
fprintf('  %-52s %-26s %-26s\n','fondamental Bg1 (T), reseau / EF',sprintf('%.4f / %.4f (%+.1f %%)',FLD(1,1),REF(1,1),100*(FLD(1,1)-REF(1,1))/REF(1,1)),sprintf('%.4f / %.4f (%+.1f %%)',FLD(1,3),REF(2,1),100*(FLD(1,3)-REF(2,1))/REF(2,1)));
fprintf('  %-52s %-26s %-26s\n','Bt rms (T), reseau / EF',sprintf('%.4f / %.4f (%+.1f %%)',FLD(1,2),REF(1,2),100*(FLD(1,2)-REF(1,2))/REF(1,2)),sprintf('%.4f / %.4f (%+.1f %%)',FLD(1,4),REF(2,2),100*(FLD(1,4)-REF(2,2))/REF(2,2)));
fprintf('  %-52s %-26s %-26s\n','RMSE Br, position du reseau (T, %% du rms EF)',sprintf('%.3f T, %.1f %%',R(1,1),100*R(1,1)/rr0(1)),sprintf('%.3f T, %.1f %%',R(1,2),100*R(1,2)/rr0(2)));
fprintf('  %-52s %-26s %-26s\n','RMSE Br, position appariee (%% du rms EF)',sprintf('%.1f %% (phi %.3f deg)',RESM(1,2),POS(1)),sprintf('%.1f %% (phi %.3f deg)',RESM(2,2),POS(2)));
fprintf('  %-52s %-26s %-26s\n','dispersion Bg1, Bt rms sur un pas rotor',sprintf('%.2f %%, %.1f %%',SPR(1),SPR(2)),sprintf('%.2f %%, %.1f %%',SPR(3),SPR(4)));
fprintf('  %-52s %-26s %-26s\n','energie d''erreur fond./stator/rotor/commun/reste',sprintf('%.1f/%.1f/%.1f/%.1f/%.1f %%',ENE(1,1,:)),sprintf('%.1f/%.1f/%.1f/%.1f/%.1f %%',ENE(1,2,:)));
save('Z4B_field_err_cav.mat','phis','R','RT','FLD','DL','EMAX','ENE','DPH','XM0','PTS','POS','RMIN','SEP','SPR','RESM','ENEM','FLDM','PARS');
fprintf('\n  duree %.0f s\n=== Z4B termine ===\n',toc(t0));
diary off;

% ======================================================================
function P=read_profile(f)
    fid=fopen(f); hdr=fgetl(fid); fclose(fid);
    tk=regexp(hdr,'"([^"]+)"','tokens'); names=cellfun(@(c)c{1},tk,'UniformOutput',false);
    D=readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
    gc=@(pat)D(:,find(contains(names,pat),1));
    dist=gc('Distance'); A=gc('Flux_Lines'); Br=gc('Br'); Bt=gc('Bt');
    C=dist(end)*1e-3; n=numel(dist)-1;
    P.th=2*pi*dist(1:n)*1e-3/C; P.A=A(1:n); P.Br=Br(1:n); P.Bt=Bt(1:n);
end
function m=arc_mask(th,dth,isFace,thq)
    th=th(:); dth=dth(:); isFace=logical(isFace(:));
    m=false(numel(thq),1); tq=mod(thq(:),2*pi);
    for k=1:numel(th)
        if ~isFace(k), continue; end
        a=mod(th(k)-dth(k)/2,2*pi); b=mod(th(k)+dth(k)/2,2*pi);
        if a<b, m=m|(tq>=a & tq<b); else, m=m|(tq>=a | tq<b); end
    end
end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
