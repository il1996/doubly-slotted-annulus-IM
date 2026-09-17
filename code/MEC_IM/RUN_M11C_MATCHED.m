%% RUN_M11C_MATCHED - M-11 : les indicateurs A LA POSITION ROTOR APPARIEE
%
%  M-11b a identifie la position mecanique du rotor de chaque cliche EF, en
%  balayant phi sur un pas dentaire rotor et en minimisant la RMSE. Ce
%  script rejoue la chaine A CES DEUX POSITIONS et emet tout ce que M-11
%  publie : RMSE, erreur ponctuelle maximale, les quatre valeurs de champ du
%  sec. 5.8, l'attribution de l'energie d'erreur par famille harmonique, et
%  la repartition ouvertures / faces de dent.
%
%  CE QUE CELA CHANGE, ET IL FAUT LE DIRE D'ABORD. Le balayage de M-11b ne
%  s'est pas contente de deplacer la RMSE : il a montre que les QUATRE
%  VALEURS DE CHAMP publiees au sec. 5.8 dependent elles aussi de phi, et
%  pour trois d'entre elles bien plus que l'accord annonce. Sur un pas
%  dentaire rotor, le fondamental en charge disperse de 5,35 % pour un
%  accord publie de -2,4 %, et la tangentielle efficace de 35,8 % a vide et
%  44,0 % en charge pour des accords publies de -6,1 % et -9,3 %. Autrement
%  dit, trois des quatre accords du sec. 5.8 sont plus petits que l'effet
%  d'un parametre que la comparaison n'appariait pas. Les valeurs emises
%  ici sont celles qu'il faut publier -- et l'ecart entre les deux jeux est
%  la mesure de ce que l'appariement corrige.
%
%  GARDE 1. A la position appariee, la RMSE doit egaler l'estimation
%           parabolique de M-11b a mieux que la resolution du balayage.
%  GARDE 2. Le recalage sur la phase du fondamental doit rester sans
%           ambiguite : les deux estimateurs de phase concordent a 1e-3 rad.
%  GARDE 3. Parseval : l'energie d'erreur reconstruite par ordre doit
%           reproduire la MSE ponctuelle a 1e-6 relatif.

clear; clc; t0=tic;
if isfile('M11C_matched_out.txt'), delete('M11C_matched_out.txt'); end
diary('M11C_matched_out.txt'); diary on;
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
S=load('M11B_rotorpos.mat');            % POS et RMIN viennent de M-11b

M=mec.machine_18_5kW(); ctx0=mec.build_context(M); G=ctx0.G; W=ctx0.W;
p=M.p; Rm=0.5*(G.Rs+G.Rr); Nr=M.Nr; Ns=M.Ns;
nT=17; nO=4; Nh=8192; s_ch=0.0188;
thq=linspace(0,2*pi,2001); thq(end)=[];

fprintf('=== M-11c : indicateurs a la position rotor appariee ===\n\n');
fprintf('  positions identifiees par M-11b (pas dentaire rotor %.4f deg) :\n',360/Nr);
fprintf('    a vide    phi = %.4f deg mec\n',S.POS(1));
fprintf('    en charge phi = %.4f deg mec\n\n',S.POS(2));

P0=read_profile(fullfile(ROOT,'transitoire','a vide','Calculator Expressions Plot 1.tab'));
P1=read_profile(fullfile(ROOT,'transitoire','en charge','Calculator Expressions Plot 1.tab'));
cf=@(y,th,k)(2/numel(y))*sum(y(:).*exp(-1j*k*th(:)));
NOM={'a vide','en charge'}; REF=[0.942 0.128; 0.920 0.131];
PUB0=[0.9452 0.1201; 0.8983 0.1188];    % les quatre valeurs a phi = 0

RES=nan(2,4); FLD=nan(2,2); ok1=true; ok2=true; ok3=true;
NG=numel(P0.th); nmax=floor(NG/2); wP=0.5*ones(nmax,1);
if mod(NG,2)==0, wP(end)=0.25; end
nn=(1:nmax).'; fS=false(nmax,1); fR=false(nmax,1);
for k=1:ceil(nmax/Ns), fS=fS|(nn==k*Ns-p)|(nn==k*Ns+p); end
for k=1:ceil(nmax/Nr), fR=fR|(nn==k*Nr-p)|(nn==k*Nr+p); end
ff0=(nn==p);

for c=1:2
    phi=S.POS(c)*pi/180;
    A=mec.airgap_dtn_tooth(M,G,phi,nT,nO,Nh,'p1');
    ctx=ctx0; ctx.AG=A; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
    r =mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
    r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
    if c==1
        Rnl=mec.magnetizing(ctx,r0.Im); U=Rnl.S.Usurf; P=P0;
    else
        p1=angle(r.I1c);
        i3=sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
        Fs=W.slotMMF(i3);
        p2=angle(r.I2c); tb=2*pi*(0:Nr-1)/Nr + phi;
        Fu=cumsum(cos(M.p*tb-p2).'); Fu=Fu-mean(Fu);
        c1=(2/Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
        Fr=-Fu*((3/2)*(4/pi)*(W.kw1*W.Nph/(2*M.p))*sqrt(2)*abs(r.I2c)/abs(c1));
        SL=mec.solve_network(ctx.net,G,ctx.BH,A,Fs,Fr,M.opt); U=SL.Usurf; P=P1;
    end
    %  --- les deux valeurs de champ du sec. 5.8 ---
    [bq,tq]=A.field(U,Rm,thq);
    FLD(c,:)=[abs((2/numel(bq))*sum(bq(:).'.*exp(-1j*p*thq))) sqrt(mean(tq(:).^2))];
    %  --- recalage exact et indicateurs ponctuels ---
    ca=cf(P.Br,P.th,p); phA=atan2(-imag(ca),real(ca));
    Uf=A.expand(U);
    [~,~,H]=A.AF.field(Uf(1:A.Msf),Uf(A.Msf+1:end),Rm,0);
    dlt=(atan2(H.Brs(p),H.Brc(p))-phA)/p;
    [Bm,Btm]=A.field(U,Rm,P.th(:)+dlt); Bm=Bm(:); Btm=Btm(:);
    cmv=cf(Bm,P.th,p); dph=mod(atan2(-imag(cmv),real(cmv))-phA+pi,2*pi)-pi;
    ok2=ok2 && abs(dph)<1e-3;
    e=Bm-P.Br(:); rmse=sqrt(mean(e.^2)); rref=sqrt(mean(P.Br.^2));
    et=Btm-P.Bt(:);
    RES(c,:)=[rmse 100*rmse/rref max(abs(e)) 100*max(abs(e))/max(abs(P.Br))];
    ok1=ok1 && abs(rmse-S.RMIN(c))<0.01*S.RMIN(c);

    fprintf('  ================ %s ================\n',upper(NOM{c}));
    fprintf('  recalage %+8.4f deg mec ; residu des deux estimateurs %.1e rad\n',dlt*180/pi,abs(dph));
    fprintf('  RMSE B_r %8.5f T (%.1f %% de la valeur efficace EF %.4f)\n', ...
        rmse,RES(c,2),rref);
    fprintf('  RMSE B_t %8.5f T (%.1f %% de %.4f)\n', ...
        sqrt(mean(et.^2)),100*sqrt(mean(et.^2))/sqrt(mean(P.Bt.^2)),sqrt(mean(P.Bt.^2)));
    fprintf('  erreur ponctuelle max %8.5f T (%.1f %% de la crete)\n',RES(c,3),RES(c,4));
    fprintf('  M-11b annoncait %.5f T par interpolation parabolique -> ecart %.1e\n', ...
        S.RMIN(c),abs(rmse-S.RMIN(c)));
    fprintf('  les deux valeurs de champ du sec. 5.8, ICI et a phi = 0 :\n');
    NF={'Bg1','Bt rms'};
    for q=1:2
        fprintf('    %-8s %8.4f (apparie, %+5.1f %% / EF %.3f) contre %8.4f a phi=0 (%+5.1f %%)\n', ...
            NF{q},FLD(c,q),100*(FLD(c,q)-REF(c,q))/REF(c,q),REF(c,q), ...
            PUB0(c,q),100*(PUB0(c,q)-REF(c,q))/REF(c,q));
    end
    %  --- energie d'erreur par ordre ---
    d=arrayfun(@(n)cf(e,P.th,n),(1:nmax).');
    Emse=mean(e.^2); Erec=sum(wP.*abs(d).^2)+mean(e)^2;
    rel=abs(Erec-Emse)/Emse; ok3=ok3 && rel<1e-6;
    fprintf('  energie d''erreur : MSE %.5e, reconstruite a %.1e relatif\n',Emse,rel);
    FAM={'fondamental',ff0; 'denture stator seule',fS&~fR&~ff0; ...
         'denture rotor seule',fR&~fS&~ff0; 'commune aux deux',fS&fR&~ff0; ...
         'reste du spectre',~fS&~fR&~ff0};
    for z=1:size(FAM,1)
        m=FAM{z,2};
        fprintf('    %-22s %5.1f %% de la MSE (%d ordres)\n', ...
            FAM{z,1},100*sum(wP(m).*abs(d(m)).^2)/Emse,sum(m));
    end
    %  --- ouvertures contre faces de dent ---
    mS=arc_mask(A.ths,A.dths,A.isFace(1:A.Msf),P.th(:));
    mR=arc_mask(A.thr,A.dthr,A.isFace(A.Msf+1:end),P.th(:));
    tot=sum(e.^2);
    Z={'ouverture stator',~mS;'face de dent stator',mS; ...
       'ouverture rotor',~mR;'face de dent rotor',mR};
    for z=1:size(Z,1)
        m=Z{z,2};
        fprintf('    %-22s rms %.4f | %5.1f %% de la MSE sur %d points\n', ...
            Z{z,1},sqrt(mean(e(m).^2)),100*sum(e(m).^2)/tot,sum(m));
    end
    fprintf('\n');
end

fprintf('  ---- GARDES ----\n');
fprintf('    1 (accord avec l''estimation de M-11b)  %s\n',tern(ok1,'PASSEE','ECHOUEE'));
fprintf('    2 (recalage sans ambiguite)             %s\n',tern(ok2,'PASSEE','ECHOUEE'));
fprintf('    3 (Parseval)                            %s\n',tern(ok3,'PASSEE','ECHOUEE'));
save('M11C_matched.mat','RES','FLD','ok1','ok2','ok3');
fprintf('  duree %.0f s\n=== M-11c termine ===\n',toc(t0));
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
