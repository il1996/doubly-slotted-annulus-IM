%% RUN_Z6_COUPLED_CONV_CAV - convergence en pavage de la SOLUTION COUPLEE
%  avec la fermeture A CAVITES (extension de RUN_M5_COUPLED_CONV, qui ne
%  traitait que la condition Phi_O = 0 aux pavages (9,4), (17,4), (33,4)).
%
%  Pour chaque pavage : X_m0 (I_m = 0,2 A), point nominal s = 0,0188 (couple,
%  I1, pertes fer, rendement) et point a vide (I0), avec (a) Phi_O = 0 et
%  (b) cavites.  Base chapeau, N_h = 8192 fixe, vrillage harmonique
%  neutralise, references EF relues des .tab (fenetre t >= 1 s).
%  Sortie : Z6_coupled_conv_cav_out.txt, Z6_coupled_conv_cav.mat
clear; clc; t0=tic;
if isfile('Z6_coupled_conv_cav_out.txt'), delete('Z6_coupled_conv_cav_out.txt'); end
diary('Z6_coupled_conv_cav_out.txt'); diary on;
Nh=8192; s_ch=0.0188;
TIL=[9 4; 17 4; 33 4; 17 8; 33 8; 17 16; 33 16; 65 16];
ROOT=fullfile('..','..','reference','ANSYS_18_5kW');
rd=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
rmsw=@(A,c,tw)sqrt(mean(A(A(:,1)>=tw,c).^2,'omitnan'));
avgw=@(A,c,tw)mean(A(A(:,1)>=tw,c),'omitnan');

fprintf('=== Z6 : convergence de la solution couplee en pavage, Phi_O = 0 et cavites ===\n\n');
fprintf('  base chapeau (P1) | N_h = %d fixe | vrillage harmonique neutralise | s = %.4f\n',Nh,s_ch);
DC=fullfile(ROOT,'transitoire','en charge'); D0=fullfile(ROOT,'transitoire','a vide'); tw=1.0;
TF =avgw(rd(fullfile(DC,'Plot 1.tab')),3,tw);
A=rd(fullfile(DC,'Winding Plot 4.tab')); I1F=mean([rmsw(A,2,tw) rmsw(A,3,tw) rmsw(A,4,tw)]);
PfF=avgw(rd(fullfile(DC,'Loss Plot 1.tab')),2,tw);
A=rd(fullfile(D0,'Winding Plot 4.tab')); I0F=mean([rmsw(A,2,tw) rmsw(A,3,tw) rmsw(A,4,tw)]);
fprintf('  references EF relues : couple %.4f N.m | I1 %.4f A | pertes fer %.2f W | I0 %.4f A\n\n',TF,I1F,PfF,I0F);

n=size(TIL,1); V=nan(n,6,2); NC=nan(n,1);
for k=1:n
    nT=TIL(k,1); nO=TIL(k,2);
    for c=1:2
        M=mec.machine_18_5kW(); M.opt.skew_harm=0; M.opt.verbose=false;
        ctx=mec.build_context(M); G=ctx.G; ctx.M.opt.skew_harm=0;
        if c==1
            ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,nT,nO,Nh,'p1',[]);
        else
            cav=load(sprintf('cavity_nO%d.mat',nO));
            ctx.AG=mec.airgap_dtn_tooth_cav(M,G,0,nT,nO,Nh,'p1',cav);
        end
        NC(k)=ctx.AG.Msf+ctx.AG.Mrf;
        ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
        r=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
        r0=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
        V(k,:,c)=[r.Tem, r.I1, r.Pfe, r.eta, ctx.Xm0, r0.I1];
        fprintf('  pavage (%2d,%2d) %-8s : T %.4f | I1 %.4f | Pfe %.2f | eta %.4f | Xm0 %.3f | I0 %.4f  (%.0f s)\n', ...
            nT,nO,tern(c==1,'Phi_O=0','cavites'),V(k,:,c),toc(t0));
    end
end

LAB={'Phi_O = 0','cavites'};
for c=1:2
    fprintf('\n  ---- SOLUTION COUPLEE, fermeture %s ----\n',LAB{c});
    fprintf('  %8s %8s %12s %12s %12s %12s %12s %10s\n','pavage','colonnes','couple N.m','I1 (A)','P_fe (W)','rendement','X_m0 (ohm)','I0 (A)');
    for k=1:n
        fprintf('  (%2d,%2d) %8d %12.4f %12.4f %12.2f %12.4f %12.3f %10.4f\n',TIL(k,1),TIL(k,2),NC(k),V(k,:,c));
    end
    fprintf('  ecarts a la reference (%%) : couple | I1 | P_fe | I0\n');
    for k=1:n
        fprintf('  (%2d,%2d) %+11.2f %% %+11.2f %% %+11.2f %% %+11.2f %%\n',TIL(k,1),TIL(k,2), ...
            100*(V(k,1,c)-TF)/TF,100*(V(k,2,c)-I1F)/I1F,100*(V(k,3,c)-PfF)/PfF,100*(V(k,6,c)-I0F)/I0F);
    end
    lab={'couple','courant I1','pertes fer','rendement','X_m0','I0'};
    fprintf('  dispersion (etendue/moyenne) sur les %d pavages :',n);
    for j=1:6
        fprintf(' %s %.3f %% |',lab{j},100*(max(V(:,j,c))-min(V(:,j,c)))/mean(V(:,j,c)));
    end
    fprintf('\n');
    i3=find(TIL(:,2)==4);
    fprintf('  dispersion sur les trois pavages (9,4),(17,4),(33,4) :');
    for j=1:6
        fprintf(' %s %.3f %% |',lab{j},100*(max(V(i3,j,c))-min(V(i3,j,c)))/mean(V(i3,j,c)));
    end
    fprintf('\n');
end

%% ---- GARDE : (17,4) Phi_O=0 reproduit-il M5 / Table 10 ? ----
kp=find(TIL(:,1)==17 & TIL(:,2)==4,1);
PUB=[114.31 19.106 227.57 0.9153 61.02];
nm={'couple','I1','pertes fer','rendement','X_m0'}; ok=true;
fprintf('\n  ---- GARDE : (17,4) Phi_O = 0 contre les valeurs publiees ----\n');
for j=1:5
    e=abs(V(kp,j,1)-PUB(j))/abs(PUB(j)); if e>5e-4, ok=false; end
    fprintf('    %-12s calcule %10.4f | publie %10.4f | ecart %.2e\n',nm{j},V(kp,j,1),PUB(j),e);
end
%  et (33,16) cavites contre Z1 (Table 8) : I0 9.153, T 115.53, I1 18.675, Xm0 70.029
k16=find(TIL(:,1)==33 & TIL(:,2)==16,1);
PUB2=[115.53 18.675 70.029 9.153]; got2=[V(k16,1,2) V(k16,2,2) V(k16,5,2) V(k16,6,2)];
nm2={'couple','I1','X_m0','I0'};
fprintf('  ---- GARDE : (33,16) cavites contre Z1 / Table 8 ----\n');
for j=1:4
    e=abs(got2(j)-PUB2(j))/abs(PUB2(j)); if e>5e-4, ok=false; end
    fprintf('    %-12s calcule %10.4f | publie %10.4f | ecart %.2e\n',nm2{j},got2(j),PUB2(j),e);
end
fprintf('    GARDE %s\n',tern(ok,'PASSEE','ECHOUEE -- chaine differente du manuscrit'));
save('Z6_coupled_conv_cav.mat','TIL','V','NC','TF','I1F','PfF','I0F','ok');
fprintf('\n  duree %.0f s\n=== Z6 termine ===\n',toc(t0));
diary off;

function s=tern(c,a,b), if c, s=a; else, s=b; end, end
