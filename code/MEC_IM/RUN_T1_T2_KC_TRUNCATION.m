%% RUN_T1_T2 - stabilite de k_C et convergence de X_m en troncature
%
%  T1. k_C(N_h) a PAVAGE FIXE (384+352 colonnes, nT=6 nO=2), seul N_h varie.
%      k_C = X_m(surfaces lisses) / X_m(geometrie encochee).
%      L'argument de compensation numerateur/denominateur de la §6.4 est
%      suspect : le cas lisse n'a AUCUN saut de potentiel de surface, le cas
%      encoche en a un par dent. Les deux queues harmoniques n'ont donc
%      aucune raison de se compenser. Ce script le tranche par la mesure.
%
%  T2. X_m au-dela du regime de couplage. 1/X_g = 337 : tout N_h < 337 est
%      SOUS la resolution du couplage, aucun harmonique retenu n'etant dans
%      le regime ou 1/sinh(n*X_g) se distingue du quasi-uniforme. On reprend
%      donc le balayage a partir de 337 et on regarde X_m vs ln(N_h) : si
%      l'alignement est lineaire, la queue est LOGARITHMIQUE et
%      l'extrapolation geometrique de la §6.5 est invalide.
clear; clc;
diary('T1_T2_out.txt'); diary on;
t0=tic;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G;
Rs=G.Rs; Rr=G.Rr; L=M.L; Xg=log(Rs/Rr);
nT=6; nO=2;
thsT=2*pi*(0:M.Ns-1)/M.Ns; thrT=2*pi*(0:M.Nr-1)/M.Nr;
taus=2*pi/M.Ns; taur=2*pi/M.Nr;

fprintf('=== T1/T2 : troncature harmonique ===\n');
fprintf('  X_g = ln(Rs/Rr) = %.4e   1/X_g = %.1f\n',Xg,1/Xg);
fprintf('  pavage FIXE nT=%d nO=%d  (%d + %d colonnes)\n\n', ...
    nT,nO,M.Ns*(nT+nO),M.Nr*(nT+nO));

Nlist = unique([512 1024 2048 3088 6176 337 674 1348 2696 5392 10784]);
res = nan(numel(Nlist),4);   % Nh | Xm_dente | Xm_lisse | kC

fprintf('  %7s %12s %12s %10s %10s\n', ...
    'N_h','Xm dente','Xm lisse','k_C','t [s]');
for k=1:numel(Nlist)
    Nh=Nlist(k); tk=tic;
    %  --- geometrie encochee : operateur condense, pavage fixe ---
    A=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh);
    cd_=ctx; cd_.AG=A; Rd=mec.magnetizing(cd_,0.2);
    %  --- surfaces LISSES : arcs pavant 100 % de l'alesage, meme N_h ---
    AL=mec.airgap_fourier(thsT,repmat(taus,1,M.Ns), ...
                          thrT,repmat(taur,1,M.Nr),Rs,Rr,L,Nh);
    cl_=ctx; cl_.AG=AL; Rl=mec.magnetizing(cl_,0.2);
    res(k,:)=[Nh, Rd.Xm, Rl.Xm, Rl.Xm/Rd.Xm];
    fprintf('  %7d %12.3f %12.3f %10.4f %10.1f\n', ...
        Nh,Rd.Xm,Rl.Xm,Rl.Xm/Rd.Xm,toc(tk));
end

%% ---- T1 : verdict sur la stabilite de k_C ---------------------------
sel=ismember(res(:,1),[512 1024 2048 3088 6176]);
kc=res(sel,4);
dk=100*(max(kc)-min(kc))/mean(kc);
fprintf('\nT1  k_C sur [512 1024 2048 3088 6176]\n');
fprintf('    min %.4f  max %.4f  moyenne %.4f\n',min(kc),max(kc),mean(kc));
fprintf('    dispersion %.2f %%   [Carter classique %.4f]\n', ...
    dk,ctx.AGcarter.kC);
if dk<0.5
    fprintf('    STABLE (<0,5 %%) : le resultat k_C emergent tient.\n');
else
    fprintf('    DERIVE (>=0,5 %%) : a declarer AVANT le rapporteur.\n');
end

%% ---- T2 : nature de la queue ----------------------------------------
s2=res(:,1)>=337;
x=log(res(s2,1)); y=res(s2,2);
p1=polyfit(x,y,1); yfit=polyval(p1,x);
R2=1-sum((y-yfit).^2)/sum((y-mean(y)).^2);
fprintf('\nT2  X_m vs ln(N_h), N_h >= 1/X_g = 337\n');
fprintf('    pente %.4f ohm par e-fold   ordonnee %.3f\n',p1(1),p1(2));
fprintf('    R^2 = %.5f\n',R2);
if R2>0.99
    fprintf(['    ALIGNEMENT LINEAIRE : la queue est LOGARITHMIQUE, donc\n' ...
             '    X_m ne converge pas et l''extrapolation geometrique de la\n' ...
             '    §6.5 est INVALIDE. A chiffrer dans le manuscrit.\n']);
else
    fprintf('    pas d''alignement logarithmique net (R^2 = %.4f).\n',R2);
end
%  Rapports successifs : une queue geometrique donne un rapport constant.
d=diff(y);
fprintf('\n    increments successifs de X_m [ohm] et leurs rapports :\n');
for k=1:numel(d)
    if k<numel(d)
        fprintf('      %7d -> %7d : %+6.3f   rapport suivant %.3f\n', ...
            res(find(s2,1)+k-1,1),res(find(s2,1)+k,1),d(k),d(k+1)/d(k));
    else
        fprintf('      %7d -> %7d : %+6.3f\n', ...
            res(find(s2,1)+k-1,1),res(find(s2,1)+k,1),d(k));
    end
end

save('T1_T2_kc.mat','res');
fprintf('\n  duree totale %.0f s\n=== T1/T2 termine ===\n',toc(t0));
diary off;
