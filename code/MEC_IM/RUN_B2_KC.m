%% RUN_B2_KC - tableau k_C(N_h) publiable, P0 et P1 cote a cote
%
%  B2. Version publiable de T1/T16 : MEMES N_h pour les deux bases, colonnes
%  X_m lisse / X_m encoche / k_C pour chacune, plus increments et rapports.
%
%  POINT A FAIRE RESSORTIR. En P0, X_m LISSE converge (+0,19 %) pendant que
%  X_m ENCOCHE derive (+9,3 %). C'est ce contraste en deux colonnes qui
%  demontre que l'argument de compensation numerateur/denominateur de la
%  §6.4 est FAUX : le cas lisse n'a AUCUN saut de potentiel de surface, le
%  cas encoche en a un par dent, et les deux queues harmoniques n'ont donc
%  aucune raison de se compenser dans un rapport.
clear; clc; t0=tic;
diary('B2_kC_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G;
Rs=G.Rs; Rr=G.Rr; L=M.L; Xg=log(Rs/Rr);
taus=2*pi/M.Ns; taur=2*pi/M.Nr;
thsT=2*pi*(0:M.Ns-1)/M.Ns; thrT=2*pi*(0:M.Nr-1)/M.Nr;
%  Pavage quasi uniforme : condition pour que le chapeau symetrique soit la
%  fonction P1 standard (cf. T16). Rapport face/ouverture ~ 4,36.
nT=17; nO=4;

fprintf('=== B2 : k_C(N_h), bases P0 et P1 ===\n');
fprintf('  MAS 48/44 18,5 kW | X_g = %.4e | 1/X_g = %.1f\n',Xg,1/Xg);
fprintf('  pavage FIXE nT=%d nO=%d | Carter classique = %.4f\n\n',nT,nO,ctx.AGcarter.kC);

Nl=[512 1024 2048 4096 8192];
R=nan(numel(Nl),6);   % Xm_lisse_P0 Xm_dent_P0 kC_P0 | idem P1
for k=1:numel(Nl)
    Nh=Nl(k);
    for b=1:2
        bas={'p0','p1'}; bs=bas{b};
        A=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,bs);
        cd_=ctx; cd_.AG=A; Rd=mec.magnetizing(cd_,0.2);
        AL=mec.airgap_fourier(thsT,repmat(taus,1,M.Ns), ...
                              thrT,repmat(taur,1,M.Nr),Rs,Rr,L,Nh,bs);
        cl_=ctx; cl_.AG=AL; Rl=mec.magnetizing(cl_,0.2);
        R(k,3*b-2:3*b)=[Rl.Xm, Rd.Xm, Rl.Xm/Rd.Xm];
    end
end

fprintf('  %6s | %11s %11s %9s | %11s %11s %9s\n', ...
    'N_h','Xm lisse P0','Xm dent P0','k_C P0','Xm lisse P1','Xm dent P1','k_C P1');
for k=1:numel(Nl)
    fprintf('  %6d | %11.3f %11.3f %9.4f | %11.3f %11.3f %9.4f\n',Nl(k),R(k,:));
end

%% ---- LE CONTRASTE QUI DEMONTRE LA §6.4 ------------------------------
fprintf('\n  --- derive de chaque colonne sur le balayage ---\n');
nm={'X_m LISSE  P0','X_m ENCOCHE P0','k_C P0','X_m LISSE  P1','X_m ENCOCHE P1','k_C P1'};
for j=1:6
    fprintf('  %-16s %9.3f -> %9.3f   %+7.2f %%\n', ...
        nm{j},R(1,j),R(end,j),100*(R(end,j)-R(1,j))/R(1,j));
end
fprintf(['\n  EN P0 : le numerateur (lisse) ne bouge pas, le denominateur\n' ...
         '  (encoche) derive. Le rapport k_C herite INTEGRALEMENT de la\n' ...
         '  derive du denominateur. L''argument de compensation de la §6.4\n' ...
         '  est donc FAUX, et ce tableau le montre en deux colonnes.\n']);

%% ---- increments et rapports ------------------------------------------
fprintf('\n  --- increments de X_m encoche entre doublements, et rapports ---\n');
fprintf('  %8s %13s %9s %13s %9s\n','N_h','dP0','rap.P0','dP1','rap.P1');
d0=diff(R(:,2)); d1=diff(R(:,5));
for k=1:numel(d0)
    r0=NaN; r1=NaN;
    if k<numel(d0), r0=d0(k+1)/d0(k); r1=d1(k+1)/d1(k); end
    fprintf('  %8d %13.4f %9.3f %13.4f %9.3f\n',Nl(k+1),d0(k),r0,d1(k),r1);
end
fprintf(['  Rapports P0 -> 1 : queue LOGARITHMIQUE (T17 le demontre en\n' ...
         '  forme fermee). Rapports P1 ~0,3 : convergence geometrique.\n']);

d=@(x)100*(max(x)-min(x))/mean(x);
fprintf('\n  dispersion de k_C : P0 %.2f %%  ->  P1 %.2f %%\n',d(R(:,3)),d(R(:,6)));
fprintf('  valeur convergee P1 : k_C = %.4f  (Carter %.4f, ecart %+.1f %%)\n', ...
    R(end,6),ctx.AGcarter.kC,100*(R(end,6)-ctx.AGcarter.kC)/ctx.AGcarter.kC);
fprintf('  X_m encoche convergee P1 : %.3f ohm\n',R(end,5));
save('B2_kC.mat','R','Nl','nT','nO');
fprintf('\n  duree %.0f s\n=== B2 termine ===\n',toc(t0));
diary off;
