%% RUN_T16_BASE_P1 - base de surface d'ordre superieur
%
%  T16. La base P0 (constante par colonne) projette en 1/n, donc le terme
%  d'energie Kn*|W|^2 ~ n*(1/n)^2 = 1/n : SERIE HARMONIQUE, logarithmiquement
%  divergente. C'est la cause de la derive de X_m et de k_C mesuree en T1/T2.
%
%  La base P1 (chapeau) projette en 1/n^2, donc Kn*|W|^2 ~ 1/n^3 :
%  ABSOLUMENT convergente. On doit donc observer une VRAIE limite en Nh.
%
%  PAVAGE QUASI UNIFORME. Le chapeau symetrique de demi-largeur dth n'est la
%  fonction P1 standard que sur une grille uniforme. On choisit donc nT/nO de
%  sorte que colonnes de face et d'ouverture aient presque la meme largeur :
%  aS/nT ~= oS/nO. Le rapport aS/oS vaut ici ~4.3, d'ou nT=17, nO=4.
clear; clc;
diary('T16_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G;
Rs=G.Rs; Rr=G.Rr; L=M.L;
taus=2*pi/M.Ns; taur=2*pi/M.Nr;
aS=(G.taus-M.bs0)/Rs; oS=taus-aS;
nT=17; nO=4;
fprintf('=== T16 : base de surface P0 vs P1 ===\n');
fprintf('  face/ouverture stator : %.4f / %.4f rad  (rapport %.2f)\n',aS,oS,aS/oS);
fprintf('  pavage nT=%d nO=%d -> colonnes %.5f (face) / %.5f (ouv.) rad\n', ...
    nT,nO,aS/nT,oS/nO);
fprintf('  non-uniformite residuelle : %.1f %%\n\n', ...
    100*abs(aS/nT-oS/nO)/(aS/nT));

thsT=2*pi*(0:M.Ns-1)/M.Ns; thrT=2*pi*(0:M.Nr-1)/M.Nr;
Nlist=[512 1024 2048 4096 8192];
fprintf('  %7s | %10s %10s | %10s %10s\n','N_h','Xm P0','k_C P0','Xm P1','k_C P1');
R=nan(numel(Nlist),4);
for k=1:numel(Nlist)
    Nh=Nlist(k);
    for b=1:2
        bas={'p0','p1'}; bs=bas{b};
        A=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,bs);
        cx=ctx; cx.AG=A; Rd=mec.magnetizing(cx,0.2);
        AL=mec.airgap_fourier(thsT,repmat(taus,1,M.Ns), ...
                              thrT,repmat(taur,1,M.Nr),Rs,Rr,L,Nh,bs);
        cl=ctx; cl.AG=AL; Rl=mec.magnetizing(cl,0.2);
        R(k,2*b-1)=Rd.Xm; R(k,2*b)=Rl.Xm/Rd.Xm;
    end
    fprintf('  %7d | %10.3f %10.4f | %10.3f %10.4f\n',Nh,R(k,1),R(k,2),R(k,3),R(k,4));
end

%% ---- verdict ---------------------------------------------------------
d0=100*(max(R(:,2))-min(R(:,2)))/mean(R(:,2));
d1=100*(max(R(:,4))-min(R(:,4)))/mean(R(:,4));
fprintf('\n  dispersion de k_C : P0 %.2f %%  ->  P1 %.2f %%\n',d0,d1);
fprintf('  Carter classique  : %.4f\n',ctx.AGcarter.kC);
%  Increments successifs : une serie en 1/n^3 doit voir ses increments
%  s'effondrer, une serie en 1/n doit les voir se maintenir.
fprintf('\n  increments de Xm entre doublements de Nh [ohm]\n');
fprintf('  %10s %10s %10s\n','Nh','P0','P1');
for k=2:numel(Nlist)
    fprintf('  %10d %10.4f %10.4f\n',Nlist(k),R(k,1)-R(k-1,1),R(k,3)-R(k-1,3));
end
save('T16_p1.mat','R','Nlist','nT','nO');
fprintf('\n=== T16 termine ===\n');
diary off;
