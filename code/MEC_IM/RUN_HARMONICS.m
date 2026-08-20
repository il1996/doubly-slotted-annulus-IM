%% RUN_HARMONICS  -  Couples parasites asynchrones (modèle multi-harmonique)
%
%  DIAGNOSTIC quantitatif : le modèle fondamental surestime le couple au fort
%  glissement (+14 % a s=0,1 ... +27 % a s=1 vs EF). Hypothese testee ici :
%  les COUPLES PARASITES des harmoniques d'espace du bobinage.
%
%  Modele multi-harmonique classique : chaque harmonique nu voit son propre
%  champ tournant a Omega_s/nu et son propre glissement
%      s_nu = 1 - nu*(1-s)   (direct)      s_nu = 1 + nu*(1-s)   (inverse)
%  et sa branche rotorique. Les branches harmoniques sont EN SERIE dans le
%  circuit statorique (leur somme = la fuite differentielle sigma_delta*Xm que
%  le modele actuel traite comme une REACTANCE PURE, donc sans couple ni perte).
%
%  Couple d'un harmonique :  T_nu = P_ag,nu * nu*p/omega * dir_nu
%  (le facteur nu est decisif : meme une petite puissance harmonique donne un
%  couple non negligeable).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== Couples parasites harmoniques au fort glissement ===\n\n');

M=mec.machine_18_5kW(); ctx=mec.build_context(M);
W=ctx.W; G=ctx.G; ref=mec.ansys_ref();
H=mec.harmonics(M,W,49);

fprintf('Spectre du bobinage (q=%d, pas %d/%d) :\n',W.q,M.yq,M.Ns/(2*M.p));
fprintf('%4s %9s %6s %10s\n','nu','kw_nu','sens','sigma_nu');
for k=1:numel(H.nu)
    s1='direct'; if H.dir(k)<0, s1='inverse'; end
    fprintf('%4d %9.4f %6s %10.2e\n',H.nu(k),H.kw(k),s1,H.sig(k));
end
fprintf('somme sigma_nu = %.4f  (= fuite differentielle du modele actuel)\n\n',sum(H.sig));

%% Couples parasites vs glissement
fprintf('%6s | %8s %8s %8s | %8s %8s %7s\n',...
    's','T_fond','T_paras','T_net','T_MEC','T_FEM','gap');
slips=[0.05 0.10 0.20 0.50 1.00];
for s=slips
    r=mec.equivalent_circuit(ctx,s,ctx.Xm0);
    Xm=r.Xm; I1=r.I1;
    Cg=mec.cage(M,G,W,ctx.Lk,s);
    Tpar=0;
    for k=1:numel(H.nu)
        nu=H.nu(k); dir=H.dir(k);
        Xmn = Xm*H.sig(k);                       % reactance magnetisante harmonique
        if Xmn<=0, continue; end
        if dir>0, snu = 1 - nu*(1-s); else, snu = 1 + nu*(1-s); end
        if abs(snu)<1e-6, continue; end
        % rotor vu par l'harmonique : REPORT HARMONIQUE (kw_nu, anneaux en
        % sin^2(pi*nu*p/Nr)) + effet de peau a la frequence rotorique |snu|*f
        Cn=mec.cage(M,G,W,ctx.Lk,min(abs(snu),20),nu,H.kw(k));
        Rrn=Cn.Rr; Xrn=Cn.Xr;
        % branche harmonique : jXmn // (Rrn/snu + jXrn)
        Zrot=Rrn/snu + 1i*Xrn;
        Zn = 1/(1/(1i*Xmn) + 1/Zrot);
        Pag = M.m*I1^2*real(Zn);                 % puissance d'entrefer harmonique
        Tnu = Pag*nu*M.p/M.w*dir;                % couple (signe = sens)
        Tpar = Tpar + Tnu;
    end
    Tf=r.Tem; Tnet=Tf+Tpar;
    Tfem=interp1(ref.s,ref.T,s);
    fprintf('%6.2f | %8.1f %8.1f %8.1f | %8.1f %8.1f %6.1f%%\n',...
        s, Tf, Tpar, Tnet, Tf, Tfem, 100*(Tf-Tfem)/Tfem);
end

fprintf('\nLecture : T_paras = somme des couples parasites (negatif = freinant).\n');
fprintf('Si |T_paras| ~ l''ecart (T_MEC - T_FEM), les harmoniques expliquent le gap.\n');
