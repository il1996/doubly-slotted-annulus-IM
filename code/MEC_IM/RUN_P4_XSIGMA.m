%% RUN_P4_XSIGMA  -  P4 : identification conjointe {encoches FEM + extremite 3D}
%
%  DOCTRINE (etablie par RUN_SLOTLEAK/RUN_LEAKAGE) : corriger les permeances
%  d'encoche vers leurs valeurs FEM SANS identifier la fuite d'extremite 3D
%  degrade le modele (erreurs compensatoires). Ce script fait les DEUX :
%
%   A) les permeances d'encoche FEM sont desormais DANS mec.leakage
%      (stator 2,722 ; rotor corps 1,043 + isthme 0,500, terme empirique
%      supprime) — rappel et effet net sur Xsigma ;
%   B) le residu d'extremite 3D (tetes reelles, environnement des anneaux,
%      zig-zag 3D) est identifie sur les Xsigma IMPLIQUES par l'EF :
%         Rr/s|EF = T*w/(m*p*I^2)   ;   Xsig|EF = sqrt((U/I)^2-(Rs+Rr/s)^2)
%      aux glissements 0,2 / 0,5 / 1,0 (a fort glissement I2 ~ I1) ;
%   C) TEST DE COHERENCE : le delta doit etre ~CONSTANT en s (une fuite
%      d'extremite est un chemin d'air, sans effet de peau) ;
%   D) point fixe (X_end change I1 qui change l'implication) puis tableau
%      final T(s)/I(s).
%
%  La valeur identifiee est portee par Lk.Xs_end3D (mec.leakage),
%  surchargable par M.opt.Xs_end3D.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== P4 : identification conjointe encoches FEM + extremite 3D ===\n\n');

M=mec.machine_18_5kW(); ref=mec.ansys_ref();
slips=[0.2 0.5 1.0];

%% A) rappel des identifications FEM (constantes de mec.leakage)
G=mec.geometry(M); W=mec.winding(M); Lk=mec.leakage(M,G,W);
fprintf('--- A) Encoches FEM (RUN_SLOTLEAK, dans mec.leakage) ---\n');
fprintf('stator : lam_slot = %.3f (FEM 2.722 a k1=1 ; analytique 2.143)\n',...
    Lk.Xs_slot/(M.w*(4*M.m/M.Ns)*4e-7*pi*M.L*W.Nph^2));
fprintf('rotor  : corps %.3f + isthme+zigzag %.3f  (ancien : 1.427+%.3f+empirique 1.161)\n\n',...
    Lk.lam_r_slot, Lk.lam_r_tip, M.hr0/M.br0);

%% B-D) identification de Xs_end3D par point fixe
Xid=0;
for it=1:4
    Mi=M; Mi.opt.Xs_end3D=Xid;
    ctx=mec.build_context(Mi);
    dX=zeros(size(slips)); Xef=dX; Xme=dX;
    for q=1:numel(slips)
        s=slips(q);
        Ife=interp1(ref.s2,ref.I,s); Tfe=interp1(ref.s,ref.T,s);
        RrS=Tfe*M.w/(M.m*M.p*Ife^2);
        Xef(q)=sqrt(max((M.Uph/Ife)^2-(ctx.Rs+RrS)^2,0));
        r=mec.equivalent_circuit(ctx,s,ctx.Xm0);
        Cg=mec.cage(Mi,ctx.G,ctx.W,ctx.Lk,s);
        Xme(q)=sqrt(max((M.Uph/r.I1)^2-(ctx.Rs+Cg.Rr/s)^2,0));
        dX(q)=Xef(q)-Xme(q);
    end
    if it==1
        fprintf('--- B) Xsigma implique (encoches FEM, X_end3D = 0) ---\n');
        fprintf('%6s | %8s %8s %8s\n','s','X_EF','X_MEC','delta');
        for q=1:numel(slips)
            fprintf('%6.2f | %8.3f %8.3f %+8.3f\n',slips(q),Xef(q),Xme(q),dX(q));
        end
        fprintf('--- C) constance du delta : %.3f / %.3f / %.3f ohm',dX);
        fprintf('  (ecart-type %.3f = %.0f %% de la moyenne)\n\n',std(dX),100*std(dX)/abs(mean(dX)));
    end
    Xid=Xid+mean(dX);
    fprintf('   iteration %d : X_end3D = %.3f ohm (delta residuel moyen %+.3f)\n',it,Xid,mean(dX));
end

fprintf('\n=> Xs_end3D IDENTIFIE = %.3f ohm  (a reporter dans mec.leakage)\n\n',Xid);

%% D) validation finale avec la valeur identifiee
Mi=M; Mi.opt.Xs_end3D=Xid;
ctx=mec.build_context(Mi);
fprintf('--- D) T(s)/I(s) avec encoches FEM + X_end3D = %.3f ohm ---\n',Xid);
fprintf('%6s | %8s %8s %7s | %8s %8s %7s\n','s','T_MEC','T_EF','err%','I_MEC','I_EF','err%');
for s=[0.005 0.02 0.05 0.10 0.20 0.50 1.0]
    r=mec.equivalent_circuit(ctx,s,ctx.Xm0);
    Tf=interp1(ref.s,ref.T,s); If=interp1(ref.s2,ref.I,s);
    fprintf('%6.3f | %8.1f %8.1f %6.1f | %8.1f %8.1f %6.1f\n',...
        s,r.Tem,Tf,100*(r.Tem-Tf)/Tf,r.I1,If,100*(r.I1-If)/If);
end
c=mec.power_balance(mec.equivalent_circuit(ctx,0.02,ctx.Xm0),Mi);
fprintf('\nBilan B6 au nominal : %.2e\n=== Termine ===\n',c.err_global);
