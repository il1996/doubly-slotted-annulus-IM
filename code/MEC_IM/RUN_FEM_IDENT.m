%% RUN_FEM_IDENT  -  Identification FEM de la largeur d'entrefer sigma0 (A1)
%
%  La permeance d'entrefer C1 du maillage raffine (mec.mesh_refined) est un
%  noyau gaussien de largeur PHYSIQUE sigma0 — le parametre qui gouverne
%  l'ondulation de couple (RUN_RIPPLE). Plutot que de le supposer, on
%  l'IDENTIFIE par elements finis (PDE Toolbox) sur le banc canonique a
%  deux dents de Gyselinck [29]/Silva :
%    * banc deroule : dent stator (V=1 via culasse) / entrefer g /
%      dent rotor (V=0 via culasse), analogie electrostatique du potentiel
%      scalaire magnetique ;
%    * balayage du decalage tangentiel sigma, extraction de Lambda(sigma) ;
%    * ajustement gaussien  Lambda = Lfloor + Lmax*exp(-(sigma/sigma0)^2).
%
%  Le sigma0 identifie est le DEFAUT de mesh_refined (M.opt.gap_sigma0).
%  Resultat attendu : sigma0 ~ 0,54*tau (fit) / ~0,56*tau (largeur 1/e),
%  coherent avec le recouvrement des becs.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== Identification FEM de sigma0 (banc 2 dents, A1) ===\n\n');

M=mec.machine_18_5kW(); G=mec.geometry(M);

id=mec.fem_airgap_ident(M,G);

fprintf('Banc : entrefer g=%.3f mm, dents w_s=%.2f / w_r=%.2f mm, tau=%.3f mm\n',...
    M.g*1e3,(G.taus-M.bs0)*1e3,(G.taur-M.br0)*1e3,id.tau*1e3);
fprintf('Maillage FEM : %d elements\n\n',id.nelem);

fprintf('%12s %14s\n','sigma/tau','Lambda [H/m]');
for q=1:numel(id.sigma)
    fprintf('%12.3f %14.4e\n',id.sigma(q)/id.tau,id.Lambda(q));
end

fprintf('\n--- Ajustement gaussien ---\n');
fprintf('sigma0 (fit)      = %.3f*tau  (= %.3f mm)\n',id.ratio,id.sigma0*1e3);
fprintf('sigma0 (1/e)      = %.3f*tau\n',id.ratio_1e);
fprintf('Lmax  (couplage)  = %.4e H/m\n',id.Lmax);
fprintf('Lfloor (plancher) = %.4e H/m  (= %.1f %% de Lmax)\n',...
    id.Lfloor,100*id.Lfloor/id.Lmax);
fprintf('\n=> M.opt.gap_sigma0 (defaut de mesh_refined) = sigma0 identifie.\n');

%% Figure
sfit=linspace(0,max(id.sigma),200);
Lfit=id.Lfloor+id.Lmax*exp(-(sfit/id.sigma0).^2);
figure('Name','Identification sigma0','Position',[120 120 640 420]);
plot(id.sigma/id.tau,id.Lambda,'bo','LineWidth',1.4,'DisplayName','FEM'); hold on;
plot(sfit/id.tau,Lfit,'r-','LineWidth',1.4,'DisplayName',...
    sprintf('fit gaussien (\\sigma_0=%.2f\\tau)',id.ratio));
grid on; xlabel('\sigma / \tau'); ylabel('\Lambda [H/m]');
legend('Location','best');
title('Permeance d''entrefer dent-a-dent \Lambda(\sigma) — banc FEM 2 dents');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'fem_ident_sigma0.png'));
fprintf('\nFigure enregistree : fem_ident_sigma0.png\n');
