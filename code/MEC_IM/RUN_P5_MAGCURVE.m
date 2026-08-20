%% RUN_P5_MAGCURVE  -  P5 : caracteristique de magnetisation + sensibilites de I0
%
%  E1 est exact (+0,2 %) mais I0 est fort de +9,4 % : la PENTE de la
%  caracteristique de magnetisation est en cause, pas le niveau de flux.
%  Un seul point EF a vide (U = 398,2 V) ne peut pas separer les causes.
%  Ce script :
%    1) produit la courbe a vide MEC I0(U), E1(U), Bg1(U), Xm(U) sur
%       0,5-1,1 pu — la table a confronter au BALAYAGE ANSYS propose
%       (protocole imprime en fin de script) ;
%    2) quantifie la SENSIBILITE de I0 aux trois causes candidates :
%       entrefer g (+/-0,01 mm), courbe B(H) (H x (1+/-5 %)), et
%       discretisation (reseau 1 element/dent vs maillage raffine, RUN_MESH) ;
%    3) en deduit ce que chaque cause DEVRAIT valoir pour expliquer +9,4 %.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== P5 : caracteristique de magnetisation a vide ===\n\n');

M0=mec.machine_18_5kW();
I0_EF=8.489; U_EF=398.17;                     % ancre EF (RUN_VALIDATION)

%% 1. Courbe a vide MEC
pu=0.5:0.1:1.1;
I0=zeros(size(pu)); E1=I0; Bg1=I0; Xm=I0;
for k=1:numel(pu)
    Mv=M0; Mv.Uph=pu(k)*M0.Uph;
    ctx=mec.build_context(Mv);
    r=mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
    I0(k)=r.I1; E1(k)=r.E1; Bg1(k)=r.Bg1; Xm(k)=r.Xm;
end
fprintf('--- 1. Courbe a vide MEC (protocole du balayage ANSYS) ---\n');
fprintf('%6s %8s | %8s %8s %8s %8s\n','U/Un','U [V]','I0 [A]','E1 [V]','Bg1 [T]','Xm [ohm]');
for k=1:numel(pu)
    fprintf('%6.2f %8.1f | %8.3f %8.1f %8.3f %8.1f\n',...
        pu(k),pu(k)*M0.Uph,I0(k),E1(k),Bg1(k),Xm(k));
end
fprintf('Ancre EF disponible : U=%.1f V -> I0=%.3f A (MEC %.3f : %+.1f %%)\n\n',...
    U_EF,I0_EF,interp1(pu*M0.Uph,I0,U_EF),...
    100*(interp1(pu*M0.Uph,I0,U_EF)-I0_EF)/I0_EF);

%% 2. Sensibilites de I0 au point 1 pu
ctx0=mec.build_context(M0);
r0=mec.equivalent_circuit(ctx0,1e-4,ctx0.Xm0); I00=r0.I1;
% (a) entrefer
dg=0.01e-3;
Mg=M0; Mg.g=M0.g+dg; ctxg=mec.build_context(Mg);
rg=mec.equivalent_circuit(ctxg,1e-4,ctxg.Xm0);
% (b) courbe B(H) : H -> 1.05*H (tole 5 % plus dure)
ctxh=ctx0; math=ctx0.mat; math.BH(:,1)=1.05*math.BH(:,1); ctxh.BH=mec.bh(math);
rh=mec.equivalent_circuit(ctxh,1e-4,ctxh.Xm0);
fprintf('--- 2. Sensibilites de I0 (U = 1 pu, I0_base = %.3f A) ---\n',I00);
Sg=100*(rg.I1-I00)/I00;
Sh=100*(rh.I1-I00)/I00;
fprintf('entrefer g +%.0f um (%.3f -> %.3f mm) : I0 %+0.1f %%\n',dg*1e6,M0.g*1e3,Mg.g*1e3,Sg);
fprintf('courbe B(H) : H x 1.05               : I0 %+0.1f %%\n',Sh);
fprintf('discretisation (RUN_MESH)            : Xm raffine -7..-8 %% vs EF => I0 ~ +8 %% (mesure)\n\n');

%% 3. Etat post-P4 et protocole
fprintf('--- 3. Etat et protocole ---\n');
fprintf('Le +9,4 %% HISTORIQUE de I0 a ete resolu par P4 (fuite d''extremite\n');
fprintf('identifiee : la chute jXs*I0 supplementaire abaisse E de ~1 %%, donc la\n');
fprintf('saturation, donc I0 sur la partie raide de la courbe) : residu %+.1f %%.\n',...
    100*(interp1(pu*M0.Uph,I0,U_EF)-I0_EF)/I0_EF);
fprintf('Les sensibilites ci-dessus donnent les BARRES D''ERREUR du point :\n');
fprintf('+/-10 um d''entrefer <=> %+.1f %% ; +/-5 %% sur H(B) <=> %+.1f %%.\n',Sg,Sh);
fprintf('=> le balayage U dans ANSYS reste LE test discriminant de la PENTE\n');
fprintf('   (un seul point ne valide pas Xm(U) hors du point nominal) :\n');
fprintf('   PROTOCOLE : essais a vide U = [0.5:0.1:1.1]*690 V, extraire I_rms\n');
fprintf('   et E (tension induite) en regime etabli ; superposer a la table 1.\n');
fprintf('   La cause « g » decale TOUTE la courbe ; la cause « B(H) » ne joue\n');
fprintf('   qu''au-dela du coude (U >= 0.8 pu) ; la discretisation est\n');
fprintf('   independante de U. Trois signatures distinctes.\n');

%% Figure
figure('Name','P5 courbe de magnetisation','Position',[100 100 900 400]);
subplot(1,2,1);
plot(I0,pu*M0.Uph,'b-o','LineWidth',1.4); hold on;
plot(I0_EF,U_EF,'rs','MarkerSize',10,'LineWidth',2);
grid on; xlabel('I_0 [A]'); ylabel('U_{ph} [V]');
legend('MEC','EF (ancre unique)','Location','southeast');
title('Caracteristique a vide U(I_0)');
subplot(1,2,2);
plot(pu,Xm,'b-o','LineWidth',1.4); grid on;
xlabel('U/U_n'); ylabel('X_m [\Omega]'); title('X_m(U) (saturation)');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'magcurve_P5.png'));
fprintf('\nFigure enregistree : magcurve_P5.png\n=== Termine ===\n');
