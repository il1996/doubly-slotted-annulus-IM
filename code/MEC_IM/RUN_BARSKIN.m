%% RUN_BARSKIN  -  Effet de peau d'une barre TRAPÉZOÏDALE (residu fort glissement)
%
%  Mes kR/kX (mec.cage) viennent des formules de Field-Emde, qui supposent une
%  barre RECTANGULAIRE. L'encoche rotorique reelle est TRAPEZOIDALE
%  (br1=5,80 mm en haut -> br2=2,04 mm au fond). On resout ici le probleme de
%  diffusion 1D exact pour b(y) quelconque (mec.bar_skin) :
%    A) VALIDATION : profil rectangulaire -> doit retomber sur Field-Emde.
%    B) Geometrie REELLE -> vrai kR/kX, et effet sur le couple a fort glissement.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== Effet de peau : barre trapezoidale vs rectangulaire ===\n\n');

M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G; W=ctx.W;
fprintf('Encoche rotor : br1(haut)=%.2f mm -> br2(fond)=%.2f mm, hauteur %.2f mm\n',...
    M.br1*1e3, M.br2*1e3, M.hr1*1e3);
fprintf('  => rapport haut/fond = %.2f  (une barre rectangulaire aurait 1,00)\n\n',M.br1/M.br2);

% formules exactes de Field-Emde (reference du modele actuel)
fe = @(xi) deal( xi*(sinh(2*xi)+sin(2*xi))/(cosh(2*xi)-cos(2*xi)), ...
                 (3/(2*xi))*(sinh(2*xi)-sin(2*xi))/(cosh(2*xi)-cos(2*xi)) );

%% A) Validation du solveur 1D : profil rectangulaire vs Field-Emde
fprintf('--- A) Validation : profil RECTANGULAIRE vs Field-Emde exact ---\n');
fprintf('%6s %6s | %8s %8s | %8s %8s | %7s %7s\n',...
    's','xi','kR_1D','kR_FE','kX_1D','kX_FE','eR%','eX%');
for s=[0.05 0.1 0.25 0.5 1.0]
    R=mec.bar_skin(M,G,s,struct('shape','rect'));
    [kRfe,kXfe]=fe(R.xi);
    fprintf('%6.2f %6.2f | %8.4f %8.4f | %8.4f %8.4f | %6.2f %6.2f\n',...
        s,R.xi,R.kR,kRfe,R.kX,kXfe,100*(R.kR-kRfe)/kRfe,100*(R.kX-kXfe)/kXfe);
end

%% B) Geometrie reelle (trapezoidale)
fprintf('\n--- B) Geometrie REELLE (trapezoidale) vs hypothese rectangulaire ---\n');
fprintf('%6s | %8s %8s %7s | %8s %8s\n','s','kR_trap','kR_rect','ecart%','kX_trap','kX_rect');
slips=[0.05 0.1 0.2 0.5 1.0];
kRt=zeros(size(slips)); kXt=kRt; kRr=kRt; kXr=kRt;
for i=1:numel(slips)
    T=mec.bar_skin(M,G,slips(i),struct('shape','trapz'));
    R=mec.bar_skin(M,G,slips(i),struct('shape','rect'));
    kRt(i)=T.kR; kXt(i)=T.kX; kRr(i)=R.kR; kXr(i)=R.kX;
    fprintf('%6.2f | %8.4f %8.4f %6.1f | %8.4f %8.4f\n',...
        slips(i),T.kR,R.kR,100*(T.kR-R.kR)/R.kR,T.kX,R.kX);
end

%% C) Effet sur le couple
fprintf('\n--- C) Effet sur le couple (Rr recalcule avec kR trapezoidal) ---\n');
ref=mec.ansys_ref();
fprintf('%6s | %8s %8s | %8s %8s %7s\n','s','T_rect','T_trap','T_EF','gap_r%','gap_t%');
for i=1:numel(slips)
    s=slips(i);
    r0=mec.equivalent_circuit(ctx,s,ctx.Xm0);            % modele actuel (rect)
    % Rr avec kR trapezoidal : on remplace la part de barre
    Cg=mec.cage(M,G,W,ctx.Lk,s);
    Rr_rect=Cg.Rr;
    kref=(4*M.m/M.Nr)*(W.kw1*W.Nph)^2/Cg.ksq^2;
    Rbar_dc=Cg.Rbar_dc; Rring_eq=Cg.Rr/kref - Cg.Rbar_dc*Cg.kR;
    Rr_trap=kref*(Rbar_dc*kRt(i) + Rring_eq);
    % couple ~ proportionnel a (Rr/s)/((Rs+Rr/s)^2+Xs^2)  -> recalcul direct
    Xs=ctx.Lk.Xs_leak+ctx.Lk.sigd_s*r0.Xm; Xr=Cg.Xr+Cg.sigd_r*r0.Xm;
    Xt=Xs+Xr; Rs=ctx.Rs;
    tf=@(Rr) (M.m*M.p/M.w)*M.Uph^2*(Rr/s)/((Rs+Rr/s)^2+Xt^2);
    Tr=tf(Rr_rect); Tt=tf(Rr_trap);
    Tfem=interp1(ref.s,ref.T,s);
    fprintf('%6.2f | %8.1f %8.1f | %8.1f %6.1f %6.1f\n',...
        s,Tr,Tt,Tfem,100*(Tr-Tfem)/Tfem,100*(Tt-Tfem)/Tfem);
end
fprintf('\n(couples estimes par la formule de Kloss simplifiee : sert a isoler\n');
fprintf(' l''effet de kR seul, pas a remplacer RUN_MEC_IM)\n');
