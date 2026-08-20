%% RUN_RIPPLE  -  Ondulation de couple : entrefer C1 analytique (A1 + A5)
%
%  Couple électromagnétique en charge vs position rotorique, sur le maillage
%  polaire 2D, par TRAVAUX VIRTUELS ANALYTIQUES :
%       T = 1/2 * sum_entrefer (Delta U)^2 * dLambda/dtheta ,
%  avec Lambda(sigma) perméance d'entrefer C1 (noyau gaussien de largeur
%  PHYSIQUE sigma0) et dLambda/dtheta ANALYTIQUE (une seule résolution par
%  position). C'est la voie Gyselinck [29]/Lannoo [34] recommandée par le
%  skill (A1+A5), en remplacement du tenseur de Maxwell monocouche (M6).
%
%  RÉSULTAT CLÉ : la largeur de couplage sigma0 est un paramètre PHYSIQUE
%  (largeur de la perméance dent-à-dent, ~ pas d'encoche, à identifier par
%  FEM). Trop étroit -> ondulation sur-estimée (artefact M6 du réseau
%  discret) ; calibré ~ pas d'encoche -> ondulation PHYSIQUE (~1 %), avec
%  couple moyen et Xm préservés.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== Ondulation de couple (entrefer C1 analytique, A1+A5) ===\n\n');

M=mec.machine_18_5kW(); ctx=mec.build_context(M);
G=ctx.G; W=ctx.W; BH=ctx.BH; p=M.p; Nr=M.Nr;

s=0.02; rr=mec.equivalent_circuit(ctx,s,ctx.Xm0);
fprintf('Point nominal s=%.3f : I1=%.1f A, I2=%.1f A, Tem(circuit)=%.1f N.m\n\n',...
    s, rr.I1, rr.I2, rr.Tem);

psi1=angle(rr.I1c); psi2=angle(rr.I2c);
Ibar_rms=rr.I2*(2*M.m*W.kw1*W.Nph)/Nr;
thbar=2*pi*(0:Nr-1)'/Nr;

nc=6; nr=3;
Npos=37; th0=linspace(0, 2*(2*pi/M.Ns), Npos);
tau=0.5*(G.taus+G.taur);                      % pas d'encoche moyen

% --- Sensibilité à la largeur de couplage physique sigma0 ---
sig_list=[0.5 0.75 1.0]*tau;
Tall=cell(numel(sig_list),1); lab=cell(numel(sig_list),1);
fprintf('%10s | %6s %8s %9s\n','sigma0','Xm','T_moyen','ondul.');
for j=1:numel(sig_list)
    M.opt.gap_sigma0=sig_list(j);
    Xm=noloadXm(M,G,W,BH,nc,nr,p);
    T=zeros(Npos,1);
    for k=1:Npos
        ph=p*th0(k);
        i3k=sqrt(2)*rr.I1*[cos(psi1+ph);cos(psi1+ph-2*pi/3);cos(psi1+ph+2*pi/3)];
        ibk=sqrt(2)*Ibar_rms*cos(p*thbar-psi2);
        me=mec.mesh_refined(M,G,nc,i3k,th0(k),nr,ibk);
        Se=mec.solve_mesh(me,BH,M.opt);
        gi=me.gapfirst:numel(me.a); dU=Se.U(me.a(gi))-Se.U(me.b(gi));
        T(k)=0.5*sum(dU.^2.*me.gap_dLdth);
    end
    Tall{j}=T; ond=100*(max(T)-min(T))/abs(mean(T));
    lab{j}=sprintf('\\sigma_0=%.2f\\tau (ond %.1f%%)',sig_list(j)/tau,ond);
    fprintf('%6.2f*tau | %6.1f %8.1f %8.1f%%\n',sig_list(j)/tau,Xm,mean(T),ond);
end

Tc=Tall{end};                                 % sigma0 = pas d'encoche (calibré)
fprintf('\nCALIBRÉ (sigma0 = pas d''encoche) : couple moyen = %.1f N.m ',mean(Tc));
fprintf('(circuit %.1f, EF %.1f), ondulation = %.1f %%\n',...
    rr.Tem, interp1(mec.ansys_ref().s,mec.ansys_ref().T,s),...
    100*(max(Tc)-min(Tc))/abs(mean(Tc)));

% --- Figure ---
figure('Name','Ondulation de couple','Position',[120 120 820 460]); hold on;
col=lines(numel(sig_list));
for j=1:numel(sig_list)
    plot(th0*180/pi, Tall{j},'-o','Color',col(j,:),'LineWidth',1.4,'DisplayName',lab{j});
end
yline(rr.Tem,'k--','circuit (moyen)');
grid on; xlabel('position rotorique [deg mec]'); ylabel('Couple [N.m]');
legend('Location','best');
title('Ondulation de couple — sensibilité à la largeur d''entrefer C1 (nc=6, nr=3)');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'ondulation_couple.png'));
fprintf('\nFigure enregistree : ondulation_couple.png\n');

% ================= fonction locale =================
function Xm=noloadXm(M,G,W,BH,nc,nr,p)
%  CORRIGE le 6 aout 2026. Ce script codait en dur 8.32 A, valeur de
%  l'ancre R.anchor.I0 de mec.ansys_ref, qui est FAUSSE : la mesure sur
%  transitoire\a vide\Winding Plot 4.tab (2001 points, vitesse etablie
%  1499,942 tr/min) donne 8,4986 A sur la moyenne des trois phases.
%  Voir A8_verif_ancres_out.txt. Script de diagnostic, non publie, mais
%  la valeur ne doit pas survivre.
I0ref=8.4986;
i3=sqrt(2)*I0ref*[1;-0.5;-0.5]; tp=pi*G.Ds/(2*p);
mm=mec.mesh_refined(M,G,nc,i3,0,nr,[]); S=mec.solve_mesh(mm,BH,M.opt);
isg=false(numel(mm.a),1); isg(mm.gapfirst:end)=true;
phig=accumarray(mm.a(isg),S.Phi(isg),[mm.Ms 1]);
Bcol=phig./(mm.dth_s*mm.Rg*M.L);
Bg1=abs((1/pi)*sum(Bcol.*exp(-1j*p*mm.th_s).*mm.dth_s));
Xm=M.w*W.kw1*W.Nph*((2/pi)*Bg1*tp*M.L)/sqrt(2)/I0ref;
end
