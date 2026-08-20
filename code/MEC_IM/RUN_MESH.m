%% RUN_MESH  -  Convergence du maillage 2D raffinable (A2/A3)
%
%  Etude de convergence de la reactance magnetisante Xm a vide sur le
%  maillage polaire CONFORME raffinable (mec.mesh_refined + mec.solve_mesh) :
%    * RADIALE     : nc fixe, nr = 1..5 (couches par region) ;
%    * TANGENTIELLE: nr fixe, nc = 1..8 (colonnes par dent ET par encoche).
%
%  CRITERE (A1/A2, invariance au maillage) : l'entrefer etant parametre en
%  grandeurs PHYSIQUES (Carter + noyau gaussien de largeur sigma0 identifiee
%  par FEM), Xm doit CONVERGER (<1 % entre les deux derniers niveaux) et non
%  deriver avec la discretisation — c'est la correction du defaut M1 de
%  Silva (ou gamma_m/gamma_z sont des arcs d'ELEMENTS).
%
%  La valeur convergee est comparee a l'EF (ANSYS ~46 ohm).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== Convergence du maillage 2D raffinable (A2/A3) ===\n\n');

M=mec.machine_18_5kW(); ctx=mec.build_context(M);
G=ctx.G; W=ctx.W; BH=ctx.BH; p=M.p;
I0=8.32;                                     % courant magnetisant d'essai [A rms]
XmEF=mec.ansys_ref().anchor.Xm;

%% Convergence RADIALE (nc fixe)
nc0=5; nr_list=1:7;
Xr=zeros(numel(nr_list),1);
fprintf('--- RADIALE (nc=%d) ---\n',nc0);
fprintf('%6s %10s %10s\n','nr','Xm [ohm]','var %');
for k=1:numel(nr_list)
    Xr(k)=noloadXm(M,G,W,BH,nc0,nr_list(k),p,I0);
    if k==1, v=''; else, v=sprintf('%+.2f',100*(Xr(k)-Xr(k-1))/Xr(k-1)); end
    fprintf('%6d %10.2f %10s\n',nr_list(k),Xr(k),v);
end
vr=100*(Xr(end)-Xr(end-1))/Xr(end-1);

%% Convergence TANGENTIELLE (nr fixe)
nr0=3; nc_list=1:8;
Xt=zeros(numel(nc_list),1);
fprintf('\n--- TANGENTIELLE (nr=%d) ---\n',nr0);
fprintf('%6s %10s %10s\n','nc','Xm [ohm]','var %');
for k=1:numel(nc_list)
    Xt(k)=noloadXm(M,G,W,BH,nc_list(k),nr0,p,I0);
    if k==1, v=''; else, v=sprintf('%+.2f',100*(Xt(k)-Xt(k-1))/Xt(k-1)); end
    fprintf('%6d %10.2f %10s\n',nc_list(k),Xt(k),v);
end
vt=100*(Xt(end)-Xt(end-1))/Xt(end-1);

%% Bilan
fprintf('\n--- BILAN ---\n');
fprintf('Xm (nc=%d, nr=%d) = %.2f ohm ; Xm (nc=%d, nr=%d) = %.2f ohm   (EF ~%.0f)\n',...
    nc0,nr_list(end),Xr(end),nc_list(end),nr0,Xt(end),XmEF);
fprintf('ecart vs EF : %+.1f %% (radial) / %+.1f %% (tangentiel)\n',...
    100*(Xr(end)-XmEF)/XmEF,100*(Xt(end)-XmEF)/XmEF);
fprintf('variation radiale     nr %d->%d : %+.2f %%\n',nr_list(end-1),nr_list(end),vr);
fprintf('variation tangentielle nc %d->%d : %+.2f %%\n',nc_list(end-1),nc_list(end),vt);
ok = abs(vr)<1 && abs(vt)<1;
if ok, fprintf('=> CONVERGE (<1 %% dans les deux directions) : critere A2 satisfait.\n');
else,  fprintf('=> NON CONVERGE : raffiner davantage.\n'); end

%% Figure
figure('Name','Convergence maillage','Position',[100 100 900 380]);
subplot(1,2,1);
plot(nr_list,Xr,'b-o','LineWidth',1.4); hold on; yline(XmEF,'r--','EF');
grid on; xlabel('nr (couches radiales/region)'); ylabel('X_m [\Omega]');
title(sprintf('Convergence RADIALE (nc=%d)',nc0));
subplot(1,2,2);
plot(nc_list,Xt,'b-o','LineWidth',1.4); hold on; yline(XmEF,'r--','EF');
grid on; xlabel('nc (colonnes/dent et /encoche)'); ylabel('X_m [\Omega]');
title(sprintf('Convergence TANGENTIELLE (nr=%d)',nr0));
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'convergence_maillage.png'));
fprintf('\nFigure enregistree : convergence_maillage.png\n');

% ================= fonction locale =================
function Xm=noloadXm(M,G,W,BH,nc,nr,p,I0)
%  Xm a vide sur le maillage (nc,nr) : champ a l'instant phase A maximale,
%  fondamental spatial de l'induction d'entrefer -> E1/I0.
i3=sqrt(2)*I0*[1;-0.5;-0.5]; tp=pi*G.Ds/(2*p);
mm=mec.mesh_refined(M,G,nc,i3,0,nr,[]); S=mec.solve_mesh(mm,BH,M.opt);
isg=false(numel(mm.a),1); isg(mm.gapfirst:end)=true;
phig=accumarray(mm.a(isg),S.Phi(isg),[mm.Ms 1]);
Bcol=phig./(mm.dth_s*mm.Rg*M.L);
Bg1=abs((1/pi)*sum(Bcol.*exp(-1j*p*mm.th_s).*mm.dth_s));
Xm=M.w*W.kw1*W.Nph*((2/pi)*Bg1*tp*M.L)/sqrt(2)/I0;
end
