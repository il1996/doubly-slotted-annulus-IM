%% RUN_A1_RIPPLE_BASIS - l'ondulation sous changement de base P0 -> P1
%
%  ITEM A-1 des CONSIGNES_COWORK_07AOUT. Le handoff §1.9 annonce que P1
%  degrade l'ondulation de +2,6 % a +20,9 %. Or ces deux chiffres viennent
%  de DEUX SCRIPTS DIFFERENTS : le 108,2 de RUN_ARTICLE (base P0 par
%  defaut), le 127,5 d'une reconstruction en P1. Les comparer, c'est
%  comparer deux chaines et non deux bases.
%
%  On pilote donc LA MEME chaine -- celle de RUN_ARTICLE:324-347, maillage
%  mesh_refined + gapF.torque -- avec l'option basis, conformement a la
%  regle du dossier : appeler la chaine existante, ne pas la refaire.
%
%  GARDE EXIGEE PAR LES CONSIGNES. L'ondulation est une grandeur de
%  GRADIENT, portee par la variation du potentiel de surface le long du
%  bore, donc de meme nature que B_t. Or B_t s'AMELIORE sous P1
%  (-12,6 -> -9,3 %). Deux grandeurs de gradient qui evoluent en sens
%  OPPOSES signalent que l'une des deux n'est pas ce qu'on croit. On
%  mesure donc les DEUX dans la meme execution, sous la meme base.
%
%  PIEGE N.12, a verifier de quel cote il se situe : comparer une carte
%  BRUTE a une grandeur TRANSITOIRE. Facteur d'attenuation mesure 0,96.
%  Les deux sont rapportees ici, jamais melangees.
clear; clc; t0=tic;
diary('A1_ripple_basis_out.txt'); diary on;
Npos=19; nT=17; nO=4; Nh=8192;
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
Dch=fullfile(ROOT,'transitoire','en charge');
rdt=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
trq=rdt(fullfile(Dch,'Plot 1.tab'));
fc=trq(:,1)>=1.90;
TppF=max(trq(fc,3))-min(trq(fc,3));      % reference CALCULEE
Mt=mec.machine_18_5kW(); n1=1469.772776;
s_ch=(Mt.ns*60-n1)/(Mt.ns*60);

fprintf('=== A-1 : ondulation sous changement de base ===\n');
fprintf('  MEME chaine (mesh_refined + gapF.torque), seule la base varie\n');
fprintf('  %d positions sur 2 pas d''encoche | s = %.6f\n',Npos,s_ch);
fprintf('  charge MESUREE, fenetre t >= 1.90 (RUN_ARTICLE:354)\n');
fprintf('  reference EF CALCULEE de Plot 1.tab : %.4f N.m c-c\n\n',TppF);

bas={'p0','p1'};
lbl={'base P0 (defaut publie)','base P1 (nT=17 nO=4 Nh=8192)'};
Tpp=nan(1,2); Tdq=nan(1,2); Btl=nan(1,2); Bg1=nan(1,2); Xm0=nan(1,2); Tsc=nan(1,2);
for k=1:2
    M=mec.machine_18_5kW();
    ctx=mec.build_context(M); G=ctx.G; W=ctx.W; BH=ctx.BH;
    if k==2
        A=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
        ctx.AG=A; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
    end
    Xm0(k)=ctx.Xm0;
    r=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0); Tsc(k)=r.Tem;

    %  ---- carte d'ondulation : chaine de RUN_ARTICLE, inchangee ----
    Pth=2*(2*pi/M.Ns); th0=linspace(0,Pth,Npos); Tk=nan(1,Npos);
    for q=1:Npos
        [i3k,ibk]=loc_currents(M,W,r);
        me=mec.mesh_refined(M,G,6,i3k,th0(q),3,ibk);
        Se=mec.solve_mesh(me,BH,M.opt);
        Us=Se.U(me.gapF.ids(1:me.Ms)); Ur=Se.U(me.gapF.ids(me.Ms+1:end));
        Tk(q)=me.gapF.torque(Us,Ur);
    end
    Tpp(k)=max(Tk)-min(Tk);

    %  ---- injection mecanique, charge MESUREE ----
    DT=@(th) interp1(th0,Tk-mean(Tk),mod(th,Pth),'linear');
    oc=mec.dq_startup(ctx,struct('tend',2.0,'TL',[trq(:,1),trq(:,2)],'Trip',DT));
    wc=oc.t>=1.90;
    Tdq(k)=max(oc.Tem(wc))-min(oc.Tem(wc));

    %  ---- GARDE : B_t, grandeur de gradient de MEME nature ----
    thq=linspace(0,2*pi,2001); thq(end)=[]; Rm=0.5*(G.Rs+G.Rr);
    p1=angle(r.I1c);
    i3=sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
    Fs=W.slotMMF(i3);
    p2=angle(r.I2c); tb=2*pi*(0:M.Nr-1)/M.Nr;
    Fu=cumsum(cos(M.p*tb-p2).'); Fu=Fu-mean(Fu);
    c1=(2/M.Nr)*sum(Fu.'.*exp(-1j*M.p*tb));
    Fr=-Fu*((3/2)*(4/pi)*(W.kw1*W.Nph/(2*M.p))*sqrt(2)*abs(r.I2c)/abs(c1));
    S=mec.solve_network(ctx.net,G,BH,ctx.AG,Fs,Fr,M.opt);
    if isfield(ctx.AG,'expand')
        [Br,Bt]=ctx.AG.field(S.Usurf,Rm,thq);
        Br=Br(:).'; Bt=Bt(:).';
        Btl(k)=sqrt(mean(Bt.^2));
        Bg1(k)=abs((2/numel(thq))*sum(Br.*exp(-1j*M.p*thq)));
    end
    fprintf('  %-30s X_m0 %7.3f | carte %8.3f | post-dq %8.3f | B_t %7.4f\n', ...
        lbl{k},Xm0(k),Tpp(k),Tdq(k),Btl(k));
end

%% ---- LE TABLEAU QUI TRANCHE -----------------------------------------
fprintf('\n  ---- ondulation : carte BRUTE et grandeur TRANSITOIRE ----\n');
fprintf('  %-30s %11s %11s %11s %11s\n', ...
    'base','carte (N.m)','post-dq','vs EF','attenuation');
for k=1:2
    fprintf('  %-30s %11.3f %11.3f %10.2f %% %11.4f\n', ...
        lbl{k},Tpp(k),Tdq(k),100*(Tdq(k)-TppF)/TppF,Tdq(k)/Tpp(k));
end
fprintf('\n  effet du changement de base P0 -> P1 :\n');
fprintf('    sur la carte brute : %+7.2f %%\n',100*(Tpp(2)-Tpp(1))/Tpp(1));
fprintf('    sur le post-dq     : %+7.2f %%\n',100*(Tdq(2)-Tdq(1))/Tdq(1));

%% ---- LA GARDE : deux grandeurs de gradient -------------------------
fprintf('\n  ---- GARDE : B_t et ondulation sont-elles de meme nature ? ----\n');
BtF=0.131;
fprintf('  %-30s %11s %11s %11s\n','grandeur','P0','P1','evolution');
fprintf('  %-30s %11.4f %11.4f %10.2f %%\n','B_t rms (T)',Btl(1),Btl(2), ...
    100*(Btl(2)-Btl(1))/Btl(1));
fprintf('  %-30s %11.2f %11.2f %10.2f %%\n','ecart de B_t a l''EF (%)', ...
    100*(Btl(1)-BtF)/BtF,100*(Btl(2)-BtF)/BtF, ...
    100*(Btl(2)-BtF)/BtF-100*(Btl(1)-BtF)/BtF);
fprintf('  %-30s %11.2f %11.2f %10.2f %%\n','ecart d''ondulation a l''EF (%)', ...
    100*(Tdq(1)-TppF)/TppF,100*(Tdq(2)-TppF)/TppF, ...
    100*(Tdq(2)-TppF)/TppF-100*(Tdq(1)-TppF)/TppF);
fprintf('  %-30s %11.3f %11.3f %10.2f %%\n','X_m0 (ohm)',Xm0(1),Xm0(2), ...
    100*(Xm0(2)-Xm0(1))/Xm0(1));
fprintf('  %-30s %11.3f %11.3f %10.2f %%\n','couple du schema (N.m)', ...
    Tsc(1),Tsc(2),100*(Tsc(2)-Tsc(1))/Tsc(1));

fprintf(['\n  LECTURE. Si l''ondulation suit le COUPLE (donc les courants) et\n' ...
         '  non B_t, elle n''est pas une grandeur de gradient du potentiel de\n' ...
         '  surface : elle est proportionnelle a l''excitation. Sa degradation\n' ...
         '  sous P1 serait alors le report de la hausse de X_m0 sur les\n' ...
         '  courants, non un defaut de la base sur les gradients.\n']);
save('A1_ripple_basis.mat','Tpp','Tdq','Btl','Bg1','Xm0','Tsc','TppF');
fprintf('\n  duree %.0f s\n=== A-1 termine ===\n',toc(t0));
diary off;

% ======================================================================
function [i3,ib]=loc_currents(M,W,r)
    psi1=angle(r.I1c); psi2=angle(r.I2c);
    i3=sqrt(2)*r.I1*[cos(psi1);cos(psi1-2*pi/3);cos(psi1+2*pi/3)];
    Ibar=r.I2*(2*M.m*W.kw1*W.Nph)/M.Nr;
    ib=sqrt(2)*Ibar*cos(M.p*(2*pi*(0:M.Nr-1)'/M.Nr)-psi2);
end
