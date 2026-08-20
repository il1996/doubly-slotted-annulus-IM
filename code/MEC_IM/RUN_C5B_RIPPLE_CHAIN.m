%% RUN_C5B_RIPPLE_CHAIN - calibration C5 sur LA CHAINE QUI PRODUIT 101,4 N.m
%
%  LECON DE T18, QUE J'AI REENFREINTE EN C5 : piloter la chaine existante,
%  ne pas la refaire. RUN_C5_RIPPLE_SKEWOFF reconstruisait la carte sur le
%  reseau dentaire et donnait 146-158 N.m ; la valeur publiee (101,4 N.m)
%  vient du MAILLAGE, RUN_ARTICLE.m:324-347 :
%      mesh_refined(M,G,6,i3,theta,3,ib) + solve_mesh + me.gapF.torque
%      19 positions sur DEUX pas d'encoche statoriques.
%  On rejoue CETTE chaine, a l'identique, avec et sans vrillage.
%
%  OU LE VRILLAGE ENTRE-T-IL ICI ? Uniquement par les COURANTS : la carte
%  est magnetostatique et n'utilise que le FONDAMENTAL (inst_currents lit
%  r.I1c et r.I2c). M.opt.skew_harm agit sur les branches de cage
%  HARMONIQUES, qui n'alimentent pas cette carte. On s'attend donc a un
%  effet de l'ordre du demi-pour-cent -- et c'est cela qu'il faut verifier,
%  car B5 annonce un facteur plusieurs centaines sur les branches de
%  denture du SCHEMA, ce qui ne se transporte pas necessairement ici.
clear; clc; t0=tic;
diary('C5b_ripple_chain_out.txt'); diary on;
Npos=19; nT=17; nO=4; Nh=8192;
%  ---- ENTREE DE CHARGE MESUREE, comme RUN_ARTICLE:199 et 350 ----
%  TL n'est PAS une constante : c'est la SERIE TEMPORELLE du couple mesure
%  en EF, qui porte sa propre ondulation. La machine la suit, donc
%  l'ondulation observee sur oc.Tem est la COMPOSITION de la carte et de la
%  charge -- pas la carte seule. C'est la derniere hypothese restante pour
%  expliquer l'ecart aux 101,4 N.m publies.
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
Dch=fullfile(ROOT,'transitoire','en charge');
rdt=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
trq=rdt(fullfile(Dch,'Plot 1.tab'));
fc=trq(:,1)>=1.90;
TppF=max(trq(fc,3))-min(trq(fc,3));   % reference CALCULEE, non transcrite
%  GLISSEMENT. RUN_ARTICLE ne prend pas 0,0188 mais le DEDUIT de la vitesse
%  mesuree du transitoire en charge : s1 = (n_s - n1)/n_s avec
%  n1 = 1469,772776 tr/min. C'est ce glissement-la qui alimente la carte
%  publiee ; imposer 0,0188 decalait deja les courants.
Mtmp=mec.machine_18_5kW(); n1=1469.772776;
s_ch=(Mtmp.ns*60-n1)/(Mtmp.ns*60);
fprintf('=== C5b : ondulation sur la chaine du MAILLAGE ===\n');
fprintf('  %d positions sur 2 pas d''encoche | reference EF %.1f N.m c-c\n',Npos,TppF);
fprintf('  glissement DEDUIT de n1 = %.6f tr/min : s = %.6f\n',n1,s_ch);
fprintf('  operateur : base P1, pavage nT=%d nO=%d, N_h = %d\n\n',nT,nO,Nh);

nm={'vrillage ON (etat publie)','vrillage OFF (voie a)'};
Tpp=nan(1,2); Tmoy=nan(1,2); Tsch=nan(1,2);
for k=1:2
    M=mec.machine_18_5kW();
    if k==2, M.opt.skew_harm=0; end
    ctx=mec.build_context(M); G=ctx.G; W=ctx.W; BH=ctx.BH;
    %  Configuration P1 convergee (B2) au lieu du defaut 6/2 P0.
    A1=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
    ctx.AG=A1; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
    r=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
    Tsch(k)=r.Tem;
    Pth=2*(2*pi/M.Ns); th0=linspace(0,Pth,Npos); Tk=nan(1,Npos);
    for q=1:Npos
        [i3k,ibk]=loc_currents(M,W,r);
        me=mec.mesh_refined(M,G,6,i3k,th0(q),3,ibk);
        Se=mec.solve_mesh(me,BH,M.opt);
        Us=Se.U(me.gapF.ids(1:me.Ms));
        Ur=Se.U(me.gapF.ids(me.Ms+1:end));
        Tk(q)=me.gapF.torque(Us,Ur);
    end
    Tpp(k)=max(Tk)-min(Tk); Tmoy(k)=mean(Tk);
    fprintf('  %-26s carte BRUTE : Tmoy %8.3f (schema %7.3f, %+6.2f %%) | Tpp %8.3f\n', ...
        nm{k},Tmoy(k),Tsch(k),100*(Tmoy(k)-Tsch(k))/Tsch(k),Tpp(k));

    %  ---- INJECTION DANS L'EQUATION MECANIQUE ----
    %  C'est ICI que se forme le 101,4 N.m du manuscrit. RUN_ARTICLE:362
    %  imprime ppv(oc.Tem(wc)) -- l'ondulation du TRANSITOIRE SIMULE, non
    %  celle de la carte brute. L'inertie (J = 0,17 kg.m2) filtre la carte :
    %  comparer la carte brute aux 105,5 N.m de la reference, qui est elle
    %  aussi une mesure transitoire, etait une erreur de grandeur.
    DT=@(th) interp1(th0,Tk-mean(Tk),mod(th,Pth),'linear');
    oc=mec.dq_startup(ctx,struct('tend',2.0,'TL',[trq(:,1),trq(:,2)],'Trip',DT));
    wc=oc.t>=1.90;                          % MEME fenetre que RUN_ARTICLE:354
    Tpp_dq(k)=max(oc.Tem(wc))-min(oc.Tem(wc)); %#ok<SAGROW>
    Tm_dq(k)=mean(oc.Tem(wc));                 %#ok<SAGROW>
    fprintf('  %-26s APRES dq_startup : Tmoy %8.3f | Tpp %8.3f N.m\n', ...
        '',Tm_dq(k),Tpp_dq(k));
end

fprintf('\n  ---- ondulation crete-a-crete : LA GRANDEUR COMPARABLE ----\n');
fprintf('  (celle du manuscrit, 101.4 N.m, est POST-dq_startup)\n');
fprintf('  %-26s %12s %12s %11s %11s\n', ...
    'etat','carte brute','apres dq','EF (N.m)','ecart/EF');
for k=1:2
    fprintf('  %-26s %12.3f %12.3f %11.1f %10.2f %%\n', ...
        nm{k},Tpp(k),Tpp_dq(k),TppF,100*(Tpp_dq(k)-TppF)/TppF);
end
fprintf('\n  effet de la neutralisation :\n');
fprintf('    sur la carte brute  : %+.3f %%\n',100*(Tpp(2)-Tpp(1))/Tpp(1));
fprintf('    APRES dq_startup    : %+.3f %%\n',100*(Tpp_dq(2)-Tpp_dq(1))/Tpp_dq(1));
fprintf('  valeur publiee 101.4 N.m -- ecart de ce run (post-dq) %+.2f %%\n', ...
    100*(Tpp_dq(1)-101.4)/101.4);
fprintf('  facteur d''attenuation par l''inertie : %.3f\n',Tpp_dq(1)/Tpp(1));
save('C5b_ripple.mat','Tpp','Tpp_dq','Tmoy','Tm_dq','Tsch','Npos');
fprintf('\n  duree %.0f s\n=== C5b termine ===\n',toc(t0));
diary off;

% ======================================================================
function [i3,ib]=loc_currents(M,W,r)
%  Copie de la recette de RUN_ARTICLE (inst_currents) : FONDAMENTAL seul.
    psi1=angle(r.I1c); psi2=angle(r.I2c);
    i3=sqrt(2)*r.I1*[cos(psi1);cos(psi1-2*pi/3);cos(psi1+2*pi/3)];
    Ibar=r.I2*(2*M.m*W.kw1*W.Nph)/M.Nr;
    ib=sqrt(2)*Ibar*cos(M.p*(2*pi*(0:M.Nr-1)'/M.Nr)-psi2);
end
