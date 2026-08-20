%% RUN_C1  -  Courants de Foucault DANS le réseau (capacité magnétique, C1)
%
%  Implémente et valide C1 : la capacité magnétique de Perho (2002, ch. 5)
%  embarque les courants de Foucault directement dans le réseau de
%  réluctances, sans circuit électrique additionnel — ce que Silva déclare
%  impossible (« une erreur existera toujours dans les modèles MEC »),
%  alors que sa PROPRE référence [16] le fait. C1 lève donc le défaut P6.
%
%  A) Validation INDÉPENDANTE de l'élément capacitif contre la solution
%     analytique exacte 1D d'une tôle : mu_eff = mu*tanh(g*d/2)/(g*d/2).
%  B) Application machine : rétroaction des Foucault sur le flux + pertes,
%     calculées PAR LE RÉSEAU (et non en post-traitement).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
fprintf('=== C1 : courants de Foucault dans le reseau (capacite magnetique) ===\n\n');

M=mec.machine_18_5kW(); G=mec.geometry(M); W=mec.winding(M);
BH=mec.bh(mec.mat_M800_50A());
mu0=4*pi*1e-7; sig=M.fe.sigma; d=M.fe.d;
fprintf('Tole M800-50A : sigma=%.2e S/m, epaisseur d=%.2f mm\n', sig, d*1e3);

%% A) Validation de l'element capacitif vs solution exacte 1D
fprintf('\n--- A) Element capacitif vs solution exacte 1D de la tole ---\n');
fprintf('Modele C1 : mu_eff = mu - j*w*Cm*l/A = mu*(1 - j*w*mu*sigma*d^2/12)\n');
fprintf('Exact     : mu_eff = mu*tanh(gamma*d/2)/(gamma*d/2), gamma=sqrt(j*w*mu*sigma)\n\n');

fN=M.f; flist=[10 25 50 75 100];            % 0 -> 2*f_N
murs=[1000 3000 5000];
fprintf('%6s |', 'f [Hz]');
for mr=murs, fprintf(' mu_r=%-5d(err%%) |',mr); end
fprintf('\n');
for f=flist
    w=2*pi*f; fprintf('%6.0f |',f);
    for mr=murs
        mu=mr*mu0;
        gam=sqrt(1i*w*mu*sig); x=gam*d/2;
        mu_ex = mu*tanh(x)/x;                       % exact
        % modele C1 : Cm d'une region (l,A quelconques -> mu_eff independant)
        l=1; A=1; Cm=mec.iron_capacitance(M,l,A,mu);
        mu_c1 = mu - 1i*w*Cm*l/A;                   % = mu*(1 - j w mu sig d^2/12)
        err=100*abs(mu_c1-mu_ex)/abs(mu_ex);
        fprintf(' %14.3f |',err);
    end
    fprintf('\n');
end
fprintf('\n(erreur = ecart du modele C1 au 1er ordre vs la solution exacte)\n');

% profondeur de peau dans la tole (domaine de validite)
for mr=murs
    dl=sqrt(2/(2*pi*fN*mr*mu0*sig));
    fprintf('  mu_r=%4d : delta_tole(50Hz)=%.2f mm  -> d/delta=%.2f\n',mr,dl*1e3,d/dl);
end

%% B) Application machine : retroaction des Foucault sur le flux
fprintf('\n--- B) Machine : retroaction des Foucault (reseau complexe) ---\n');
nc=4; nr=3; Im=8.32; i3=sqrt(2)*Im*[1;-0.5;-0.5];
mesh=mec.mesh_refined(M,G,nc,i3,0,nr,[]);
S0=mec.solve_mesh(mesh,BH,M.opt);          % solve non lineaire (point de fonctionnement)

% perméabilité locale, réluctance et capacité par branche de fer
Nb=numel(mesh.a); isFe=mesh.iron(:);
B=abs(S0.B(:)); Hh=BH.Hof(B);
mu=zeros(Nb,1);
small=B<1e-6; mu(small)=BH.mu_init; mu(~small)=B(~small)./max(Hh(~small),1e-12);
Rm=zeros(Nb,1); Cm=zeros(Nb,1);
Rm(isFe)=mesh.l(isFe)./(mu(isFe).*mesh.A(isFe));
Cm(isFe)=mec.iron_capacitance(M, mesh.l(isFe), mesh.A(isFe), mu(isFe));
mesh.Rm=Rm; mesh.Cm=Cm;

w=2*pi*M.f;
Sc = mec.solve_mesh_complex(mesh, w, M.opt);          % AVEC Foucault
mesh0=mesh; mesh0.Cm=zeros(Nb,1);
Sn = mec.solve_mesh_complex(mesh0, w, M.opt);         % SANS Foucault (reference)

% flux fondamental d'entrefer dans les deux cas
    function bg1=gapfund(mesh,S,M,G)   % renvoie le fondamental COMPLEXE
        isg=false(numel(mesh.a),1); isg(mesh.gapfirst:end)=true;
        isg=isg & (mesh.a(:)<=mesh.Ms);
        phig=accumarray(mesh.a(isg),S.Phi(isg),[mesh.Ms 1]);
        Bc=phig./(mesh.dth_s*mesh.Rg*M.L);
        bg1=(1/pi)*sum(Bc.*exp(-1i*M.p*mesh.th_s).*mesh.dth_s);
    end
bg_c=gapfund(mesh,Sc,M,G); bg_n=gapfund(mesh0,Sn,M,G);
fprintf('Bg1 SANS Foucault = %.4f T  (phase %+.3f deg)\n',abs(bg_n),angle(bg_n)*180/pi);
fprintf('Bg1 AVEC Foucault = %.4f T  (phase %+.3f deg)\n',abs(bg_c),angle(bg_c)*180/pi);
fprintf('  -> effet sur l''AMPLITUDE : %+.3f %%  (2e ordre : negligeable)\n',...
    100*(abs(bg_c)-abs(bg_n))/abs(bg_n));
fprintf('  -> effet sur la PHASE     : %+.3f deg  (1er ordre : C''EST la perte)\n',...
    (angle(bg_c)-angle(bg_n))*180/pi);
fprintf('Pertes Foucault calculees PAR LE RESEAU = %.1f W\n',Sc.Peddy_tot);

% comparaison au post-traitement classique (Foucault seuls, meme sigma/d)
kc_cl = pi^2*sig*d^2/(6*M.fe.rho);          % coeff. Foucault CLASSIQUE [W/kg/(T.Hz)^2]
mass = M.fe.rho*mesh.A(isFe).*mesh.l(isFe);
Bfe = abs(S0.B(isFe));
Ppp = sum(kc_cl*(M.f*Bfe).^2 .* mass);
fprintf('Pertes Foucault en post-traitement classique  = %.1f W  (ecart %+.1f %%)\n',...
    Ppp, 100*(Sc.Peddy_tot-Ppp)/Ppp);
fprintf('\n(kc classique = pi^2*sigma*d^2/(6*rho) = %.2e ; kc Bertotti ajuste = %.2e\n',...
    kc_cl, M.iron.kc);
fprintf(' -> le kc de Bertotti est ~%.0fx plus grand : il absorbe hysteresis+exces,\n',M.iron.kc/kc_cl);
fprintf('    qui ne produisent PAS de retroaction. C1 ne modelise que les Foucault.)\n');

fprintf('\n=== C1 : Foucault embarques dans le reseau -> defaut P6 leve ===\n');
