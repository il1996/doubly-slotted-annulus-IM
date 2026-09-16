function H = harmonics(M, W, numax)
%HARMONICS  Spectre d'harmoniques d'espace du bobinage statorique.
%
%   H = mec.harmonics(M, W, numax) renvoie les harmoniques d'espace ν du
%   bobinage (ν = 6k±1, incluant les harmoniques de DENTURE ν = (Ns/p)k ± 1),
%   avec pour chacun :
%       H.nu    : ordre
%       H.kw    : facteur de bobinage de l'harmonique
%       H.dir   : sens de rotation (+1 direct, -1 inverse)
%       H.sig   : (kw_nu/(nu*kw1))^2  — poids de fuite différentielle
%
%   Convention triphasée : ν = 6k+1 tournent dans le sens DIRECT,
%   ν = 6k−1 dans le sens INVERSE. Les harmoniques de denture
%   (ν = 24k±1 ici : 23, 25, 47, 49) ont |kw_nu| = kw1 — ce sont les plus
%   nocifs (couples parasites asynchrones).
%
%   Sert au modèle multi-harmonique des couples parasites (mec.harmonic_torque).

if nargin<3 || isempty(numax), numax=49; end
p=M.p; Ns=M.Ns; q=W.q; kw1=W.kw1;
beta=M.yq/(Ns/(2*p));
alse=2*pi*p/Ns;                       % angle électrique entre encoches

nus=[]; kws=[]; dirs=[];
for k=1:ceil(numax/6)
    for nu=[6*k-1, 6*k+1]
        if nu>numax, continue; end
        kdn=sin(nu*q*alse/2)/(q*sin(nu*alse/2));
        kpn=sin(nu*beta*pi/2);
        nus(end+1)=nu;            %#ok<AGROW>
        kws(end+1)=kdn*kpn;       %#ok<AGROW>
        dirs(end+1)=mod(nu,6)==1; %#ok<AGROW>  1 si 6k+1 (direct)
    end
end
dirs=double(dirs); dirs(dirs==0)=-1;   % 6k-1 -> inverse

H.nu=nus(:); H.kw=kws(:); H.dir=dirs(:);
H.sig=(H.kw./(H.nu*kw1)).^2;
H.kw1=kw1;
end
