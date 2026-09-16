function id = fem_airgap_ident(M, G, opt)
%FEM_AIRGAP_IDENT  Identification FEM de la perméance d'entrefer dent-à-dent.
%
%   id = mec.fem_airgap_ident(M,G,opt) identifie, par éléments finis
%   magnétostatiques 2D (PDE Toolbox), la perméance d'entrefer Lambda(sigma)
%   entre une dent statorique et une dent rotorique en regard, en fonction
%   de leur décalage tangentiel sigma — démarche de Gyselinck [29]/Silva
%   (banc à deux dents). On en extrait la LARGEUR PHYSIQUE sigma0 du noyau
%   de couplage C1 du MEC (paramètre gouvernant l'ondulation de couple).
%
%   Méthode : analogie électrostatique du potentiel scalaire magnétique
%   (div(mu grad Omega)=0  <->  div(eps grad V)=0). Banc unrolled : dent
%   stator (V=1 via culasse) / entrefer g / dent rotor (V=0 via culasse),
%   parois latérales à flux nul. Perméance (par unité de longueur axiale) :
%       Lambda(sigma) = mu0 * integrale( eps_r |grad V|^2 ) dA     (V=1)
%   calculée pour une série de sigma, puis ajustement
%       Lambda(sigma) = Lfloor + Lmax * exp(-(sigma/sigma0)^2).
%
%   Sortie id : .sigma,.Lambda (FEM), .sigma0,.Lmax,.Lfloor (fit),
%   .tau, .ratio=sigma0/tau, .nelem.
%
%   Voir aussi : mec.mesh_refined, RUN_FEM_IDENT.

if nargin<3, opt=struct(); end
mu0=4*pi*1e-7;
g   = M.g;
w_s = G.taus - M.bs0;
w_r = G.taur - M.br0;
ht  = getdef(opt,'ht', 6e-3);
ty  = getdef(opt,'ty', 4e-3);
tau = 0.5*(G.taus+G.taur);
W   = getdef(opt,'W', 4*tau);
epsFe = getdef(opt,'epsFe',5000);
Hmax  = getdef(opt,'Hmax', g);

sig = linspace(0, 0.9*tau, getdef(opt,'Nsig',13));
Lam = zeros(size(sig));  nel = 0;
for q=1:numel(sig)
    [Lam(q), nel] = solve_bench(sig(q));
end

% Retrait de la fuite lointaine (plancher) : on soustrait le min (couplage
% dent-à-dent = partie variable), puis ajustement GAUSSIEN PUR sur la partie
% de couplage normalisée.
Lfloor=min(Lam);
Lc=Lam-Lfloor;                                   % couplage (au-dessus du plancher)
p0=[max(Lc), 0.6*tau];
ff=@(p) (p(1)*exp(-(sig(:)/p(2)).^2)-Lc(:));     % gaussienne pure, résidu colonne
popt=lsqfit(ff,p0);
% largeur à 1/e (mesure robuste indépendante du modèle)
s1e=interp1(Lc/Lc(1), sig, 1/exp(1), 'linear');
id.sigma=sig; id.Lambda=Lam; id.nelem=nel;
id.Lmax=popt(1); id.sigma0=abs(popt(2)); id.Lfloor=Lfloor;
id.sigma0_1e=s1e; id.tau=tau; id.ratio=id.sigma0/tau; id.ratio_1e=s1e/tau;

% ------------------------------------------------------------------
    function [L,nelem] = solve_bench(sigma)
        model=createpde('electromagnetic','electrostatic');
        model.VacuumPermittivity=8.854e-12;
        rect=@(xc,w,y0,y1)[3;4;xc-w/2;xc+w/2;xc+w/2;xc-w/2;y0;y0;y1;y1];
        gd=[rect(0,W,-ht-ty,g+ht+ty), rect(0,W,g+ht,g+ht+ty), ...
            rect(0,w_s,g,g+ht), rect(sigma,w_r,-ht,0), rect(0,W,-ht-ty,-ht)];
        ns=char('BOX','SY','ST','RT','RY')';
        dl=decsg(gd,'BOX+SY+ST+RT+RY',ns);
        geometryFromEdges(model,dl);
        generateMesh(model,'Hmax',Hmax,'GeometricOrder','linear');

        % matériaux : par face, testés au centroïde d'un élément de la face
        me=model.Mesh; Nodes=me.Nodes;
        for f=1:model.Geometry.NumFaces
            el=findElements(me,'region','Face',f);
            if isempty(el), continue; end
            nn=me.Elements(:,el(1));
            xc=mean(Nodes(1,nn)); yc=mean(Nodes(2,nn));
            er = infe(xc,yc,sigma)*epsFe + (~infe(xc,yc,sigma))*1;
            electromagneticProperties(model,'RelativePermittivity',er,'Face',f);
        end
        % BC : arêtes dont tous les nœuds sont a y=ytop (V=1) ou y=ybot (V=0)
        ytop=g+ht+ty; ybot=-ht-ty;
        etop=[]; ebot=[];
        for e=1:model.Geometry.NumEdges
            ny=Nodes(2, findNodes(me,'region','Edge',e));
            if all(abs(ny-ytop)<1e-7), etop(end+1)=e; end %#ok<AGROW>
            if all(abs(ny-ybot)<1e-7), ebot(end+1)=e; end %#ok<AGROW>
        end
        electromagneticBC(model,'Voltage',1,'Edge',etop);
        electromagneticBC(model,'Voltage',0,'Edge',ebot);

        R=solve(model);
        Ex=R.ElectricField.Ex(:); Ey=R.ElectricField.Ey(:);
        P=R.Mesh.Nodes; T=R.Mesh.Elements(1:3,:); nelem=size(T,2);
        x=reshape(P(1,T),3,[]); y=reshape(P(2,T),3,[]);
        Ar=0.5*abs((x(2,:)-x(1,:)).*(y(3,:)-y(1,:))-(x(3,:)-x(1,:)).*(y(2,:)-y(1,:)));
        xc=mean(x,1); yc=mean(y,1);
        er=infe(xc,yc,sigma)*epsFe + (~infe(xc,yc,sigma))*1;
        E2=mean(Ex(T).^2+Ey(T).^2,1);
        L=mu0*sum(er.*E2.*Ar);
    end

    function b=infe(x,y,sigma)
        b = (y>=g+ht & y<=g+ht+ty) ...                  % culasse stator
          | (y<=-ht & y>=-ht-ty) ...                    % culasse rotor
          | (y>=g & y<=g+ht & abs(x)<=w_s/2) ...        % dent stator
          | (y<=0 & y>=-ht & abs(x-sigma)<=w_r/2);      % dent rotor
    end
end

% ======================================================================
function v=getdef(s,f,d), if isfield(s,f)&&~isempty(s.(f)),v=s.(f);else,v=d;end, end

function p=lsqfit(ff,p0)
p=p0(:);
for it=1:60
    r=ff(p); J=zeros(numel(r),numel(p));
    for j=1:numel(p)
        dp=p; h=max(1e-7,abs(p(j))*1e-4); dp(j)=dp(j)+h;
        J(:,j)=(ff(dp)-r)/h;
    end
    d=-(J.'*J+1e-9*eye(numel(p)))\(J.'*r);
    p=p+d; if norm(d)<1e-13, break; end
end
end
