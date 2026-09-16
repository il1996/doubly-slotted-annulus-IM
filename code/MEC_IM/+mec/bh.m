function BH = bh(mat)
%BH  Construit les interpolants magnétiques à partir d'une table B(H).
%
%   BH = mec.bh(mat) prend la structure matériau (mat.BH = [H B]) et
%   renvoie une structure de fonctions vectorisées :
%
%       BH.Hof(B)    : champ H [A/m] pour une induction B [T]
%       BH.Bof(H)    : induction B [T] pour un champ H [A/m]
%       BH.dHdB(B)   : dérivée dH/dB [A/m/T]  (réluctance différentielle)
%       BH.nu(B)     : réluctivité (chord) nu = H/B [A/m/T]  (mu = B/H)
%       BH.Bsat, BH.Hsat : dernier point tabulé
%
%   Conventions et robustesse :
%     * B et H strictement croissants (courbe de 1re aimantation) ;
%     * au-delà du dernier point, extrapolation AIR : dB/dH = mu0
%       (H = Hsat + (B-Bsat)/mu0), ce qui garantit la monotonie et un
%       comportement physique en saturation profonde ;
%     * près de B=0, mu tend vers la perméabilité initiale (pente de la
%       table), évitant toute division par zéro dans nu = H/B.
%
%   Ce module est le SEUL point de dépendance au matériau : changer de
%   nuance de tôle = fournir une autre table à mec.mat_*. La disponibilité
%   de dH/dB (réluctance différentielle) permet un solveur de Newton exact
%   (cf. mec.solve_network), au lieu du point fixe sous-relaxé.

mu0 = 4*pi*1e-7;

H = mat.BH(:,1);
B = mat.BH(:,2);

% Nettoyage : suppression d'éventuels doublons pour une interpolation propre
[B, iu] = unique(B, 'stable');
H = H(iu);

BH.Hsat = H(end);
BH.Bsat = B(end);
BH.mu0  = mu0;

% Perméabilité initiale (pente à l'origine) pour nu près de zéro
mu_init = B(2)/H(2);

% --- H(B) : interpolation linéaire par morceaux + extrapolation air ---
    function Hout = Hof(Bin)
        Bin = Bin(:);
        Hout = zeros(size(Bin));
        in  = abs(Bin) <= B(end);
        out = ~in;
        s   = sign(Bin); s(s==0) = 1;
        % branche interne (table)
        Hout(in) = s(in) .* interp1(B, H, abs(Bin(in)), 'linear');
        % branche externe (air au-delà de la saturation)
        Hout(out) = s(out) .* ( H(end) + (abs(Bin(out)) - B(end))/mu0 );
    end

% --- B(H) : réciproque, extrapolation air ---
    function Bout = Bof(Hin)
        Hin = Hin(:);
        Bout = zeros(size(Hin));
        in  = abs(Hin) <= H(end);
        out = ~in;
        s   = sign(Hin); s(s==0) = 1;
        Bout(in)  = s(in)  .* interp1(H, B, abs(Hin(in)), 'linear');
        Bout(out) = s(out) .* ( B(end) + (abs(Hin(out)) - H(end))*mu0 );
    end

% --- dH/dB : réluctance différentielle (pente locale de H(B)) ---
    function d = dHdB(Bin)
        Bin = abs(Bin(:));
        d = zeros(size(Bin));
        in  = Bin <= B(end);
        out = ~in;
        % pente locale par différences finies sur la table
        dHt = diff(H)./diff(B);          % pente sur chaque segment
        % indice de segment pour chaque point interne
        idx = discretize(Bin(in), B);
        idx(isnan(idx)) = 1;             % B=0 -> premier segment
        d(in)  = dHt(min(idx, numel(dHt)));
        d(out) = 1/mu0;                  % au-delà : pente air
    end

% --- nu(B) = H/B (réluctivité corde) ---
    function n = nu(Bin)
        Bin = Bin(:);
        aB  = abs(Bin);
        n   = zeros(size(Bin));
        small = aB < 1e-6;
        n(small)  = 1/mu_init;           % perméabilité initiale
        n(~small) = Hof(aB(~small)) ./ aB(~small);
    end

% --- densité de co-énergie wco(B) = B*H - w_energie(B) ---
% w_energie(B) = integrale_0^B H dB' (densité d'énergie magnétique)
Btab = B; Htab = H;                        % la table inclut déjà (0,0)
Wtab = cumtrapz(Btab, Htab);               % energie cumulée sur la table
    function wc = wco(Bin)
        Bin = abs(Bin(:));
        we = zeros(size(Bin));
        in = Bin <= B(end);
        we(in)  = interp1(Btab, Wtab, Bin(in), 'linear');
        % au-delà : energie table + partie air (½/mu0 (B^2-Bsat^2))
        we(~in) = Wtab(end) + H(end)*(Bin(~in)-B(end)) ...
                  + 0.5*(Bin(~in)-B(end)).^2/mu0;
        wc = Bin.*Hof(Bin) - we;           % co-énergie = B*H - energie
    end

BH.Hof   = @Hof;
BH.Bof   = @Bof;
BH.dHdB  = @dHdB;
BH.nu    = @nu;
BH.wco   = @wco;
BH.mu_init = mu_init;

end
