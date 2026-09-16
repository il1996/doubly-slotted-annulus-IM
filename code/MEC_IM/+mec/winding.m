function W = winding(M)
%WINDING  Description du bobinage statorique et opérateurs de FMM.
%
%   W = mec.winding(M) construit :
%       W.C     : matrice conducteurs signés par encoche (Ns x 3),
%                 C(k,ph) = somme algébrique des spires de la phase ph
%                 logées dans l'encoche k (aller +, retour -).
%       W.kw1   : facteur de bobinage fondamental (distribution x raccourci)
%       W.kd1, W.kp1 : facteurs de distribution et de raccourcissement
%       W.Nph   : nombre de spires en série par phase
%       W.q     : encoches par pôle et par phase
%       W.slotMMF(i3) : renvoie le vecteur (Ns x 1) de FMM cumulée par
%                       dent (ampères-tours), moyenne retirée, pour un
%                       jeu de courants triphasés i3 = [iA;iB;iC].
%
%   La matrice C reproduit exactement le bobinage double-couche 5/6 du
%   modèle EF de référence. La FMM par dent (escalier de FMM d'entrefer)
%   sert de source dans les branches de dents statoriques du réseau
%   (formulation nodale, cf. mec.build_network).
%
%   Voir aussi : mec.build_network, mec.magnetizing.

Ns = M.Ns;  m = M.m;  p = M.p;

% ------------------------------------------------------------------
%  1. Matrice conducteurs / encoche (machine complète)
% ------------------------------------------------------------------
C = zeros(Ns, m);
w = @(x) mod(x-1,Ns)+1;

% Phase A : couche du bas et couche du haut
for k = 1:numel(M.wind.pDown)
    C(w(M.wind.pDown(k)),1) = C(w(M.wind.pDown(k)),1) + M.wind.sDown(k)*M.Ntc;
end
for k = 1:numel(M.wind.pUp)
    C(w(M.wind.pUp(k)),1) = C(w(M.wind.pUp(k)),1) + M.wind.sUp(k)*M.Ntc;
end

% Phases B et C : décalage d'encoches ; C inversée (convention modèle EF)
shiftB = M.wind.shift(2);
shiftC = M.wind.shift(3);
sgnC   = 1; if M.wind.Cinv, sgnC = -1; end
C(:,2) = circshift(C(:,1),  shiftB);
C(:,3) = sgnC * circshift(C(:,1), shiftC);

W.C = C;

% ------------------------------------------------------------------
%  2. Facteurs de bobinage fondamentaux
% ------------------------------------------------------------------
q    = Ns/(2*p*m);                 % encoches par pôle et par phase
alse = 2*pi*p/Ns;                  % angle électrique entre encoches [rad]
kd1  = sin(q*alse/2)/(q*sin(alse/2));
beta = M.yq/(Ns/(2*p));            % pas relatif (y/pas_polaire)
kp1  = sin(beta*pi/2);
W.q   = q;
W.kd1 = kd1;
W.kp1 = kp1;
W.kw1 = kd1*kp1;

% ------------------------------------------------------------------
%  3. Spires en série par phase
% ------------------------------------------------------------------
% Conducteurs par encoche (2 couches) = layers*Ntc ; en série par phase :
zQ  = M.layers*M.Ntc;              % conducteurs par encoche
W.Nph = zQ*Ns/(2*M.a*m);           % spires en série par phase

% ------------------------------------------------------------------
%  4. Opérateur FMM par dent (escalier d'entrefer)
% ------------------------------------------------------------------
    function F = slotMMF(i3)
        % courant par CONDUCTEUR = courant de phase / voies parallèles
        icond = i3(:) / M.a;
        % ampères-tours par encoche (C = nb de conducteurs signés)
        A = C * icond;             % (Ns x 1)
        % FMM cumulée le long de la périphérie (escalier)
        Fc = cumsum(A);
        % retrait de la moyenne (potentiel homopolaire indéterminé)
        F = Fc - mean(Fc);
    end
W.slotMMF = @slotMMF;

end
