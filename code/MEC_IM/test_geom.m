clear; clc; addpath(fileparts(mfilename('fullpath')));
M=mec.machine_18_5kW(); G=mec.geometry(M);
fprintf('--- Cotes reelles utilisees ---\n');
fprintf('entrefer g   = %.3f mm\n', M.g*1e3);
fprintf('Rs           = bs2/3 = %.4f mm\n', M.Rs*1e3);
fprintf('\n--- Hauteurs d''encoche (formules fournies) ---\n');
fprintf('G.hs = hs0+hs1+hs2+Rs                       = %.4f mm\n', G.hs*1e3);
fprintf('G.hr = hr0+hr01+br1/2+hr1+br2/2 (figure ANSYS) = %.4f mm\n', G.hr*1e3);
fprintf('  (hc_r = br1/2 = %.4f mm ; ancienne variante recoupee %.4f mm)\n', ...
        G.hc_r*1e3, G.hc_r_cut*1e3);
fprintf('\n--- Dents ---\n');
fprintf('bts (stator) = %.4f mm   [erreur de parallelisme = %.2e]\n', ...
        G.bts*1e3, G.bts_parallel_err);
fprintf('btr (rotor)  = %.4f mm   [min %.4f (haut) .. max %.4f (fond)]\n', ...
        G.btr*1e3, G.btr_min*1e3, G.btr_max*1e3);
fprintf('\n--- Culasses ---\n');
fprintf('hys = %.4f mm     hyr = %.4f mm\n', G.hys*1e3, G.hyr*1e3);
fprintf('\n--- Sections ---\n');
fprintf('Aslot_s = %.2f mm2    Abar = %.2f mm2  (segment ampute %.3f mm2)\n', ...
        G.Aslot_s*1e6, G.Abar*1e6, G.Aseg_r*1e6);
fprintf('Pas d''encoche : taus = %.3f mm   taur = %.3f mm\n', G.taus*1e3, G.taur*1e3);
