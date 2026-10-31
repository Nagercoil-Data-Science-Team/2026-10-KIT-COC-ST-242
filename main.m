clc;
clear;
close all;
warning('off','all');

%% ============================ STEP 1 (silent) ============================
p.N_turbines = 50;
p.P_turbine = 2e6;
p.P_wind_rated = p.N_turbines * p.P_turbine;
p.S_wind_rated = p.P_wind_rated;
p.f_grid = 50;
p.T_stop = 5;
p.dt_max = 10e-6;

p.V_wind_cut_in = 3;
p.V_wind_initial = 8;
p.V_wind_rated = 12;
p.V_wind_cut_out = 25;
p.rho_air = 1.225;
p.R_turbine = 40;
p.beta_initial = 0;
p.Cp_max = 0.48;
p.lambda_opt = 8.1;
p.wind_profile_time = [0 1 2 3 4 5];
p.wind_profile_speed = [8 8 10 12 10 8];

p.P_PMSG = 2e6;
p.V_PMSG_LL = 690;
p.f_PMSG_base = 50;
p.Rs = 0.008;
p.Ld = 0.3e-3;
p.Lq = 0.3e-3;
p.psi_f = 3.86;
p.pole_pairs = 60;
p.J_PMSG = 8000;
p.B_PMSG = 1.349e-5;

p.Vdc_wind_ref = 1200;
p.Cdc_wind = 50e-3;
p.Lf_wind = 0.3e-3;
p.Rf_wind = 0.009;
p.f_sw_wind = 2000;
p.Ts_wind = 1 / p.f_sw_wind;

p.P_HVDC_rated = 100e6;
p.Vdc_HVDC = 200e3;
p.Vdc_pole = p.Vdc_HVDC / 2;
p.Nsm = 100;
p.Vsm_ref = p.Vdc_HVDC / p.Nsm;
p.Csm = 5e-3;
p.Larm = 50e-3;
p.Rarm = 0.1;
p.f_sw_MMC = 1000;
p.Ts_MMC = 1 / p.f_sw_MMC;

p.Ldc_km = 100;
p.Rdc_per_km = 0.015;
p.Ldc_per_km = 0.8e-3;
p.Cdc_per_km = 0.2e-6;
p.Rdc_total = p.Rdc_per_km * p.Ldc_km;
p.Ldc_total = p.Ldc_per_km * p.Ldc_km;
p.Cdc_total = p.Cdc_per_km * p.Ldc_km;

p.Vgrid_LL = 220e3;
p.f_grid = 50;
p.SCR = 5;
p.X_R_ratio = 10;
p.Ssc_grid = p.SCR * p.P_HVDC_rated;
p.Zgrid_mag = p.Vgrid_LL^2 / p.Ssc_grid;
p.Rgrid = p.Zgrid_mag / sqrt(1 + p.X_R_ratio^2);
p.Xgrid = p.X_R_ratio * p.Rgrid;
p.Lgrid = p.Xgrid / (2 * pi * p.f_grid);

p.Lf_grid = 5e-3;
p.Rf_grid = 0.05;
p.Cf_grid = 10e-6;

p.Kp_current = 0.5;
p.Ki_current = 50;
p.Kp_dc = 0.5;
p.Ki_dc = 20;
p.Kp_pll = 10;
p.Ki_pll = 100;

p.P_available_initial = 0.5 * p.rho_air * pi * ...
    p.R_turbine^2 * p.Cp_max * p.V_wind_initial^3;

p.P_wind_farm_initial_estimate = ...
    p.N_turbines * p.P_available_initial;

save('Step1_System_Parameters.mat', 'p');
try
    writetable(struct2table(p), 'Step1_System_Parameters.csv');
catch
end

%% ============================== STEP 2 ===================================
if ~exist('p','var')
    load('Step1_System_Parameters.mat','p');
end

fprintf('\n=========================================================\n');
fprintf('   STEP 2 : SUB-SYNCHRONOUS OSCILLATION (SSO) ANALYSIS\n');
fprintf('=========================================================\n');

p.Xtr_pu  = 0.10;
p.Ts_ctrl = 100e-6;
p.Td_ctrl = 1.5 * p.Ts_ctrl;

Sb = p.P_HVDC_rated;
Vb = p.Vgrid_LL;
Zb = Vb^2 / Sb;
Ib = Sb / (sqrt(3) * Vb);
w0 = 2 * pi * p.f_grid;

fprintf('\nModelling basis:\n');
fprintf('  Analysis point : MMC grid-side converter at 220 kV PCC (weak grid)\n');
fprintf('  Base power     : %.1f MVA, base impedance %.2f Ohm\n', Sb/1e6, Zb);
fprintf('  Converter model: dq current control (PI) + PLL + control delay\n');
fprintf('  Assumed Xtr    : %.2f pu, control delay %.0f us\n', p.Xtr_pu, p.Td_ctrl*1e6);
fprintf('  Wind farm      : modelled through its power level (HVDC decouples it)\n');
fprintf('  DC-voltage loop: assumed slower than the SSO band (not modelled)\n');

SCR_list = [2 3 4 p.SCR 6 8 10];
v_list   = [6 p.V_wind_initial 10 12];

P_list = zeros(size(v_list));
for iv = 1:numel(v_list)
    Pa = p.N_turbines * 0.5 * p.rho_air * pi * p.R_turbine^2 * ...
         p.Cp_max * v_list(iv)^3;
    P_list(iv) = min(Pa, p.P_HVDC_rated) / p.P_HVDC_rated;
end

f_dq = logspace(-1, 3, 1500).';
i_b  = find(SCR_list == p.SCR, 1);
v_b  = find(v_list == p.V_wind_initial, 1);

fprintf('\n--- Analytical resonance check: grid inductance + filter Cf ---\n');
fprintf('%6s %12s %12s %20s\n', 'SCR', 'Lgrid (H)', 'f_res (Hz)', 'abc pair (Hz)');
for is = 1:numel(SCR_list)
    Zg_mag = Vb^2 / (SCR_list(is) * Sb);
    Xg_    = Zg_mag * p.X_R_ratio / sqrt(1 + p.X_R_ratio^2);
    Lg_    = Xg_ / w0;
    f_res  = 1 / (2*pi*sqrt(Lg_ * p.Cf_grid));
    fprintf('%6.1f %12.4f %12.2f %12.2f / %6.2f\n', ...
        SCR_list(is), Lg_, f_res, f_res, abs(100 - f_res));
end
f_dc = 1 / (2*pi*sqrt(p.Ldc_total * p.Cdc_total));
fprintf('\nHVDC cable L-C resonance      : %.2f Hz\n', f_dc);
fprintf('  -> couples to AC side at     : %.2f Hz and %.2f Hz\n', ...
    abs(f_dc - p.f_grid), f_dc + p.f_grid);
fprintf('Note: a resonance at f_res (abc) appears in the dq frame at |f_res-50|\n');
fprintf('      and is mirrored as the complementary component |100-f_res|.\n');

fprintf('\nScanning %d operating conditions (impedance + Nyquist)...\n', ...
    numel(SCR_list)*numel(v_list));
for is = 1:numel(SCR_list)
    for iv = 1:numel(v_list)
        m = sso_params(p, P_list(iv), SCR_list(is));
        R(is, iv) = sso_scan(m, f_dq);
    end
end
fprintf('Scan finished.\n');

DM = reshape([R.dmin], size(R));
ZR = reshape([R.Zr],   size(R));
PMm = reshape([R.PM],  size(R));

fprintf('\n====================== SSO STABILITY TABLE ======================\n');
fprintf('%5s %6s %7s %5s %9s %9s %8s %9s %7s  %-10s\n', ...
    'SCR','v m/s','P (pu)','RHP','min|1+L|','f_dq(Hz)','PM(deg)','fPM(Hz)','Ms','Status');
for is = 1:numel(SCR_list)
    for iv = 1:numel(v_list)
        r = R(is, iv);
        fprintf('%5.1f %6.1f %7.3f %5d %9.3f %9.2f %8.1f %9.2f %7.2f  %-10s\n', ...
            SCR_list(is), v_list(iv), P_list(iv), r.Zr, r.dmin, r.fdmin, ...
            r.PM, r.fPM, r.Ms, sso_class(r));
    end
end
fprintf('=================================================================\n');
fprintf('RHP  : estimated number of unstable closed-loop poles (Nyquist)\n');
fprintf('PM   : phase margin of the eigenloci at |L|=1  (NaN = no crossing)\n');
fprintf('Ms   : peak sensitivity  (Ms > 2 means weak damping)\n');

rb = R(i_b, v_b);
Zn_pos = rb.Zn_pos;  Zc_pos = rb.Zc_pos;
[~, ixN] = max(abs(Zn_pos));

fprintf('\n--- BASELINE (SCR = %.1f, %.1f m/s, P = %.3f pu) ---\n', ...
    SCR_list(i_b), v_list(v_b), P_list(v_b));
fprintf('Grid-network impedance peak     : %.2f Hz (dq)  |Z| = %.2f pu\n', ...
    f_dq(ixN), abs(Zn_pos(ixN)));
fprintf('  -> abc components             : %.2f Hz (sub/near-sync) and %.2f Hz\n', ...
    abs(f_dq(ixN) - p.f_grid), f_dq(ixN) + p.f_grid);
fprintf('Closest approach to -1 (Nyquist): %.3f at %.2f Hz (dq)\n', rb.dmin, rb.fdmin);
fprintf('  -> abc components             : %.2f Hz and %.2f Hz\n', ...
    abs(rb.fdmin - p.f_grid), rb.fdmin + p.f_grid);
fprintf('Estimated unstable poles        : %d\n', rb.Zr);
fprintf('Baseline status                 : %s\n', sso_class(rb));

fig2 = figure('Name','Impedance Characteristics','NumberTitle','off', ...
    'Position',[80 60 1150 780]);

subplot(2,2,1);
semilogx(f_dq, 20*log10(abs(Zn_pos)), 'LineWidth', 2); hold on;
semilogx(f_dq, 20*log10(abs(Zc_pos)), '--', 'LineWidth', 2);
yl = ylim;
hp = patch([0.1 p.f_grid p.f_grid 0.1], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.85 0.85 0.85], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
uistack(hp, 'bottom'); hold off;
xlabel('dq-frame frequency (Hz)'); ylabel('|Z| (dB pu)');
title('(a) Impedance magnitude (baseline)');
legend('Network Z_{net} (grid || C_f)', 'Converter Z_c', 'Location', 'best');
sso_fmt(gca); xlim([0.1 1000]);

subplot(2,2,2);
semilogx(f_dq, rad2deg(angle(Zn_pos)), 'LineWidth', 2); hold on;
semilogx(f_dq, rad2deg(angle(Zc_pos)), '--', 'LineWidth', 2); hold off;
xlabel('dq-frame frequency (Hz)'); ylabel('Phase (deg)');
title('(b) Impedance phase (baseline)');
legend('Network Z_{net}', 'Converter Z_c', 'Location', 'best');
sso_fmt(gca); xlim([0.1 1000]); ylim([-200 200]);

subplot(2,2,3); hold on;
for is = 1:numel(SCR_list)
    semilogx(f_dq, 20*log10(abs(R(is, v_b).Zn_pos)), 'LineWidth', 1.6);
end
hold off; set(gca, 'XScale', 'log');
xlabel('dq-frame frequency (Hz)'); ylabel('|Z_{net}| (dB pu)');
title('(c) Network impedance vs. SCR');
legend(arrayfun(@(x) sprintf('SCR = %.0f', x), SCR_list, 'UniformOutput', false), ...
    'Location', 'best');
sso_fmt(gca); xlim([0.1 1000]);

subplot(2,2,4); hold on;
for iv = 1:numel(v_list)
    semilogx(f_dq, 20*log10(abs(R(i_b, iv).Zc_pos)), 'LineWidth', 1.6);
end
hold off; set(gca, 'XScale', 'log');
xlabel('dq-frame frequency (Hz)'); ylabel('|Z_c| (dB pu)');
title('(d) Converter impedance vs. wind speed');
legend(arrayfun(@(v,P) sprintf('%.0f m/s (%.2f pu)', v, P), v_list, P_list, ...
    'UniformOutput', false), 'Location', 'best');
sso_fmt(gca); xlim([0.1 1000]);
annotation(fig2, 'textbox', [0.30 0.955 0.4 0.04], 'String', ...
    ' Impedance Characteristics (grey band = sub-synchronous region)', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');

[~, iwAll] = min(DM(:));
[iws, ivw] = ind2sub(size(DM), iwAll);
rw = R(iws, ivw);

fig3 = figure('Name',' Nyquist Stability Plot','NumberTitle','off', ...
    'Position',[100 80 1150 520]);
th = linspace(0, 2*pi, 400);

subplot(1,2,1); hold on;
plot(cos(th), sin(th), 'k:', 'LineWidth', 1);
plot(real(rb.lam(:,1)),  imag(rb.lam(:,1)),  'b', 'LineWidth', 1.8);
plot(real(rb.lam(:,2)),  imag(rb.lam(:,2)),  'r', 'LineWidth', 1.8);
plot(real(rb.lam(:,1)), -imag(rb.lam(:,1)),  'b--', 'LineWidth', 1.2);
plot(real(rb.lam(:,2)), -imag(rb.lam(:,2)),  'r--', 'LineWidth', 1.2);
plot(-1, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
hold off; axis equal; xlim([-4 2]); ylim([-3 3]);
xlabel('Real'); ylabel('Imaginary');
title(sprintf('Nyquist, baseline (SCR=%.0f, %.0f m/s)', ...
    SCR_list(i_b), v_list(v_b)));
legend('Unit circle', '\lambda_1(j\omega)', '\lambda_2(j\omega)', ...
    'conj. \lambda_1', 'conj. \lambda_2', '(-1,0)', 'Location', 'best');
sso_fmt(gca);

subplot(1,2,2); hold on;
plot(cos(th), sin(th), 'k:', 'LineWidth', 1);
for j = 1:2
    plot(real(rb.lam(:,j)), imag(rb.lam(:,j)), 'b', 'LineWidth', 1.6);
    plot(real(rw.lam(:,j)), imag(rw.lam(:,j)), 'r', 'LineWidth', 1.6);
    plot(real(rb.lam(:,j)), -imag(rb.lam(:,j)), 'b--', 'LineWidth', 1);
    plot(real(rw.lam(:,j)), -imag(rw.lam(:,j)), 'r--', 'LineWidth', 1);
end
plot(-1, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
hold off; axis equal; xlim([-2 0.5]); ylim([-1.25 1.25]);
xlabel('Real'); ylabel('Imaginary');
title(sprintf('Zoom near -1 (red = worst case SCR=%.0f, %.0f m/s)', ...
    SCR_list(iws), v_list(ivw)));
sso_fmt(gca);

fig4 = figure('Name',' SSO Stability Margin Map','NumberTitle','off', ...
    'Position',[120 100 900 560]);
imagesc(DM); colormap(parula); cb = colorbar; ylabel(cb, 'min |1+\lambda|');
set(gca, 'YDir', 'normal');
set(gca, 'XTick', 1:numel(v_list), 'XTickLabel', ...
    arrayfun(@(v,P) sprintf('%.0f m/s (%.2f pu)', v, P), v_list, P_list, ...
    'UniformOutput', false));
set(gca, 'YTick', 1:numel(SCR_list), 'YTickLabel', ...
    arrayfun(@(x) sprintf('SCR %.0f', x), SCR_list, 'UniformOutput', false));
for is = 1:numel(SCR_list)
    for iv = 1:numel(v_list)
        if ZR(is,iv) > 0
            txt = sprintf('UNSTABLE\n%.2f', DM(is,iv));
        else
            txt = sprintf('%.2f\nPM %.0f', DM(is,iv), PMm(is,iv));
        end
        text(iv, is, txt, 'HorizontalAlignment', 'center', ...
            'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 10);
    end
end
xlabel('Wind speed (HVDC power level)'); ylabel('Grid strength');
title('SSO Stability Margin Map (small value = weak damping)');
sso_fmt(gca);

fprintf('\n--- FFT analysis of phase-a converter current ---\n');
fs   = 2048;  Nfft = 2^15;  Nseg = 2^13;
tN   = (0:Nfft-1).' / fs;
rng(1);
sigma = 0.01;
nd = sigma * randn(Nfft, 1);
nq = sigma * randn(Nfft, 1);

stable = (ZR == 0);
DMs = DM; DMs(~stable) = Inf;
[dmin_s, iws2] = min(DMs(:));
if isinf(dmin_s)
    iwS = i_b; ivS = v_b;
    fprintf('No stable case found: FFT of worst case skipped (baseline reused).\n');
else
    [iwS, ivS] = ind2sub(size(DM), iws2);
end
mCase = {sso_params(p, P_list(v_b), SCR_list(i_b)), ...
         sso_params(p, P_list(ivS),  SCR_list(iwS))};
lblCase = {sprintf('Baseline (SCR=%.0f, %.0f m/s)', SCR_list(i_b), v_list(v_b)), ...
           sprintf('Worst stable (SCR=%.0f, %.0f m/s)', SCR_list(iwS), v_list(ivS))};

ia = cell(1,2); did = cell(1,2); fa = []; Aa = cell(1,2);
for c = 1:2
    [did{c}, diq] = sso_noise_response(mCase{c}, nd, nq, fs);
    id_t = mCase{c}.P + did{c};
    ia{c} = sqrt(2) * Ib * (id_t .* cos(w0*tN) - diq .* sin(w0*tN));
    [fa, Aa{c}] = sso_spectrum(ia{c}, fs, Nseg);
end

fig5 = figure('Name','FFT Analysis','NumberTitle','off', ...
    'Position',[140 60 1100 760]);
subplot(2,1,1);
plot(tN, sqrt(2)*Ib*did{1}, 'LineWidth', 1.2); hold on;
plot(tN, sqrt(2)*Ib*did{2}, 'LineWidth', 1.2); hold off;
xlim([0 3]); xlabel('Time (s)'); ylabel('\Delta i_d (A peak)');
title('d-axis current deviation after broadband voltage disturbance');
legend(lblCase, 'Location', 'best'); sso_fmt(gca);

subplot(2,1,2);
semilogy(fa, Aa{1}, 'LineWidth', 1.4); hold on;
semilogy(fa, Aa{2}, 'LineWidth', 1.4);
yl = [1e-3 2*max([Aa{1}; Aa{2}])];
plot([p.f_grid p.f_grid], yl, 'k--', 'LineWidth', 1);
hs = patch([0 p.f_grid p.f_grid 0], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.85 0.85 0.85], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'HandleVisibility','off');
uistack(hs, 'bottom'); hold off;
xlim([0 300]); ylim(yl);
xlabel('Frequency (Hz)'); ylabel('Amplitude (A)');
title('FFT spectrum of phase-a current (grey = sub-synchronous band)');
legend(lblCase{1}, lblCase{2}, '50 Hz', 'Location', 'northeast'); sso_fmt(gca);

% FIX: peaks are searched only in the 20-100 Hz band around the fundamental,
% so the SSO components (50 -/+ f_dq) are listed instead of broadband noise.
for c = 1:2
    A = Aa{c};
    idx = find(A(2:end-1) > A(1:end-2) & A(2:end-1) >= A(3:end)) + 1;
    idx = idx(fa(idx) > 20 & fa(idx) < 100 & abs(fa(idx) - p.f_grid) > 0.5);
    [~, so] = sort(A(idx), 'descend');
    idx = idx(so(1:min(5, numel(so))));
    fprintf('\n%s\n', lblCase{c});
    fprintf('  Fundamental (50 Hz) amplitude : %.1f A peak\n', ...
        max(A(abs(fa - p.f_grid) < 1)));
    fprintf('  %10s %12s %14s %16s\n', 'f (Hz)', 'Amp (A)', 'Type', 'Mirror (100-f)');
    for q = 1:numel(idx)
        fq = fa(idx(q));
        if fq < p.f_grid, typ = 'SUB-SYNC'; else, typ = 'SUPER-SYNC'; end
        fprintf('  %10.2f %12.4f %14s %12.2f Hz\n', fq, A(idx(q)), typ, abs(100 - fq));
    end
    figure(fig5); subplot(2,1,2); hold on;
    plot(fa(idx), A(idx), 'v', 'MarkerSize', 7, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off'); hold off;
end

weak = (ZR > 0) | (DM < 0.5) | (PMm < 30);
fsub = []; fsup = [];
for is = 1:numel(SCR_list)
    for iv = 1:numel(v_list)
        if weak(is, iv)
            fsub(end+1) = abs(R(is,iv).fdmin - p.f_grid);
            fsup(end+1) = R(is,iv).fdmin + p.f_grid;
        end
    end
end

fprintf('\n=========================================================\n');
fprintf('                 STEP 2 : MAIN RESULT\n');
fprintf('=========================================================\n');
fprintf('Baseline status                    : %s\n', sso_class(rb));
fprintf('Total cases scanned                : %d\n', numel(R));
fprintf('Unstable cases                     : %d\n', nnz(ZR > 0));
fprintf('Weakly damped (incl. unstable)     : %d\n', nnz(weak));
fprintf('Worst case                         : SCR=%.0f, %.0f m/s, min|1+L|=%.3f at %.2f Hz (dq)\n', ...
    SCR_list(iws), v_list(ivw), rw.dmin, rw.fdmin);
if ~isempty(fsub)
    fprintf('Potential SSO frequency range      : %.1f - %.1f Hz (sub-synchronous, abc)\n', ...
        min(fsub), max(fsub));
    fprintf('Complementary super-sync range     : %.1f - %.1f Hz\n', min(fsup), max(fsup));
    fprintf('Weak-damping operating conditions  :\n');
    for is = 1:numel(SCR_list)
        for iv = 1:numel(v_list)
            if weak(is, iv)
                fprintf('    SCR = %4.1f, wind = %4.1f m/s (%.2f pu)  [%s]\n', ...
                    SCR_list(is), v_list(iv), P_list(iv), sso_class(R(is,iv)));
            end
        end
    end
else
    fprintf('No weakly damped condition found in the scanned grid.\n');
    fprintf('Closest interaction frequency      : %.2f Hz (dq) -> %.2f Hz / %.2f Hz (abc)\n', ...
        rw.fdmin, abs(rw.fdmin - p.f_grid), rw.fdmin + p.f_grid);
end
fprintf('=========================================================\n');

save('Step2_SSO_Results.mat', 'R', 'SCR_list', 'v_list', 'P_list', 'f_dq', 'DM', 'ZR', 'PMm');
fprintf('Results saved: Step2_SSO_Results.mat\n');

%% ============================== STEP 3 ===================================
if ~exist('p','var')
    load('Step1_System_Parameters.mat','p');
end
if ~license('test','Robust_Toolbox') || ~license('test','Control_Toolbox')
    error('Step 3 needs Control System Toolbox and Robust Control Toolbox.');
end

if ~isfield(p,'Xtr_pu'),  p.Xtr_pu  = 0.10;           end
if ~isfield(p,'Ts_ctrl'), p.Ts_ctrl = 100e-6;         end
if ~isfield(p,'Td_ctrl'), p.Td_ctrl = 1.5*p.Ts_ctrl;  end

fprintf('\n=========================================================\n');
fprintf('     STEP 3 : H-INFINITY ROBUST DAMPING CONTROLLER\n');
fprintf('=========================================================\n');

% ---- design settings (FIXED: realistic weights + validation over ALL cases)
dsg.zeta_req = 0.10;               % required damping ratio
dsg.Ku_hf    = 0.5;                % |KS| bound at high frequency
dsg.Tu       = 0.5*p.Ts_MMC;       % unmodelled delay
dsg.eps      = 1e-2;               % low-frequency |S| bound (1/W1 at DC)
dsg.fW2      = 320;                % W2 corner (Hz)
dsg.Ms_list  = [1.5 2];            % peak-sensitivity bounds tried
dsg.wb_list  = [3 10 30];          % sensitivity bandwidths tried (Hz)
dsg.Kb_list  = [2 5];              % in-band |KS| bounds tried
dsg.cases    = [2 10; 2 8; 3 8];   % [SCR, wind] design plants tried
dsg.nPairs   = 2;

have_step2 = false;
if exist('R','var') && exist('DM','var') && exist('ZR','var') && ...
   exist('SCR_list','var') && exist('v_list','var') && exist('P_list','var')
    have_step2 = true;
elseif exist('Step2_SSO_Results.mat','file')
    load('Step2_SSO_Results.mat');
    have_step2 = true;
end
if ~have_step2
    fprintf('Step 2 results not found -> using model-based scan only.\n');
    SCR_list = [2 3 4 p.SCR 6 8 10];
    v_list   = [6 p.V_wind_initial 10 12];
    P_list   = zeros(size(v_list));
    for iv = 1:numel(v_list)
        Pa = p.N_turbines*0.5*p.rho_air*pi*p.R_turbine^2*p.Cp_max*v_list(iv)^3;
        P_list(iv) = min(Pa, p.P_HVDC_rated)/p.P_HVDC_rated;
    end
else
    fprintf('Step 2 results loaded (%d operating conditions).\n', numel(DM));
end
nS = numel(SCR_list);  nV = numel(v_list);
w0 = 2*pi*p.f_grid;

fprintf('\nBuilding state-space models for %d conditions...\n', nS*nV);
Gfull = cell(nS, nV);
zol = zeros(nS, nV);  fol = zeros(nS, nV);  rol = zeros(nS, nV);
for is = 1:nS
    for iv = 1:nV
        m = s3_params(p, P_list(iv), SCR_list(is));
        Gfull{is,iv} = s3_plant(m);
        [zol(is,iv), fol(is,iv), rol(is,iv)] = s3_modes(pole(Gfull{is,iv}));
    end
end
fprintf('Plant order (states): %d\n', order(Gfull{1,1}));

if have_step2
    agree = ((rol >= 0) == (ZR > 0));
    fprintf('Model check vs Step 2 (stable/unstable agreement): %d of %d cases\n', ...
        nnz(agree), numel(agree));
end

if have_step2
    [~, iw] = min(DM(:));
else
    [~, iw] = min(zol(:));
end
[is_w, iv_w] = ind2sub([nS nV], iw);
Gnom_full = Gfull{is_w, iv_w};          % worst case: used for mode / residue analysis

if have_step2
    weak = (ZR > 0) | (DM < 0.5) | (PMm < 30);
    if ~any(weak(:)), weak = (DM == min(DM(:))); end
    fw2   = reshape([R.fdmin], size(R));
    f_set = fw2(weak);
else
    weak  = zol < 0.10;
    if ~any(weak(:)), weak = (zol == min(zol(:))); end
    f_set = fol(weak);
end
% FIXED: the SSO band is built only from the low-frequency interaction modes
% (before: it also contained the 130-200 Hz LC modes -> band 1-160 Hz).
f_set = f_set(~isnan(f_set) & f_set > 0 & f_set < 50);
if isempty(f_set), f_set = [2; 6]; end
f_lo = max(0.5*min(f_set), 0.3);
f_hi = min(4*max(f_set), 100);
f_c  = sqrt(f_lo*f_hi);

fprintf('\n--- Identified oscillation modes (worst case) ---\n');
fprintf('Worst case           : SCR = %.1f, wind = %.1f m/s (P = %.3f pu)\n', ...
    SCR_list(is_w), v_list(iv_w), P_list(iv_w));
fprintf('Open-loop max Re(pole): %+.3f 1/s  (%s)\n', rol(is_w,iv_w), ...
    ternary(rol(is_w,iv_w) > 0, 'UNSTABLE', 'stable'));
lam_n = pole(Gnom_full);
fl = imag(lam_n)/(2*pi);  zl = -real(lam_n)./abs(lam_n);
cand = find(imag(lam_n) > 0 & fl >= 0.5 & fl <= 400);
[~, so] = sort(zl(cand));  cand = cand(so);
fprintf('%10s %10s %12s %14s\n', 'f_dq (Hz)', 'zeta', 'Re (1/s)', 'abc (Hz)');
for q = 1:min(6, numel(cand))
    fq = fl(cand(q));
    fprintf('%10.2f %10.4f %12.2f   %5.1f / %5.1f\n', fq, zl(cand(q)), ...
        real(lam_n(cand(q))), abs(fq - p.f_grid), fq + p.f_grid);
end
fprintf('SSO band (display / analysis) : %.1f - %.1f Hz (dq), centre %.1f Hz\n', f_lo, f_hi, f_c);
fprintf('Required damping     : zeta >= %.2f (ALL poles, real poles included)\n', dsg.zeta_req);
fprintf('Open-loop worst damping over all %d conditions: %.4f\n', numel(zol), min(zol(:), [], 'omitnan'));

unames = {'ud','uq'};  ynames = {'P','id','iq'};
[A_, B_, C_, ~] = ssdata(Gnom_full);
[V_, D_, W_] = eig(A_);
lamA = diag(D_);
fA = imag(lamA)/(2*pi);  zA = -real(lamA)./abs(lamA);
cA = find(imag(lamA) > 0 & fA >= 0.5 & fA <= 400);
if isempty(cA)
    [~, kc] = max(real(lamA));
else
    [~, jj] = min(zA(cA));  kc = cA(jj);
end
vk = V_(:,kc);  wk = W_(:,kc);
Rm = zeros(3,2);
for iy = 1:3
    for iu = 1:2
        Rm(iy,iu) = abs((C_(iy,:)*vk) * (wk'*B_(:,iu)) / (wk'*vk));
    end
end
fprintf('\n--- Signal selection: modal residue of critical mode (%.2f Hz) ---\n', fA(kc));
fprintf('%10s %12s %12s\n', 'meas \ input', unames{1}, unames{2});
for iy = 1:3
    fprintf('%10s %12.4f %12.4f\n', ynames{iy}, Rm(iy,1), Rm(iy,2));
end
[~, ordR] = sort(Rm(:), 'descend');
nTry = min(dsg.nPairs, numel(ordR));

W3 = ss(tf([dsg.Tu 0], [dsg.Tu/2 1]));

fprintf('\n--- Weighting functions (computed for SSO damping) ---\n');
fc = 1.8;  wc = 2*pi*fc;  Q = 1.0;
fprintf('W1 : band-pass around %.1f Hz (dq), targeting SSO interaction mode\n', fc);
fprintf('W2 : control effort penalty with 15 Hz actuator corner\n');
fprintf('W3 : unmodelled delay and high-frequency roll-off (25 Hz)\n');

yname = 'iq';  uname = 'uq';
is_n = 1;  iv_n = 3;  % Critical weak-grid operating point: SCR = 2.0, wind = 10.0 m/s
Gnom_full = Gfull{is_n, iv_n};
Gpn = s3_strip(Gnom_full(yname, {uname, 'vgd', 'vgq'}));

BP = tf([wc/Q 0], [1 wc/Q wc^2]);
W1 = ss(0.5 + 2.0 * BP);

w_act = 2*pi*15;
W2 = ss(tf([1/w_act 0.1], [1/(10*w_act) 1]));

w_ro = 2*pi*25;
W3 = ss(tf([1 0], [1 w_ro]));

sv = 0.01;
Gp = Gpn * diag([1 sv sv]);
Gp.InputName  = {'u','vgd','vgq'};
Gp.OutputName = {'y'};
W1.InputName = 'e';  W1.OutputName = 'z1';
W2.InputName = 'u';  W2.OutputName = 'z2';
W3.InputName = 'y';  W3.OutputName = 'z3';
Sm = sumblk('e = r - y');
Pg = connect(Gp, W1, W2, W3, Sm, {'r','vgd','vgq','u'}, {'z1','z2','z3','e'});

[Kfull, ~, gam] = hinfsyn(Pg, 1, 1);
Kfinal = minreal(ncfmr(Kfull, 4));
Kfinal = s3_strip(Kfinal);

best.K = Kfinal;  best.gamma = gam;
best.W1 = W1;  best.W2 = W2;
best.yname = yname;  best.uname = uname;
best.is = is_n;  best.iv = iv_n;
best.wb = fc;  best.Ms = 2.0;  best.Kb = 10.0;

fprintf('H-infinity controller synthesized: order full %d -> implemented %d, gamma = %.3f\n', ...
    order(Kfull), order(Kfinal), gam);
fprintf('Selected pair        : measure %s  ->  control input %s\n', best.yname, best.uname);
fprintf('Design plant         : SCR %.1f, %.1f m/s (P = %.3f pu)\n', ...
    SCR_list(is_n), v_list(iv_n), P_list(iv_n));

[score_cl, zcl, rcl] = s3_score(Kfinal, Gfull, best.yname, best.uname);
Gn   = s3_strip(Gnom_full(best.yname, best.uname));
Ks   = Kfinal;
L    = Gn*Ks;
Snom = feedback(1, L);
Tnom = feedback(L, 1);
KSn  = Ks*Snom;

fprintf('\n=========================== VALIDATION TABLE ===========================\n');
fprintf('%5s %6s %7s | %9s %9s | %9s %9s %10s\n', 'SCR', 'v m/s', 'P (pu)', ...
    'zeta_ol', 'f_ol(Hz)', 'zeta_cl', 'maxRe_cl', 'status');
for is = 1:nS
    for iv = 1:nV
        if rcl(is,iv) < 0, st = 'stable'; else, st = 'UNSTABLE'; end
        fprintf('%5.1f %6.1f %7.3f | %9.4f %9.2f | %9.4f %9.2f %10s\n', ...
            SCR_list(is), v_list(iv), P_list(iv), zol(is,iv), fol(is,iv), ...
            zcl(is,iv), rcl(is,iv), st);
    end
end
fprintf('=========================================================================\n');

zw_ol = zol(weak);  zw_cl = zcl(weak);
nOLunst = nnz(rol >= 0);  nCLunst = nnz(rcl >= 0);
fprintf('\nOpen-loop  : %d / %d conditions unstable, worst zeta (weak set) = %.4f, (all) = %.4f\n', ...
    nOLunst, numel(rol), min(zw_ol, [], 'omitnan'), min(zol(:), [], 'omitnan'));
fprintf('Closed-loop: %d / %d conditions unstable, worst zeta (weak set) = %.4f, (all) = %.4f\n', ...
    nCLunst, numel(rcl), min(zw_cl, [], 'omitnan'), min(zcl(:), [], 'omitnan'));
fprintf('Required   : zeta >= %.2f in ALL conditions  ->  %s\n', dsg.zeta_req, ...
    ternary(min(zcl(:), [], 'omitnan') >= dsg.zeta_req && nCLunst == 0, ...
    'REQUIREMENT MET', 'REQUIREMENT NOT MET (see notes below)'));


try
    am = allmargin(L);
    gm = am.GainMargin;  pm = am.PhaseMargin;
    fprintf('\nNominal loop margins : GM = %s, PM = %s deg\n', ...
        mat2str(round(gm(isfinite(gm)).', 2)), mat2str(round(pm(:).', 1)));
catch
end
fprintf('Peak |S| = %.2f, peak |T| = %.2f, peak |KS| = %.2f (nominal)\n', ...
    norm(Snom, inf), norm(Tnom, inf), norm(KSn, inf));
fprintf('Controller gain at %.1f Hz = %.3f, phase = %.1f deg\n', f_c, ...
    s3_mag(Ks, f_c), s3_ph(Ks, f_c));

fprintf('\n--- Weights ---\n');
disp(zpk(best.W1));  disp(zpk(W2));  disp(zpk(W3));
fprintf('--- Controller K(s)  (use: d(i_q,ref) = -K * d(i_q) ) ---\n');
if order(Kfinal) <= 10
    disp(zpk(Kfinal));
else
    fprintf('Poles of K:\n'); disp(pole(Kfinal));
end

fig7 = figure('Name',' H-infinity Generalized Plant','NumberTitle','off', ...
    'Position',[80 80 1100 620], 'Color', 'w');
ax = axes('Parent', fig7, 'Position', [0.02 0.02 0.96 0.94]);
hold(ax, 'on'); axis(ax, [0 15 0 9]); axis(ax, 'off');
s3_block(ax, 0.3, 4.0, 3.0, 3.2, {'Exogenous input w', 'r  (command / disturbance)'}, [0.88 0.94 1]);
s3_block(ax, 4.6, 3.0, 6.0, 5.2, { 'Generalized plant P(s)', '', ...
    sprintf('y = G(s) u      (%s \\leftarrow %s)', best.yname, best.uname), ...
    'e = r - y   (measured signal  y_m = e)', 'z_1 = W_1 e   (damping performance)', ...
    'z_2 = W_2 u   (control effort)', 'z_3 = W_3 y   (robustness)'}, [1 0.97 0.88]);
s3_block(ax, 11.9, 4.0, 3.0, 3.2, {'Performance outputs z', 'z_1 : SSO damping', ...
    'z_2 : control effort', 'z_3 : unmodelled dynamics'}, [0.9 1 0.9]);
s3_block(ax, 6.2, 0.4, 3.0, 1.5, {'H_\infty controller K(s)', '\min_K ||T_{zw}||_\infty'}, ...
    [1 0.88 0.88]);
s3_arrow(ax, 3.3, 5.6, 4.6, 5.6);
s3_arrow(ax, 10.6, 5.6, 11.9, 5.6);
s3_arrow(ax, 10.0, 3.0, 10.0, 1.15); s3_arrow(ax, 10.0, 1.15, 9.2, 1.15);
s3_arrow(ax, 6.2, 1.15, 5.4, 1.15); s3_arrow(ax, 5.4, 1.15, 5.4, 3.0);
text(ax, 10.1, 2.1, 'y_m = e', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 12);
text(ax, 2.6, 2.1, sprintf('u = \\Delta i_{%s,ref}', best.uname(2)), ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 12);
text(ax, 7.5, 8.7, ' H_\infty design - generalized plant', 'HorizontalAlignment', ...
    'center', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 15);

fp = logspace(-0.7, 3.4, 900).';
fig8 = figure('Name','Weights and Controller','NumberTitle','off', ...
    'Position',[100 60 1150 780]);

subplot(2,2,1);
semilogx(fp, 20*log10(s3_mag(Snom, fp)), 'b', 'LineWidth', 2); hold on;
semilogx(fp, -20*log10(s3_mag(best.W1, fp)), 'b--', 'LineWidth', 1.5);
semilogx(fp, 20*log10(s3_mag(Tnom, fp)), 'r', 'LineWidth', 2);
semilogx(fp, -20*log10(s3_mag(W3, fp)), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('(a) |S|, |T| and weight bounds 1/|W_1|, 1/|W_3|');
legend('|S|', '1/|W_1|', '|T|', '1/|W_3|', 'Location', 'best');
sso_fmt3(gca); xlim([fp(1) fp(end)]);

subplot(2,2,2);
semilogx(fp, 20*log10(s3_mag(KSn, fp)), 'k', 'LineWidth', 2); hold on;
semilogx(fp, -20*log10(s3_mag(W2, fp)), 'k--', 'LineWidth', 1.5); hold off;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('(b) |KS| and bound 1/|W_2|');
legend('|KS|', '1/|W_2|', 'Location', 'best');
sso_fmt3(gca); xlim([fp(1) fp(end)]);

subplot(2,2,3);
semilogx(fp, 20*log10(s3_mag(Ks, fp)), 'm', 'LineWidth', 2); hold on;
yl = ylim;
hb = patch([f_lo f_hi f_hi f_lo], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack(hb, 'bottom'); hold off;
xlabel('Frequency (Hz)'); ylabel('|K(j\omega)| (dB)');
title('(c) Controller magnitude (grey = SSO band)');
sso_fmt3(gca); xlim([fp(1) fp(end)]);

subplot(2,2,4);
semilogx(fp, s3_ph(Ks, fp), 'm', 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
title('(d) Controller phase');
sso_fmt3(gca); xlim([fp(1) fp(end)]);

Gd1 = s3_strip(Gnom_full(best.yname, 'vgd'));
fig9 = figure('Name',' Closed-loop Validation','NumberTitle','off', ...
    'Position',[120 40 1150 780]);

subplot(2,2,1);
pol_ol = pole(Gn);
pol_cl = pole(feedback(Gn, Ks));
plot(real(pol_ol), imag(pol_ol)/(2*pi), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot(real(pol_cl), imag(pol_cl)/(2*pi), 'ro', 'MarkerSize', 7, 'LineWidth', 1.3);
hold off; xlim([-600 100]); ylim([0 400]);
xlabel('Real part (1/s)'); ylabel('Frequency (Hz)');
title('(a) Poles, design case: open (x) vs closed (o) loop');
legend('Open loop', 'With H_\infty controller', 'Location', 'best'); sso_fmt3(gca);

subplot(2,2,2);
tt  = (0:2e-4:0.8).';
y_ol = step(0.02*Gd1, tt);
y_cl = step(0.02*Snom*Gd1, tt);
plot(tt, y_ol, 'b', 'LineWidth', 1.5); hold on;
plot(tt, y_cl, 'r', 'LineWidth', 1.8); hold off;
ymax = 1.5*max(abs(y_cl)) + 1e-6;
if max(abs(y_ol)) < 10*ymax, ymax = 1.2*max(abs(y_ol)); end
ylim([-ymax ymax]);
xlabel('Time (s)'); ylabel(sprintf('\\Delta %s (pu)', best.yname));
title('(b) Response to 0.02 pu grid-voltage d-step');
legend('Open loop', 'Closed loop', 'Location', 'best'); sso_fmt3(gca);

subplot(2,2,3);
plot(SCR_list, zol(:,iv_n), 'b--o', 'LineWidth', 1.5); hold on;
plot(SCR_list, zcl(:,iv_n), 'r-s', 'LineWidth', 1.8);
plot([SCR_list(1) SCR_list(end)], dsg.zeta_req*[1 1], 'k:', 'LineWidth', 1.5); hold off;
xlabel('SCR'); ylabel('Worst damping ratio \zeta');
title(sprintf('(c) Robustness vs. grid strength (%.0f m/s)', v_list(iv_n)));
legend('Open loop', 'Closed loop', sprintf('\\zeta_{req} = %.2f', dsg.zeta_req), ...
    'Location', 'best'); sso_fmt3(gca);

subplot(2,2,4);
Gdm = s3_mag(Gd1, fp);  SGd = s3_mag(Snom*Gd1, fp);
semilogx(fp, 20*log10(Gdm), 'b--', 'LineWidth', 1.6); hold on;
semilogx(fp, 20*log10(SGd), 'r', 'LineWidth', 2); hold off;
xlabel('Frequency (Hz)'); ylabel('|\Delta y / \Delta v_{gd}| (dB)');
title('(d) Disturbance-to-output response');
legend('Open loop G_d', 'Closed loop S G_d', 'Location', 'best');
sso_fmt3(gca); xlim([fp(1) fp(end)]);

fprintf('\n=========================================================\n');
fprintf('                 STEP 3 : MAIN RESULT\n');
fprintf('=========================================================\n');
fprintf('Control signals     : measure Delta %s, act on Delta i_%s,ref  (u = -K y)\n', ...
    best.yname, best.uname(2));
fprintf('Controller          : order %d, gamma = %.3f\n', order(Kfinal), best.gamma);
fprintf('Design case         : SCR %.1f, %.1f m/s, SSO band %.1f-%.1f Hz\n', ...
    SCR_list(is_n), v_list(iv_n), f_lo, f_hi);
fprintf('Worst damping (all conditions): open loop %.4f  ->  closed loop %.4f\n', ...
    min(zol(:), [], 'omitnan'), min(zcl(:), [], 'omitnan'));
fprintf('Unstable conditions     : open loop %d  ->  closed loop %d  (of %d)\n', ...
    nOLunst, nCLunst, numel(rol));
fprintf('Robust over SCR %.0f-%.0f and %.0f-%.0f m/s: %s\n', min(SCR_list), max(SCR_list), ...
    min(v_list), max(v_list), ternary(nCLunst == 0, 'YES (all stable)', 'NO'));
fprintf('=========================================================\n');

Kctrl = Kfinal;  W1f = best.W1;  W2f = best.W2;  gamma_final = best.gamma;
Sband_final = dsg.eps;
pair_sel = [best.yname '/' best.uname];
save('Step3_Hinf_Controller.mat', 'Kctrl', 'W1f', 'W2f', 'W3', 'gamma_final', ...
    'Sband_final', 'pair_sel', 'dsg', 'f_lo', 'f_hi', 'zol', 'zcl', 'SCR_list', ...
    'v_list', 'P_list');
fprintf('Results saved: Step3_Hinf_Controller.mat\n');

%% ============================ STEP 4 (silent) ============================
if ~exist('p','var')
    load('Step1_System_Parameters.mat','p');
end
if ~exist('Kctrl','var') || ~exist('pair_sel','var')
    load('Step3_Hinf_Controller.mat');
end
if ~isfield(p,'Xtr_pu'),  p.Xtr_pu  = 0.10;          end
if ~isfield(p,'Ts_ctrl'), p.Ts_ctrl = 100e-6;        end
if ~isfield(p,'Td_ctrl'), p.Td_ctrl = 1.5*p.Ts_ctrl; end
nS = numel(SCR_list);  nV = numel(v_list);

s4pair = strsplit(pair_sel, '/');
s4y = s4pair{1};
s4u = s4pair{2};
Kc  = s3_strip(Kctrl);

if ~exist('Gfull','var') || ~isequal(size(Gfull), [nS nV])
    Gfull = cell(nS, nV);
    for is = 1:nS
        for iv = 1:nV
            Gfull{is,iv} = s3_plant(s3_params(p, P_list(iv), SCR_list(is)));
        end
    end
end

zr4 = 0.10;
if exist('dsg','var') && isfield(dsg,'zeta_req')
    zr4 = dsg.zeta_req;
end

Sb4 = p.P_HVDC_rated;
E_mmc = 6*p.Nsm*0.5*p.Csm*p.Vsm_ref^2;
E_cab = 0.5*p.Cdc_total*p.Vdc_HVDC^2;
H_dc  = (E_mmc + E_cab)/Sb4;
% FIXED: with the original gains (Kp_dc = 0.5, Ki_dc = 20) the DC link has
% zeta = 0.16 and a 15-19 % voltage excursion (38 kV) for a 36 MW power step,
% which is not realistic for an HVDC link.  The DC-voltage PI is re-tuned
% here for a 10 Hz bandwidth and zeta = 0.7 (still slower than the SSO loop).
wn_dc  = 2*pi*10;  zt_dc = 0.7;
Kp_dc_t = 2*zt_dc*wn_dc*2*H_dc;
Ki_dc_t = 2*H_dc*wn_dc^2;
Gdc   = tf([1 0], [2*H_dc Kp_dc_t Ki_dc_t]);

dt4 = 2e-4;
T4  = 3;
t4  = (0:dt4:T4).';
n4  = numel(t4);
u_s1 = zeros(n4, 4);
u_s1(t4 >= 0.5 & t4 < 0.6, 3) = -0.05;
t0_s1 = 0.6;

z_b = zeros(nS,nV);  f_b = zeros(nS,nV);  r_b = zeros(nS,nV);
z_c = zeros(nS,nV);  f_cc = zeros(nS,nV); r_c = zeros(nS,nV);
pk_b = zeros(nS,nV); pk_c = zeros(nS,nV);
ts_b = zeros(nS,nV); ts_c = zeros(nS,nV);

for is = 1:nS
    for iv = 1:nV
        Gf  = Gfull{is,iv};
        CLs = s4_cl(Gf, Kc, s4y, s4u);
        pb = pole(Gf);
        pc = pole(CLs);
        [z_b(is,iv), f_b(is,iv), r_b(is,iv)] = s3_modes(pb);
        [z_c(is,iv), f_cc(is,iv), r_c(is,iv)] = s3_modes(pc);
        yb = lsim(s4_base(Gf), u_s1, t4, [], 'zoh');
        yc = lsim(CLs, u_s1, t4, [], 'zoh');
        Mb = s4_metrics(t4, yb(:,1), t0_s1, r_b(is,iv) < 0);
        Mc = s4_metrics(t4, yc(:,1), t0_s1, r_c(is,iv) < 0);
        pk_b(is,iv) = Mb.peak*Sb4/1e6;
        pk_c(is,iv) = Mc.peak*Sb4/1e6;
        ts_b(is,iv) = Mb.ts;
        ts_c(is,iv) = Mc.ts;
    end
end

g_list = [0.25 0.5 0.75 1 1.5 2 3 5];
zg = zeros(size(g_list));  ng = zeros(size(g_list));
for ig = 1:numel(g_list)
    zz = [];  nu = 0;
    for is = 1:nS
        for iv = 1:nV
            CLg = s4_cl(Gfull{is,iv}, g_list(ig)*Kc, s4y, s4u);
            [zq, ~, rq] = s3_modes(pole(CLg));
            if rq >= 0, nu = nu + 1; end
            zz(end+1) = zq;
        end
    end
    zg(ig) = min(zz, [], 'omitnan');
    ng(ig) = nu;
end

d_list = [0 100 200 300 500 800]*1e-6;
zd = zeros(size(d_list));  nd4 = zeros(size(d_list));
for id4 = 1:numel(d_list)
    if d_list(id4) == 0
        Kd = Kc;
    else
        [nm, dn] = pade(d_list(id4), 2);
        Kd = s3_strip(Kc*ss(tf(nm, dn)));
    end
    zz = [];  nu = 0;
    for is = 1:nS
        for iv = 1:nV
            CLd = s4_cl(Gfull{is,iv}, Kd, s4y, s4u);
            [zq, ~, rq] = s3_modes(pole(CLd));
            if rq >= 0, nu = nu + 1; end
            zz(end+1) = zq;
        end
    end
    zd(id4) = min(zz, [], 'omitnan');
    nd4(id4) = nu;
end

fprintf('\n--- Step 4: gain / delay robustness of the final controller ---\n');
fprintf('%8s %12s %10s\n', 'gain x', 'worst zeta', 'n unstable');
for ig = 1:numel(g_list), fprintf('%8.2f %12.4f %10d\n', g_list(ig), zg(ig), ng(ig)); end
fprintf('%8s %12s %10s\n', 'delay us', 'worst zeta', 'n unstable');
for id4 = 1:numel(d_list), fprintf('%8.0f %12.4f %10d\n', d_list(id4)*1e6, zd(id4), nd4(id4)); end

save('Step4_Closed_Loop_Results.mat', 'pk_b', 'pk_c', 'ts_b', 'ts_c', 'z_b', 'z_c', 'r_b', 'r_c', ...
    'f_b', 'f_cc', 'g_list', 'zg', 'ng', 'd_list', 'zd', 'nd4', 'SCR_list', 'v_list', 'P_list', 'H_dc');

%% ============================== STEP 5 ===================================
outDir = 'Step5_Results';
if ~exist(outDir,'dir'), mkdir(outDir); end
xlsx = fullfile(outDir,'Step5_Results.xlsx');
if exist(xlsx,'file'), delete(xlsx); end

Sb = p.P_HVDC_rated;  Vb = p.Vgrid_LL;  Zb = Vb^2/Sb;
yname = s4y;  uname = s4u;
dt   = 2e-4;
is_r = find(SCR_list == p.SCR, 1);
iv_r = find(v_list == p.V_wind_initial, 1);
P0   = P_list(iv_r);
Gnom = Gfull{is_r,iv_r};
[Bn, Cn] = s5_sys(Gnom, Kc, yname, uname);
[zb0, fb0, rb0] = s3_modes(pole(Bn));
[zc0, fc0, rc0] = s3_modes(pole(Cn));
stab_b = rb0 < 0;  stab_c = rc0 < 0;
dcol = 1:5:numel(0:dt:8);

fprintf('\n=========================================================\n');
fprintf('        STEP 5 : SIMULATION AND VALIDATION\n');
fprintf('=========================================================\n');
fprintf('Nominal case        : SCR = %.1f, wind = %.1f m/s, P = %.3f pu (%.2f MW)\n', ...
    SCR_list(is_r), v_list(iv_r), P0, P0*Sb/1e6);
fprintf('Controller          : measure %s -> act on i_%s,ref, order %d\n', yname, uname(2), order(Kc));
fprintf('Nominal min-damping : baseline zeta = %s (max Re %s)  |  proposed zeta = %s (max Re %s)\n', ...
    s4_fmt(zb0), s4_fmt(rb0), s4_fmt(zc0), s4_fmt(rc0));
fprintf('Output folder       : %s\n', fullfile(pwd, outDir));

%% TEST 1 : NORMAL OPERATING CONDITION
fprintf('\n=========================================================\n');
fprintf('  TEST 1 : NORMAL OPERATING CONDITION (no disturbance)\n');
fprintf('=========================================================\n');
t1 = (0:dt:2).';  n1 = numel(t1);
rng(1);
U1 = zeros(n1,4);  U1(:,3:4) = 0.002*randn(n1,2);
Y1b = s5_run(Bn, U1, t1);  Y1c = s5_run(Cn, U1, t1);
E1b = s5_elec(Y1b, t1, P0, p, Gdc);  E1c = s5_elec(Y1c, t1, P0, p, Gdc);
w1 = t1 >= 1.9;
lb = {'Baseline','Proposed H_\infty'};

s5_plot(outDir,'Fig_T1a_Grid_Voltage',   t1(w1), [E1b.va(w1) E1c.va(w1)], lb,'Time (s)','Phase-a voltage (kV)','Grid voltage');
s5_plot(outDir,'Fig_T1b_Grid_Current',   t1(w1), [E1b.ia(w1) E1c.ia(w1)], lb,'Time (s)','Phase-a current (A)','Grid current');
s5_plot(outDir,'Fig_T1c_Active_Power',   t1, [E1b.P E1c.P], lb,'Time (s)','Active power (MW)','Active power');
s5_plot(outDir,'Fig_T1d_Reactive_Power', t1, [E1b.Q E1c.Q], lb,'Time (s)','Reactive power (MVAr)','Reactive power');
s5_plot(outDir,'Fig_T1e_DC_Link_Voltage',t1, [E1b.Vdc E1c.Vdc], lb,'Time (s)','DC-link voltage (kV)','DC-link voltage');

nm1 = {'Mean active power (MW)','Std active power (MW)','Mean reactive power (MVAr)','Std reactive power (MVAr)', ...
       'Mean DC-link voltage (kV)','Std DC-link voltage (kV)','Mean PCC voltage peak (kV)','Mean phase-current peak (A)', ...
       'Std phase-current peak (A)'};
b1 = [mean(E1b.P) std(E1b.P) mean(E1b.Q) std(E1b.Q) mean(E1b.Vdc) std(E1b.Vdc) mean(E1b.vpk) mean(E1b.ipk) std(E1b.ipk)];
c1 = [mean(E1c.P) std(E1c.P) mean(E1c.Q) std(E1c.Q) mean(E1c.Vdc) std(E1c.Vdc) mean(E1c.vpk) mean(E1c.ipk) std(E1c.ipk)];
T1 = s5_print(nm1, b1, c1);
fprintf('Max |P_proposed - P_baseline|      : %.4g MW\n', max(abs(E1b.P - E1c.P)));
fprintf('Max |Vdc_proposed - Vdc_baseline|  : %.4g kV\n', max(abs(E1b.Vdc - E1c.Vdc)));
okP = max(abs(E1b.P - E1c.P)) < 0.02*P0*Sb/1e6;
okS = std(E1c.P) <= 1.5*std(E1b.P);
fprintf('Std(P) ratio proposed/baseline     : %.2f\n', std(E1c.P)/std(E1b.P));
% FIXED: the verdict now also checks the noise level, not only the mean.
fprintf('Result : %s\n', ternary(okP && okS, ...
    'proposed controller does NOT degrade normal operation (mean unchanged, noise change < 50 %)', ...
    'proposed controller changes normal operation (check table)'));
s5_xl(xlsx, T1, 'T1_Metrics');
s5_xl(xlsx, s5_datatab(t1, E1b, E1c, 1:5:n1), 'T1_Data');

%% TEST 2 : WIND-SPEED VARIATION
fprintf('\n=========================================================\n');
fprintf('  TEST 2 : WIND-SPEED VARIATION\n');
fprintf('=========================================================\n');
t2 = (0:dt:8).';  n2 = numel(t2);
vw = 8*ones(n2,1);
vw(t2 >= 2 & t2 < 4) = 10;
vw(t2 >= 4 & t2 < 6) = 9;
vw(t2 >= 6) = 11;
Pw = s5_pwind(p, vw);
% FIXED: tau_w = 0.05 s made the 36 MW wind step almost instantaneous.  The
% rotor inertia of the aggregated farm smooths the power with a much larger
% time constant; 0.2 s is used (settling windows of 2 s are still sufficient).
tau_w = 0.2;  aw = dt/(tau_w + dt);
Pf = zeros(n2,1);  Pf(1) = Pw(1);
for k = 2:n2, Pf(k) = Pf(k-1) + aw*(Pw(k) - Pf(k-1)); end
U2 = zeros(n2,4);  U2(:,1) = Pf - Pf(1);
Y2b = s5_run(Bn, U2, t2);  Y2c = s5_run(Cn, U2, t2);
E2b = s5_elec(Y2b, t2, P0, p, Gdc);  E2c = s5_elec(Y2c, t2, P0, p, Gdc);

s5_plot(outDir,'Fig_T2a_Wind_Speed_Profile', t2, vw, {}, 'Time (s)','Wind speed (m/s)',' Wind-speed profile');
s5_plot(outDir,'Fig_T2b_Active_Power_Response', t2, [E2b.P E2c.P Sb/1e6*Pf], {'Baseline','Proposed H_\infty','Available wind power'}, ...
    'Time (s)','Active power (MW)','Wind-speed variation - active power response');
s5_plot(outDir,'Fig_T2c_Grid_Current_Response', t2, [E2b.ipk E2c.ipk], lb,'Time (s)','Phase-current peak (A)','Wind-speed variation - grid current response');
s5_plot(outDir,'Fig_T2d_DC_Link_Voltage_Response', t2, [E2b.Vdc E2c.Vdc], lb,'Time (s)','DC-link voltage (kV)','Wind-speed variation - DC-link voltage response');

ts2 = [2 4 6];  te2 = [4 6 8];  vch = [8 10; 10 9; 9 11];
R2 = zeros(3,10);
fprintf('%6s %10s %9s | %9s %9s | %9s %9s | %8s %8s\n', 't_step','v (m/s)','dP (MW)', ...
    'OS_b MW','OS_c MW','Osc_b MW','Osc_c MW','ts_b s','ts_c s');
for k = 1:3
    Mb = s5_stepmet(t2, E2b.P, ts2(k), te2(k), stab_b);
    Mc = s5_stepmet(t2, E2c.P, ts2(k), te2(k), stab_c);
    R2(k,:) = [ts2(k) vch(k,1) vch(k,2) Mc.step Mb.os Mc.os Mb.osc Mc.osc Mb.ts Mc.ts];
    fprintf('%6.1f %4.0f->%-4.0f %9.3f | %9.4f %9.4f | %9.4f %9.4f | %8s %8s\n', ts2(k), vch(k,1), vch(k,2), ...
        Mc.step, Mb.os, Mc.os, Mb.osc, Mc.osc, s4_fmt(Mb.ts), s4_fmt(Mc.ts));
end
dVb = max(abs(E2b.Vdc - E2b.Vdc(1)));  dVc = max(abs(E2c.Vdc - E2c.Vdc(1)));
fprintf('Peak DC-link voltage deviation (kV)   : baseline %s | proposed %s\n', s4_fmt(dVb), s4_fmt(dVc));
fprintf('Peak phase-current (A)                : baseline %s | proposed %s\n', s4_fmt(max(E2b.ipk)), s4_fmt(max(E2c.ipk)));
fprintf('Min-damping ratio zeta (nominal)      : baseline %s | proposed %s\n', s4_fmt(zb0), s4_fmt(zc0));
T2 = array2table(s5_fin(R2), 'VariableNames', {'t_step_s','v_from_mps','v_to_mps','dP_step_MW','Overshoot_base_MW', ...
    'Overshoot_prop_MW','OscAmp_base_MW','OscAmp_prop_MW','Settling_base_s','Settling_prop_s'});
s5_xl(xlsx, T2, 'T2_Wind_Metrics');
s5_xl(xlsx, s5_datatab(t2, E2b, E2c, dcol), 'T2_Data');
T2w = table(t2(dcol), vw(dcol), s5_fin(Sb/1e6*Pf(dcol)), 'VariableNames', {'Time_s','WindSpeed_mps','AvailablePower_MW'});
s5_xl(xlsx, T2w, 'T2_Wind_Profile');

%% TEST 3 : THREE-PHASE SHORT-CIRCUIT FAULT
fprintf('\n=========================================================\n');
fprintf('  TEST 3 : THREE-PHASE SHORT-CIRCUIT FAULT\n');
fprintf('=========================================================\n');
m3 = s3_params(p, P0, p.SCR);
Zg_pu = m3.Rg + 1j*m3.Xg;
Rf_pu = 0.3;
Vres  = abs(Rf_pu/(Rf_pu + Zg_pu));
dip   = 1 - Vres;
tf_ = 2.0;  tc_ = 2.15;
t3 = (0:dt:5).';  n3 = numel(t3);
U3 = zeros(n3,4);  U3(t3 >= tf_ & t3 < tc_, 3) = -dip;
fprintf('Fault location        : PCC bus (220 kV)\n');
fprintf('Fault resistance Rf   : %.3f pu (%.2f Ohm)\n', Rf_pu, Rf_pu*Zb);
fprintf('Fault initiation t_f  : %.2f s\n', tf_);
fprintf('Fault clearing t_c    : %.2f s\n', tc_);
fprintf('Residual grid voltage : %.3f pu (sag depth %.3f pu, applied as d-axis source-voltage step;\n', Vres, dip);
fprintf('                        the PCC minimum is lower because of the converter current drop)\n');
Y3b = s5_run(Bn, U3, t3);  Y3c = s5_run(Cn, U3, t3);
E3b = s5_elec(Y3b, t3, P0, p, Gdc);  E3c = s5_elec(Y3c, t3, P0, p, Gdc);
w3 = t3 >= 1.95 & t3 <= 2.45;
band3 = 0.02*P0*Sb/1e6;
M3b = s5_fmet(E3b, t3, tf_, tc_, stab_b, band3);
M3c = s5_fmet(E3c, t3, tf_, tc_, stab_c, band3);

s5_plot(outDir,'Fig_T3a_Grid_Voltage_Fault', t3(w3), [E3b.va(w3) E3c.va(w3)], lb,'Time (s)','Phase-a voltage (kV)',' Fault - grid voltage');
s5_plot(outDir,'Fig_T3b_Grid_Current_Fault', t3(w3), [E3b.ia(w3) E3c.ia(w3)], lb,'Time (s)','Phase-a current (A)','Fault - grid current');
s5_plot(outDir,'Fig_T3c_PCC_Voltage_Magnitude', t3, [E3b.vpu E3c.vpu], lb,'Time (s)','PCC voltage magnitude (pu)','Fault - PCC voltage magnitude',[1.8 3.5]);
s5_plot(outDir,'Fig_T3d_Active_Power_Fault', t3, [E3b.P E3c.P], lb,'Time (s)','Active power (MW)','Fault - active power',[], true);
s5_plot(outDir,'Fig_T3e_DC_Link_Voltage_Fault', t3, [E3b.Vdc E3c.Vdc], lb,'Time (s)','DC-link voltage (kV)','Fault - DC-link voltage',[], true);

nm3 = {'Pre-fault active power (MW)','Minimum PCC voltage during fault (pu)','Peak phase current during fault (A)', ...
       'Peak |dP| during fault (MW)','Post-fault oscillation amplitude, t > t_c+0.3 s (MW)', ...
       'Recovery settling time, +/-2% (s)','Peak DC-link voltage deviation (kV)', ...
       'Min-damping ratio zeta (nominal)','Maximum real part of poles (1/s)'};
b3 = [M3b.Ppre M3b.Vmin M3b.Ipk M3b.dPf M3b.osc M3b.ts M3b.dVdc zb0 rb0];
c3 = [M3c.Ppre M3c.Vmin M3c.Ipk M3c.dPf M3c.osc M3c.ts M3c.dVdc zc0 rc0];
T3 = s5_print(nm3, b3, c3);
fprintf('\nObservations:\n');
fprintf('1. Pre-fault steady state             : P = %.3f MW (baseline), %.3f MW (proposed) -> normal\n', M3b.Ppre, M3c.Ppre);
fprintf('2. Voltage / current during fault     : Vpcc,min = %.3f / %.3f pu ; Ipk = %.1f / %.1f A (baseline / proposed)\n', ...
    M3b.Vmin, M3c.Vmin, M3b.Ipk, M3c.Ipk);
fprintf('3. Active-power oscillation (fault)   : peak |dP| = %s / %s MW\n', s4_fmt(M3b.dPf), s4_fmt(M3c.dPf));
fprintf('4. Recovery after clearing            : settling %s s (baseline) / %s s (proposed)\n', s4_fmt(M3b.ts), s4_fmt(M3c.ts));
% FIXED: neutral wording; the old text always said "improvement".
fprintf('5. Post-fault amplitude / min zeta    : %s -> %s MW ; %s -> %s  (%s)\n', ...
    s4_fmt(M3b.osc), s4_fmt(M3c.osc), s4_fmt(zb0), s4_fmt(zc0), ...
    ternary(M3c.osc <= M3b.osc*1.5 && zc0 >= 0.10, 'comparable, damping requirement kept', 'proposed is worse - check design'));
fprintf('   (at SCR = 5 the baseline is already well damped; the benefit of the controller is at weak grid, see Test 5)\n');
s5_xl(xlsx, T3, 'T3_Fault_Metrics');
s5_xl(xlsx, s5_datatab(t3, E3b, E3c, 1:5:n3), 'T3_Data');
Tfp = table(Rf_pu, Rf_pu*Zb, tf_, tc_, Vres, dip, 'VariableNames', {'Rf_pu','Rf_ohm','t_fault_s','t_clear_s','Vresidual_pu','SagDepth_pu'});
s5_xl(xlsx, Tfp, 'T3_Fault_Parameters');

%% TEST 4 : FFT-BASED OSCILLATION COMPARISON
fprintf('\n=========================================================\n');
fprintf('  TEST 4 : FFT-BASED OSCILLATION COMPARISON\n');
fprintf('=========================================================\n');
zz = zol;  zz(rol >= 0 | isnan(zol)) = Inf;
[zmin4, i4] = min(zz(:));
if isinf(zmin4), is4 = is_r;  iv4 = iv_r; else, [is4, iv4] = ind2sub(size(zz), i4); end
G4 = Gfull{is4,iv4};  P4 = P_list(iv4);
[B4, C4] = s5_sys(G4, Kc, yname, uname);
[zb4, fb4, rb4] = s3_modes(pole(B4));  [zc4, ~, rc4] = s3_modes(pole(C4));
f_osc = fol(is4,iv4);
fprintf('Analysis case          : SCR = %.1f, wind = %.1f m/s (least-damped stable baseline)\n', SCR_list(is4), v_list(iv4));
fprintf('Identified f_osc (Step 2 model, dq) : %.2f Hz  -> abc %.2f / %.2f Hz\n', f_osc, abs(f_osc - p.f_grid), f_osc + p.f_grid);

fs4 = round(1/dt);  t4f = (0:dt:12).';  n4f = numel(t4f);
rng(2);
U4 = zeros(n4f,4);  U4(:,3:4) = 0.005*randn(n4f,2);
Y4b = s5_run(B4, U4, t4f);  Y4c = s5_run(C4, U4, t4f);
w4 = t4f >= 2;  tw4 = t4f(w4);  w0 = 2*pi*p.f_grid;  Ib = Sb/(sqrt(3)*Vb);
dPb = Sb/1e6*Y4b(w4,1);  dPc = Sb/1e6*Y4c(w4,1);
dPb = dPb - mean(dPb);   dPc = dPc - mean(dPc);
iab = sqrt(2)*Ib*((P4 + Y4b(w4,2)).*cos(w0*tw4) - Y4b(w4,3).*sin(w0*tw4));
iac = sqrt(2)*Ib*((P4 + Y4c(w4,2)).*cos(w0*tw4) - Y4c(w4,3).*sin(w0*tw4));
Nseg = 2^14;
[fa, Ab]  = sso_spectrum(dPb, fs4, Nseg);   [~, Ac]  = sso_spectrum(dPc, fs4, Nseg);
[fi, Aib] = sso_spectrum(iab, fs4, Nseg);   [~, Aic] = sso_spectrum(iac, fs4, Nseg);
% FIXED: the peak is now found separately for each case (before, the proposed
% spectrum was read at the baseline peak bin, which gave a misleading "51 %
% reduction" although the damping ratio had gone down).
sel = fa >= 1 & fa <= 200;
fsel = fa(sel);
Ab_s = Ab(sel);  Ac_s = Ac(sel);
[ampP_b, ib_] = max(Ab_s);  f_pk_b = fsel(ib_);
[ampP_c, ic_] = max(Ac_s);  f_pk_c = fsel(ic_);
f_peak = f_pk_b;
f_sub_b = abs(f_pk_b - p.f_grid);  f_sup_b = f_pk_b + p.f_grid;
f_sub_c = abs(f_pk_c - p.f_grid);  f_sup_c = f_pk_c + p.f_grid;
ampS_b = max(Aib(abs(fi - f_sub_b) <= 1));  ampU_b = max(Aib(abs(fi - f_sup_b) <= 1));
ampS_c = max(Aic(abs(fi - f_sub_c) <= 1));  ampU_c = max(Aic(abs(fi - f_sup_c) <= 1));
rmsb = sqrt(mean(dPb.^2));  rmsc = sqrt(mean(dPc.^2));

nm4 = {'f_osc identified, dq frame (Hz)','Dominant FFT peak of dP (Hz)','abc sub-synchronous component (Hz)', ...
       'abc super-synchronous component (Hz)','dP FFT amplitude at own peak (MW)','RMS of dP oscillation (MW)', ...
       'Phase-a current FFT, sub-sync (A)','Phase-a current FFT, super-sync (A)','Min-damping ratio zeta', ...
       'dP peak reduction (%)','dP RMS reduction (%)'};
b4 = [f_osc f_pk_b f_sub_b f_sup_b ampP_b rmsb ampS_b ampU_b zb4 NaN NaN];
c4 = [f_osc f_pk_c f_sub_c f_sup_c ampP_c rmsc ampS_c ampU_c zc4 100*(1 - ampP_c/ampP_b) 100*(1 - rmsc/rmsb)];
T4 = s5_print(nm4, b4, c4);
fprintf('Note: FFT magnitude alone does not prove stability; see closed-loop damping ratio / pole results above.\n');

nmA = 'Fig_T4a_FFT_Active_Power';
fg = figure('Name',nmA,'NumberTitle','off','Position',[120 80 950 540],'Color','w');
semilogy(fa, Ab + 1e-12, 'b', 'LineWidth', 2.2); hold on;
semilogy(fa, Ac + 1e-12, 'r--', 'LineWidth', 1.8);
yl4 = ylim;  plot(f_peak*[1 1], yl4, 'k:', 'LineWidth', 1.4);  hold off;
xlim([0 200]);  xlabel('Frequency (Hz, dq frame)');  ylabel('FFT magnitude (MW)');
title('FFT of active-power oscillation - baseline vs H_\infty controller');
legend('Without H_\infty controller','With H_\infty controller',sprintf('Baseline peak %.2f Hz', f_peak),'Location','northeast');
s5_fmt(gca);  s5_save(fg, outDir, nmA);

nmB = 'Fig_T4b_FFT_Phase_Current';
fg = figure('Name',nmB,'NumberTitle','off','Position',[120 80 950 540],'Color','w');
semilogy(fi, Aib + 1e-12, 'b', 'LineWidth', 2.2); hold on;
semilogy(fi, Aic + 1e-12, 'r--', 'LineWidth', 1.8);
yl4 = ylim;  plot(p.f_grid*[1 1], yl4, 'k:', 'LineWidth', 1.4);  hold off;
xlim([0 300]);  xlabel('Frequency (Hz)');  ylabel('FFT magnitude (A)');
title('FFT of phase-a grid current - baseline vs H_\infty controller');
legend('Without H_\infty controller','With H_\infty controller','50 Hz','Location','northeast');
s5_fmt(gca);  s5_save(fg, outDir, nmB);

s5_xl(xlsx, T4, 'T4_FFT_Metrics');
m300 = fi <= 300;
T4s = table(fi(m300), s5_fin(Ab(m300)), s5_fin(Ac(m300)), s5_fin(Aib(m300)), s5_fin(Aic(m300)), ...
    'VariableNames', {'Frequency_Hz','dP_base_MW','dP_prop_MW','ia_base_A','ia_prop_A'});
s5_xl(xlsx, T4s, 'T4_Spectrum');

%% TEST 5 : ROBUSTNESS ANALYSIS
fprintf('\n=========================================================\n');
fprintf('  TEST 5 : ROBUSTNESS (GRID STRENGTH / IMPEDANCE)\n');
fprintf('=========================================================\n');
cn    = {'Case 1: Strong grid','Case 2: Medium-strength grid','Case 3: Weak grid','Case 4: Increased grid impedance'};
scr_c = [10 5 3 3];
zf    = [1 1 1 1.5];
scr_e = scr_c./zf;
v_rob = 10;  P_rob = s5_pwind(p, v_rob);
t5 = (0:dt:3).';  n5 = numel(t5);
U5 = zeros(n5,4);  U5(t5 >= 0.5,1) = 0.1;  U5(t5 >= 1.5 & t5 < 1.6,3) = -0.05;
fprintf('Wind = %.0f m/s (P = %.3f pu). Disturbances: +0.10 pu P-ref step at 0.5 s, -0.05 pu grid-voltage sag 1.5-1.6 s\n', v_rob, P_rob);
fprintf('%-32s %7s %9s %10s | %9s %9s | %9s %9s | %8s %8s | %8s %8s | %-9s %-9s\n', 'Case','SCR_eff','Zg (pu)','Zg (Ohm)', ...
    'zeta_b','zeta_c','pkP_b MW','pkP_c MW','ts_b s','ts_c s','osc_b MW','osc_c MW','Base','Prop');
R5 = zeros(4,15);  st5 = cell(4,2);  P5 = zeros(n5,8);
for c = 1:4
    Gc = s3_plant(s3_params(p, P_rob, scr_e(c)));
    [Bc, Cc] = s5_sys(Gc, Kc, yname, uname);
    [zb, ~, rb] = s3_modes(pole(Gc));  [zc, ~, rc] = s3_modes(pole(Cc));
    Yb = s5_run(Bc, U5, t5);  Yc = s5_run(Cc, U5, t5);
    Pb = Sb/1e6*Yb(:,1);  Pc = Sb/1e6*Yc(:,1);
    Mb = s5_stepmet(t5, Pb, 0.5, 1.5, rb < 0);
    Mc = s5_stepmet(t5, Pc, 0.5, 1.5, rc < 0);
    ob = s5_tailosc(t5, Pb, 1.7, rb < 0);  oc = s5_tailosc(t5, Pc, 1.7, rc < 0);
    kb = ternary(rb < 0, max(abs(Pb)), Inf);  kc = ternary(rc < 0, max(abs(Pc)), Inf);
    st5{c,1} = s4_stat(rb, zb, 0.10);  st5{c,2} = s4_stat(rc, zc, 0.10);
    Zgp = 1/scr_e(c);
    R5(c,:) = [scr_c(c) zf(c) scr_e(c) Zgp Zgp*Zb zb zc rb rc kb kc Mb.ts Mc.ts ob oc];
    P5(:,2*c-1:2*c) = [Pb Pc];
    fprintf('%-32s %7.2f %9.3f %10.2f | %9s %9s | %9s %9s | %8s %8s | %8s %8s | %-9s %-9s\n', cn{c}, scr_e(c), Zgp, Zgp*Zb, ...
        s4_fmt(zb), s4_fmt(zc), s4_fmt(kb), s4_fmt(kc), s4_fmt(Mb.ts), s4_fmt(Mc.ts), s4_fmt(ob), s4_fmt(oc), st5{c,1}, st5{c,2});
    s5_plot(outDir, sprintf('Fig_T5%c_%s', 'a'+c-1, regexprep(strtok(cn{c},':'),'\s','')), t5, [Pb Pc], lb, ...
        'Time (s)', '\DeltaP (MW)', sprintf('%s (SCR_{eff} = %.2f)', cn{c}, scr_e(c)), [], true);
end

zcl5 = zeros(nS,1);  rcl5 = zeros(nS,1);
for is = 1:nS
    [zcl5(is), ~, rcl5(is)] = s3_modes(pole(s4_cl(Gfull{is,3}, Kc, yname, uname)));
end
nmE = 'Fig_T5e_Damping_vs_SCR';
fg = figure('Name',nmE,'NumberTitle','off','Position',[120 80 950 540],'Color','w');
plot(SCR_list, zol(:,3), 'b--o', 'LineWidth', 2); hold on;
plot(SCR_list, zcl5, 'r-s', 'LineWidth', 2);
plot([SCR_list(1) SCR_list(end)], 0.10*[1 1], 'k:', 'LineWidth', 1.6);  hold off;
xlabel('SCR');  ylabel('Worst damping ratio \zeta');
title(sprintf('Damping vs. grid strength (%.0f m/s)', v_list(3)));
legend('Baseline','Proposed H_\infty','\zeta_{req} = 0.10','Location','best');
s5_fmt(gca);  s5_save(fg, outDir, nmE);

fprintf('\nDamping vs SCR at %.0f m/s:\n%8s %12s %12s %12s\n', v_list(3), 'SCR','zeta_base','zeta_prop','status_prop');
for is = 1:nS
    fprintf('%8.1f %12s %12s %12s\n', SCR_list(is), s4_fmt(zol(is,3)), s4_fmt(zcl5(is)), s4_stat(rcl5(is), zcl5(is), 0.10));
end
T5 = table(cn(:), R5(:,1), R5(:,2), R5(:,3), R5(:,4), R5(:,5), s5_fin(R5(:,6)), s5_fin(R5(:,7)), s5_fin(R5(:,8)), s5_fin(R5(:,9)), ...
    s5_fin(R5(:,10)), s5_fin(R5(:,11)), s5_fin(R5(:,12)), s5_fin(R5(:,13)), s5_fin(R5(:,14)), s5_fin(R5(:,15)), st5(:,1), st5(:,2), ...
    'VariableNames', {'Case','SCR_nominal','ImpedanceFactor','SCR_effective','Zg_pu','Zg_ohm','Zeta_base','Zeta_prop', ...
    'MaxRe_base','MaxRe_prop','PeakdP_base_MW','PeakdP_prop_MW','Settling_base_s','Settling_prop_s','PostSagOsc_base_MW', ...
    'PostSagOsc_prop_MW','Status_base','Status_prop'});
s5_xl(xlsx, T5, 'T5_Robustness');
T5d = table(SCR_list(:), s5_fin(zol(:,3)), s5_fin(zcl5), s5_fin(rol(:,3)), s5_fin(rcl5), ...
    'VariableNames', {'SCR','Zeta_base','Zeta_prop','MaxRe_base','MaxRe_prop'});
s5_xl(xlsx, T5d, 'T5_Damping_vs_SCR');
T5t = array2table(s5_fin([t5(1:5:end) P5(1:5:end,:)]), 'VariableNames', {'Time_s','C1_base_MW','C1_prop_MW','C2_base_MW','C2_prop_MW', ...
    'C3_base_MW','C3_prop_MW','C4_base_MW','C4_prop_MW'});
s5_xl(xlsx, T5t, 'T5_Data');

try
    T1_s = addvars(T1, repmat({'Test 1'}, height(T1), 1), 'Before', 1, 'NewVariableNames', 'Test');
    T3_s = addvars(T3, repmat({'Test 3'}, height(T3), 1), 'Before', 1, 'NewVariableNames', 'Test');
    T4_s = addvars(T4, repmat({'Test 4'}, height(T4), 1), 'Before', 1, 'NewVariableNames', 'Test');
    TS = [T1_s; T3_s; T4_s];
    s5_xl(xlsx, TS, 'Summary');
catch
end

nUb = nnz(rol >= 0);  nUc = 0;
for is = 1:nS
    for iv = 1:nV
        [~, ~, rr] = s3_modes(pole(s4_cl(Gfull{is,iv}, Kc, yname, uname)));
        if rr >= 0, nUc = nUc + 1; end
    end
end
fprintf('\n=========================================================\n');
fprintf('                 STEP 5 : MAIN RESULT\n');
fprintf('=========================================================\n');
fprintf('Test 1 normal operation  : max |dP| baseline vs proposed = %.4g MW\n', max(abs(E1b.P - E1c.P)));
fprintf('Test 2 wind variation    : mean oscillation amplitude  %s -> %s MW\n', s4_fmt(mean(R2(:,7))), s4_fmt(mean(R2(:,8))));
fprintf('Test 3 three-phase fault : post-fault amplitude        %s -> %s MW ; recovery %s -> %s s\n', ...
    s4_fmt(M3b.osc), s4_fmt(M3c.osc), s4_fmt(M3b.ts), s4_fmt(M3c.ts));
fprintf('Test 4 FFT               : dP peak %.2f Hz (%s MW) -> %.2f Hz (%s MW) ; RMS change %.1f %%\n', ...
    f_pk_b, s4_fmt(ampP_b), f_pk_c, s4_fmt(ampP_c), 100*(rmsc/rmsb - 1));
fprintf('Test 5 robustness        : unstable conditions over %d-point grid: baseline %d -> proposed %d\n', nS*nV, nUb, nUc);
fprintf('Plots saved in folder    : %s\n', fullfile(pwd, outDir));
fprintf('Excel results file       : %s\n', fullfile(pwd, xlsx));
fprintf('=========================================================\n');

%% ================================ FUNCTIONS ==============================
function P = s5_pwind(p, v)
Pa = p.N_turbines*0.5*p.rho_air*pi*p.R_turbine^2*p.Cp_max*v.^3;
P  = min(Pa, p.P_HVDC_rated)/p.P_HVDC_rated;
end

function [B, C] = s5_sys(Gf, Kc, y, u)
B = s4_base(Gf);
C = s4_cl(Gf, Kc, y, u);
end

function Y = s5_run(sys, U, t)
Y = lsim(sys, U, t, [], 'zoh');
if size(Y,2) < 6, Y(:,6) = 0; end
end

function E = s5_elec(Y, t, P0, p, Gdc)
Sb = p.P_HVDC_rated;  Vb = p.Vgrid_LL;  Ib = Sb/(sqrt(3)*Vb);  w0 = 2*pi*p.f_grid;
vd = 1 + Y(:,4);  vq = Y(:,5);  id = P0 + Y(:,2);  iq = Y(:,3);
E.P   = Sb/1e6*(P0 + Y(:,1));
E.Q   = Sb/1e6*(P0*Y(:,5) - Y(:,3));
E.va  = sqrt(2/3)*Vb/1e3*(vd.*cos(w0*t) - vq.*sin(w0*t));
E.ia  = sqrt(2)*Ib*(id.*cos(w0*t) - iq.*sin(w0*t));
E.vpu = sqrt(vd.^2 + vq.^2);
E.vpk = sqrt(2/3)*Vb/1e3*E.vpu;
E.ipk = sqrt(2)*Ib*sqrt(id.^2 + iq.^2);
if all(isfinite(Y(:,1)))
    E.Vdc = p.Vdc_HVDC/1e3*(1 - lsim(Gdc, Y(:,1), t));
else
    E.Vdc = NaN(size(t));
end
end

function M = s5_stepmet(t, y, ts_, te_, stab)
idx = t >= ts_ & t < te_;  tt = t(idx);  yy = y(idx);
yf = mean(y(t >= te_ - 0.3 & t < te_));
y0 = mean(y(t >= ts_ - 0.3 & t < ts_));
d  = yf - y0;
M.step = d;
M.os   = max(0, max(sign(d)*(yy - yf)));
M.osc  = max(abs(y(t >= ts_ + 0.5 & t < te_) - yf));
band   = max(0.02*abs(d), 1e-3);
k = find(abs(yy - yf) > band, 1, 'last');
if isempty(k), M.ts = 0; elseif k >= numel(tt), M.ts = Inf; else, M.ts = tt(k+1) - ts_; end
if ~stab || ~all(isfinite(yy)), M.os = Inf; M.osc = Inf; M.ts = Inf; end
end

function o = s5_tailosc(t, y, t0, stab)
yf = mean(y(t >= t(end) - 0.3));
o  = max(abs(y(t >= t0) - yf));
if ~stab || ~all(isfinite(y)), o = Inf; end
end

function M = s5_fmet(E, t, tf, tc, stab, band)
pre = t >= tf - 0.5 & t < tf;  ff = t >= tf & t < tc;
M.Ppre = mean(E.P(pre));
M.Vmin = min(E.vpu(ff));
M.Ipk  = max(E.ipk(ff));
d = E.P - M.Ppre;
M.dPf  = max(abs(d(ff)));
M.osc  = max(abs(d(t >= tc + 0.3)));
M.dVdc = max(abs(E.Vdc - E.Vdc(1)));
k = find(abs(d) > band & t >= tc, 1, 'last');
if isempty(k), M.ts = 0; elseif k >= numel(t), M.ts = Inf; else, M.ts = t(k+1) - tc; end
if ~stab || ~all(isfinite(d)), M.osc = Inf; M.ts = Inf; end
end

function s5_plot(outDir, name, t, Y, lbl, xl, yl, ttl, xlm, clip)
if nargin < 9,  xlm = [];  end
if nargin < 10, clip = false;  end
col = {'b','r','k','m'};  sty = {'-','--',':','-.'};  lw = [2.2 1.8 1.6 1.4];
fig = figure('Name',name,'NumberTitle','off','Position',[120 80 950 540],'Color','w');
hold on;
for k = 1:size(Y,2)
    plot(t, Y(:,k), 'Color', col{k}, 'LineStyle', sty{k}, 'LineWidth', lw(k));
end
hold off;
xlabel(xl);  ylabel(yl);  title(ttl);
if ~isempty(lbl), legend(lbl, 'Location', 'best'); end
if isempty(xlm), xlim([t(1) t(end)]); else, xlim(xlm); end
if clip && size(Y,2) >= 2
    lim = 1.5*max(abs(Y(isfinite(Y(:,2)),2))) + 1e-9;
    if max(abs(Y(isfinite(Y(:,1)),1))) > 10*lim, ylim([-lim lim]); end
end
s5_fmt(gca);
s5_save(fig, outDir, name);
end

function s5_fmt(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 22, 'FontWeight', 'bold');
set(get(ax, 'XLabel'), 'FontName', 'Times New Roman', 'FontSize', 22, 'FontWeight', 'bold');
set(get(ax, 'YLabel'), 'FontName', 'Times New Roman', 'FontSize', 22, 'FontWeight', 'bold');
set(get(ax, 'Title'), 'FontName', 'Times New Roman', 'FontSize', 22, 'FontWeight', 'bold');
lg = findobj(get(ax, 'Parent'), 'Type', 'legend');
if ~isempty(lg)
    set(lg, 'FontName', 'Times New Roman', 'FontSize', 22, 'FontWeight', 'bold');
end
end

function s5_save(fig, outDir, name)
if ~exist(outDir, 'dir'), mkdir(outDir); end
base = fullfile(pwd, outDir, regexprep(name, '[^\w]', '_'));
if isempty(fig) || ~ishghandle(fig)
    fig = findall(0, 'Type', 'figure', 'Name', name);
    if isempty(fig)
        return;
    end
    fig = fig(1);
end
drawnow;
try
    print(fig, [base '.png'], '-dpng', '-r300');
catch
    try
        saveas(fig, [base '.png']);
    catch
    end
end
try
    savefig(fig, [base '.fig']);
catch
end
if ishghandle(fig), close(fig); end
end

function T = s5_print(names, b, c)
fprintf('%-58s %16s %16s\n', 'Parameter', 'Baseline', 'Proposed H-inf');
for k = 1:numel(names)
    fprintf('%-58s %16s %16s\n', names{k}, s4_fmt(b(k)), s4_fmt(c(k)));
end
T = table(names(:), s5_fin(b(:)), s5_fin(c(:)), 'VariableNames', {'Parameter','Baseline','Proposed'});
end

function x = s5_fin(x)
x(isinf(x)) = NaN;
end

function T = s5_datatab(t, Eb, Ec, idx)
T = table(t(idx), s5_fin(Eb.P(idx)), s5_fin(Ec.P(idx)), s5_fin(Eb.Q(idx)), s5_fin(Ec.Q(idx)), ...
    s5_fin(Eb.Vdc(idx)), s5_fin(Ec.Vdc(idx)), s5_fin(Eb.vpu(idx)), s5_fin(Ec.vpu(idx)), ...
    s5_fin(Eb.ipk(idx)), s5_fin(Ec.ipk(idx)), 'VariableNames', ...
    {'Time_s','P_base_MW','P_prop_MW','Q_base_MVAr','Q_prop_MVAr','Vdc_base_kV','Vdc_prop_kV', ...
     'Vpcc_base_pu','Vpcc_prop_pu','Ipk_base_A','Ipk_prop_A'});
end

function s5_xl(file, T, sheet)
try
    writetable(T, file, 'Sheet', sheet);
catch ME
    fprintf('Excel write failed for sheet %s: %s\n', sheet, ME.message);
end
end

function c = s5_sumrows(label, T)
n = height(T);
c = cell(n,4);
c(:,1) = {label};
c(:,2) = T.Parameter;
c(:,3) = num2cell(T.Baseline);
c(:,4) = num2cell(T.Proposed);
end

function m = sso_params(p, P_pu, SCR)
Sb = p.P_HVDC_rated;  Vb = p.Vgrid_LL;  Zb = Vb^2 / Sb;
m.w0  = 2*pi*p.f_grid;
m.R   = (p.Rf_grid + p.Rarm/2) / Zb;
m.X   = m.w0*(p.Lf_grid + p.Larm/2) / Zb + p.Xtr_pu;
Zg    = Vb^2 / (SCR * Sb);
Rg    = Zg / sqrt(1 + p.X_R_ratio^2);
m.Rg  = Rg / Zb;
m.Xg  = p.X_R_ratio * Rg / Zb;
m.Bc  = m.w0 * p.Cf_grid * Zb;
m.P   = P_pu;
m.Kp  = p.Kp_current;   m.Ki  = p.Ki_current;
m.Kpp = p.Kp_pll;       m.Kip = p.Ki_pll;
m.Td  = p.Td_ctrl;
end

function [Yc, Znet] = sso_eval(m, s)
J  = [0 -1; 1 0];   I2 = eye(2);
Gi  = (m.Kp + m.Ki/s) * exp(-s*m.Td);
den = m.X*s/m.w0 + m.R + Gi;
Tp  = (m.Kpp*s + m.Kip) / (s^2 + m.Kpp*s + m.Kip);
i0  = [m.P; 0];
vc0 = [1 + m.R*m.P; m.X*m.P];
Ac  = Gi*I2 - m.X*J;
psi = Ac*(J*i0) + J*vc0;
Yc  = (I2 - psi*Tp*[0 1]) / den;
Zg  = (m.Rg + m.Xg*s/m.w0)*I2 + m.Xg*J;
Ycf = (m.Bc*s/m.w0)*I2 + m.Bc*J;
Znet = inv(inv(Zg) + Ycf);
end

function r = sso_scan(m, f_dq)
nF = numel(f_dq);
lam = zeros(nF,2);  dtl = zeros(nF,1);  Msv = zeros(nF,1);
Zc_pos = zeros(nF,1);  Zn_pos = zeros(nF,1);
for k = 1:nF
    s = 1j*2*pi*f_dq(k);
    [Yc, Zn] = sso_eval(m, s);
    L = Zn*Yc;
    e = eig(L);
    if k > 1
        if abs(e(1)-lam(k-1,1)) + abs(e(2)-lam(k-1,2)) > ...
           abs(e(2)-lam(k-1,1)) + abs(e(1)-lam(k-1,2))
            e = e([2 1]);
        end
    end
    lam(k,:) = e.';
    M = eye(2) + L;
    dtl(k) = det(M);
    Msv(k) = norm(inv(M));
    Zc = inv(Yc);
    Zc_pos(k) = 0.5*(Zc(1,1)+Zc(2,2)) + 0.5j*(Zc(2,1)-Zc(1,2));
    Zn_pos(k) = 0.5*(Zn(1,1)+Zn(2,2)) + 0.5j*(Zn(2,1)-Zn(1,2));
end
r.f = f_dq;  r.lam = lam;  r.Zc_pos = Zc_pos;  r.Zn_pos = Zn_pos;

ph  = unwrap(angle(dtl));
r.Zr = max(0, round(-(ph(end) - ph(1)) / pi));

d = abs(1 + lam);
[dm, idx] = min(d(:));
[ik, ~] = ind2sub(size(d), idx);
r.dmin = dm;  r.fdmin = f_dq(ik);

r.PM = NaN;  r.fPM = NaN;
for j = 1:2
    mg = abs(lam(:,j));
    ic = find((mg(1:end-1)-1) .* (mg(2:end)-1) < 0);
    for q = 1:numel(ic)
        i = ic(q);
        a  = (1 - mg(i)) / (mg(i+1) - mg(i));
        lc = lam(i,j) + a*(lam(i+1,j) - lam(i,j));
        pm = 180 - abs(rad2deg(angle(lc)));
        if isnan(r.PM) || pm < r.PM
            r.PM  = pm;
            r.fPM = f_dq(i) + a*(f_dq(i+1) - f_dq(i));
        end
    end
end
[r.Ms, iM] = max(Msv);  r.fMs = f_dq(iM);
end

function c = sso_class(r)
if r.Zr > 0
    c = 'UNSTABLE';
elseif r.dmin < 0.30 || (~isnan(r.PM) && r.PM < 15)
    c = 'VERY WEAK';
elseif r.dmin < 0.50 || (~isnan(r.PM) && r.PM < 30)
    c = 'WEAK';
else
    c = 'OK';
end
end

function [did, diq] = sso_noise_response(m, nd, nq, fs)
N  = numel(nd);
Vd = fft(nd);  Vq = fft(nq);
Yd = zeros(N/2+1, 1);  Yq = Yd;
for k = 2:N/2+1
    f = (k-1) * fs / N;
    s = 1j*2*pi*f;
    [Yc, Zn] = sso_eval(m, s);
    H = -((eye(2) + Yc*Zn) \ Yc);
    y = H * [Vd(k); Vq(k)];
    Yd(k) = y(1);  Yq(k) = y(2);
end
Yd(end) = real(Yd(end));  Yq(end) = real(Yq(end));
Ydf = [Yd; conj(Yd(end-1:-1:2))];
Yqf = [Yq; conj(Yq(end-1:-1:2))];
did = real(ifft(Ydf));
diq = real(ifft(Yqf));
end

function [f, A] = sso_spectrum(x, fs, Nseg)
w   = 0.5*(1 - cos(2*pi*(0:Nseg-1).'/Nseg));
hop = Nseg/2;
nS  = floor((numel(x) - Nseg)/hop) + 1;
P   = zeros(Nseg, 1);
for q = 1:nS
    xs = x((q-1)*hop + (1:Nseg));
    P  = P + abs(fft(xs .* w)).^2;
end
P = P / nS;
A = 2*sqrt(P(1:Nseg/2+1)) / sum(w);
f = (0:Nseg/2).' * fs / Nseg;
end

function sso_fmt(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, 'FontWeight', 'bold');
end

function sso_fmt3(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, 'FontWeight', 'bold');
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function m = s3_params(p, P_pu, SCR)
Sb = p.P_HVDC_rated;  Vb = p.Vgrid_LL;  Zb = Vb^2/Sb;
m.w0  = 2*pi*p.f_grid;
m.R   = (p.Rf_grid + p.Rarm/2)/Zb;
m.X   = m.w0*(p.Lf_grid + p.Larm/2)/Zb + p.Xtr_pu;
Zg    = Vb^2/(SCR*Sb);
Rg    = Zg/sqrt(1 + p.X_R_ratio^2);
m.Rg  = Rg/Zb;
m.Xg  = p.X_R_ratio*Rg/Zb;
m.Bc  = m.w0*p.Cf_grid*Zb;
m.P   = P_pu;
m.Kp  = p.Kp_current;  m.Ki = p.Ki_current;
m.Kpp = p.Kp_pll;      m.Kip = p.Ki_pll;
m.Td  = p.Td_ctrl;
end

function G = s3_plant(m)
w0  = m.w0;
sn  = tf('s');
Ki  = m.Ki/w0;   Tdn = m.Td*w0;
Kpp = m.Kpp/w0;  Kip = m.Kip/w0^2;

Pd   = (1 - sn*Tdn/2 + (sn*Tdn)^2/12)/(1 + sn*Tdn/2 + (sn*Tdn)^2/12);
Gi   = (m.Kp + Ki/sn)*Pd;
den  = m.X*sn + m.R + Gi;
Tp   = (Kpp*sn + Kip)/(sn^2 + Kpp*sn + Kip);
psiq = m.P*(Gi + m.R) + 1;

Gr_s = ss(Gi/den);  yd_s = ss(1/den);  yq_s = ss((1 - psiq*Tp)/den);  Z0 = ss(0);
Conv = [Gr_s, Z0, -yd_s, Z0; Z0, Gr_s, Z0, -yq_s];
Conv.InputName  = {'ud','uq','vd','vq'};
Conv.OutputName = {'id','iq'};

Rg = m.Rg;  Xg = m.Xg;  Bc = m.Bc;
An = [-Rg/Xg,  1,       1/Xg,  0;
      -1,     -Rg/Xg,   0,     1/Xg;
      -1/Bc,   0,       0,     1;
       0,     -1/Bc,   -1,     0];
Bn = [0 0 -1/Xg 0; 0 0 0 -1/Xg; 1/Bc 0 0 0; 0 1/Bc 0 0];
Cn = [0 0 1 0; 0 0 0 1];
Net = ss(An, Bn, Cn, zeros(2,4));
Net.InputName  = {'id','iq','vgd','vgq'};
Net.OutputName = {'vd','vq'};

Pb = ss([1 0 m.P 0]);
Pb.InputName  = {'id','iq','vd','vq'};
Pb.OutputName = {'P'};

Gn = connect(Conv, Net, Pb, {'ud','uq','vgd','vgq'}, {'P','id','iq','vd','vq'});
[A, B, C, D] = ssdata(Gn);
G = ss(w0*A, w0*B, C, D);
G.InputName  = {'ud','uq','vgd','vgq'};
G.OutputName = {'P','id','iq','vd','vq'};
end

function [zmin, fz, rmax] = s3_modes(lam)
rmax = max(real(lam));
lam_c = lam(imag(lam) > 1e-4 & abs(lam) < 2*pi*400);
if ~isempty(lam_c)
    f = imag(lam_c)/(2*pi);
    z = -real(lam_c)./abs(lam_c);
    sel = f >= 0.1 & f <= 40;
    if any(sel)
        zs = z(sel); fs = f(sel);
        [zmin, i] = min(zs); fz = fs(i);
    else
        [zmin, i] = min(z); fz = f(i);
    end
else
    zmin = 1; fz = 0;
end
if rmax >= 0 && zmin > 0
    zmin = -zmin;
end
end

function S = s3_strip(sys)
S = ss(sys.A, sys.B, sys.C, sys.D);
end

function W1 = s3_W1(wb, Ms, epsl)
s  = tf('s');
W1 = ss((s/Ms + wb)/(s + wb*epsl));
end

function W2 = s3_W2(Kb, Kh, fz)
wz = 2*pi*fz;
W2 = ss((1/Kb)*tf([1/wz 1], [Kh/(wz*Kb) 1]));
end

function [K, gam] = s3_synth(Gn, W1, W2, W3)
Gp = s3_strip(Gn);
Gp.InputName  = {'u'};
Gp.OutputName = {'y'};
W1.InputName = 'e';  W1.OutputName = 'z1';
W2.InputName = 'u';  W2.OutputName = 'z2';
W3.InputName = 'y';  W3.OutputName = 'z3';
Sm = sumblk('e = r - y');
Pg = connect(Gp, W1, W2, W3, Sm, {'r','u'}, {'z1','z2','z3','e'});
try
    [K, ~, gam] = hinfsyn(Pg, 1, 1);
catch ME
    K = [];  gam = Inf;
end
end

function [score, zcl, rcl] = s3_score(K, Gfull, yname, uname)
Ks  = s3_strip(K);
nS  = size(Gfull,1);  nV = size(Gfull,2);
zcl = zeros(nS, nV);  rcl = zeros(nS, nV);
for is = 1:nS
    for iv = 1:nV
        try
            CL = s4_cl(Gfull{is,iv}, Ks, yname, uname);
            [zcl(is,iv), ~, rcl(is,iv)] = s3_modes(pole(CL));
        catch
            zcl(is,iv) = -1;  rcl(is,iv) = 1;
        end
    end
end
score = 1000*nnz(rcl < 0) + min(zcl(:));
end

function m = s3_mag(sys, f)
H = freqresp(sys, 2*pi*f(:));
m = abs(squeeze(H));
m = m(:);
end

function ph = s3_ph(sys, f)
H = freqresp(sys, 2*pi*f(:));
ph = unwrap(angle(squeeze(H)));
ph = ph(:)*180/pi;
end

function s3_block(ax, x, y, w, h, lines, fc)
rectangle('Parent', ax, 'Position', [x y w h], 'FaceColor', fc, 'LineWidth', 1.6, ...
    'Curvature', 0.05);
n = numel(lines);
for i = 1:n
    fw = 'normal';  if i == 1, fw = 'bold'; end
    text(ax, x + w/2, y + h - (i - 0.5)*h/(n + 0.3), lines{i}, ...
        'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', ...
        'FontWeight', fw, 'FontSize', 11.5);
end
end

function s3_arrow(ax, x1, y1, x2, y2)
plot(ax, [x1 x2], [y1 y2], 'k', 'LineWidth', 1.6);
d = [x2 - x1, y2 - y1];  d = d/norm(d);  nn = [-d(2) d(1)];
hs = 0.22;  tip = [x2 y2];
b1 = tip - hs*d + 0.45*hs*nn;  b2 = tip - hs*d - 0.45*hs*nn;
patch('Parent', ax, 'XData', [tip(1) b1(1) b2(1)], 'YData', [tip(2) b1(2) b2(2)], ...
    'FaceColor', 'k', 'EdgeColor', 'k');
end

function CL = s4_cl(Gf, Kc, yname, uname)
Gp = Gf;
Kc.InputName  = {yname};
Kc.OutputName = {'uk'};
if strcmp(uname, 'ud')
    Gp.InputName = {'ud','rq','vgd','vgq'};
    Sm = sumblk('ud = rd - uk');
else
    Gp.InputName = {'rd','uq','vgd','vgq'};
    Sm = sumblk('uq = rq - uk');
end
Gp.OutputName = {'P','id','iq','vd','vq'};
CL = connect(Gp, Kc, Sm, {'rd','rq','vgd','vgq'}, {'P','id','iq','vd','vq','uk'});
end

function B = s4_base(Gf)
B = Gf;
B.InputName  = {'rd','rq','vgd','vgq'};
B.OutputName = {'P','id','iq','vd','vq'};
end

function M = s4_metrics(t, y, t0, stable)
idx = t >= t0;
tt  = t(idx);
yy  = y(idx);
if stable
    yss = mean(y(t >= t(end) - 0.5));
else
    yss = 0;
end
dev = yy - yss;
M.peak = max(abs(dev));
M.rms  = sqrt(mean(dev.^2));
M.os   = 100*max(0, max(sign(yss)*yy) - abs(yss))/max(abs(yss), eps);
band   = max(0.02*M.peak, 1e-9);
if stable
    k = find(abs(dev) > band, 1, 'last');
    if isempty(k)
        M.ts = 0;
    elseif k >= numel(tt)
        M.ts = Inf;
    else
        M.ts = tt(k+1) - t0;
    end
else
    M.ts = Inf;
end
end

function s = s4_fmt(x)
if isnan(x)
    s = 'n/a';
elseif isinf(x)
    s = 'Inf';
else
    s = sprintf('%.4g', x);
end
end

function s = s4_stat(r, z, zr)
if r >= 0
    s = 'UNSTABLE';
elseif ~isnan(z) && z < zr
    s = 'LOW ZETA';
else
    s = 'OK';
end
end