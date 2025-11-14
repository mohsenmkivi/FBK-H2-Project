fprintf('=== LOADING SIMULINK DATA ===\n');

t   = out.Ia_pcc_no_filter.time(:);
ia  = out.Ia_pcc_no_filter.signals.values(:);
ib  = out.Ib_pcc_no_filter.signals.values(:);
ic  = out.Ic_pcc_no_filter.signals.values(:);

vab = out.Vab_pcc_no_filter.signals.values(:);
vbc = out.Vbc_pcc_no_filter.signals.values(:);
vca = out.Vca_pcc_no_filter.signals.values(:);

N_samples = numel(t);
signal_lengths = [numel(t), numel(ia), numel(ib), numel(ic), numel(vab), numel(vbc), numel(vca)];
signal_names = {'t', 'ia', 'ib', 'ic', 'vab', 'vbc', 'vca'};

if any(signal_lengths ~= N_samples)
    error('Signal length mismatch detected!\n%s', sprintf('%s=%d ', signal_names{:}, signal_lengths));
end

duration_total = t(end) - t(1);
fprintf('✓ Loaded %d samples, Duration: %.4f s\n', N_samples, duration_total);

va = (2*vab + vbc) / 3;
vb = (-vab + vbc) / 3;
vc = -(vab + 2*vbc) / 3;

fprintf('✓ Phase voltages reconstructed\n');

f1 = 50;
dt = diff(t);
Fs = 1 / mean(dt);
Fs_std = std(1./dt);

if Fs_std/Fs > 0.01
    warning('Variable time step detected (%.2f%% variation). Using mean Fs=%.1f Hz', 100*Fs_std/Fs, Fs);
else
    fprintf('✓ Sampling frequency: Fs = %.1f Hz\n', Fs);
end

samples_per_cycle = Fs / f1;
fprintf('✓ Samples per cycle: %.1f\n', samples_per_cycle);

fprintf('\n=== SELECTING ANALYSIS WINDOW ===\n');

target_cycles = 40;
target_samples = round(target_cycles * samples_per_cycle);

if target_samples > N_samples
    actual_cycles = floor(N_samples / samples_per_cycle);
    actual_samples = round(actual_cycles * samples_per_cycle);
    warning('Only %.1f cycles available. Using %d cycles instead of %d.', N_samples/samples_per_cycle, actual_cycles, target_cycles);
else
    actual_cycles = target_cycles;
    actual_samples = target_samples;
end

idx_start = N_samples - actual_samples + 1;
idx_end = N_samples;
idx = idx_start:idx_end;

t_win  = t(idx);
va_win = va(idx);
vb_win = vb(idx);
vc_win = vc(idx);
ia_win = ia(idx);
ib_win = ib(idx);
ic_win = ic(idx);

duration_win = t_win(end) - t_win(1);

fprintf('✓ Window: samples %d to %d (%d total)\n', idx_start, idx_end, numel(idx));
fprintf('✓ Analyzing: %d cycles (%.4f s)\n', actual_cycles, duration_win);
fprintf('✓ Window position: last %.1f%% of signal\n', 100*numel(idx)/N_samples);

fprintf('\n=== CALCULATING RMS VALUES ===\n');

Vrms_a = rms(va_win);
Vrms_b = rms(vb_win);
Vrms_c = rms(vc_win);

Irms_a = rms(ia_win);
Irms_b = rms(ib_win);
Irms_c = rms(ic_win);

fprintf('Phase A: V_rms = %.2f V, I_rms = %.2f A\n', Vrms_a, Irms_a);
fprintf('Phase B: V_rms = %.2f V, I_rms = %.2f A\n', Vrms_b, Irms_b);
fprintf('Phase C: V_rms = %.2f V, I_rms = %.2f A\n', Vrms_c, Irms_c);

fprintf('\n=== CALCULATING REAL POWER ===\n');

p_inst = va_win.*ia_win + vb_win.*ib_win + vc_win.*ic_win;
P = mean(p_inst);

fprintf('✓ Real Power P = %.2f W (%.4f kW)\n', P, P*1e-3);

fprintf('\n=== CALCULATING APPARENT POWER ===\n');

S = Vrms_a*Irms_a + Vrms_b*Irms_b + Vrms_c*Irms_c;

fprintf('✓ Apparent Power S = %.2f VA (%.4f kVA)\n', S, S*1e-3);

PF_true = P / S;
fprintf('✓ True Power Factor = %.4f\n', PF_true);

fprintf('\n=== FUNDAMENTAL COMPONENT ANALYSIS ===\n');

function X_rms = goertzel_rms(x, Fs, f_target)
    N = numel(x);
    k = round(N * f_target / Fs);
    w = 2*pi*k/N;
    coeff = 2*cos(w);
    s_prev2 = 0;
    s_prev1 = 0;
    for n = 1:N
        s = x(n) + coeff*s_prev1 - s_prev2;
        s_prev2 = s_prev1;
        s_prev1 = s;
    end
    X_rms = (s_prev1 - exp(-1j*w)*s_prev2) / N * sqrt(2);
end

Va1 = goertzel_rms(va_win, Fs, f1);
Vb1 = goertzel_rms(vb_win, Fs, f1);
Vc1 = goertzel_rms(vc_win, Fs, f1);

Ia1 = goertzel_rms(ia_win, Fs, f1);
Ib1 = goertzel_rms(ib_win, Fs, f1);
Ic1 = goertzel_rms(ic_win, Fs, f1);

Pa1 = real(Va1 * conj(Ia1));
Pb1 = real(Vb1 * conj(Ib1));
Pc1 = real(Vc1 * conj(Ic1));

Qa1 = imag(Va1 * conj(Ia1));
Qb1 = imag(Vb1 * conj(Ib1));
Qc1 = imag(Vc1 * conj(Ic1));

P1 = Pa1 + Pb1 + Pc1;
Q1 = Qa1 + Qb1 + Qc1;
S1 = sqrt(P1^2 + Q1^2);
DPF = P1 / S1;

fprintf('✓ Fundamental Active Power P1 = %.2f W (%.4f kW)\n', P1, P1*1e-3);
fprintf('✓ Fundamental Reactive Power Q1 = %.2f var (%.4f kvar)\n', Q1, Q1*1e-3);
fprintf('✓ Displacement Power Factor DPF = %.4f\n', DPF);

fprintf('\n=== THD ANALYSIS (Phase A Current) ===\n');

hann_window = 0.5 * (1 - cos(2*pi*(0:numel(ia_win)-1)/(numel(ia_win)-1)))';
ia_windowed = ia_win .* hann_window;

N_fft = numel(ia_windowed);
Ia_fft = fft(ia_windowed) / N_fft;
f_axis = (0:N_fft-1)' * Fs / N_fft;

Ia_mag = abs(Ia_fft) * sqrt(2);
harmonic_order = f_axis / f1;

max_harmonic = 50;
Irms_harmonic = zeros(max_harmonic, 1);

for h = 1:max_harmonic
    bin_mask = (harmonic_order >= (h-0.5)) & (harmonic_order < (h+0.5));
    Irms_harmonic(h) = norm(Ia_mag(bin_mask));
end

I1_rms = Irms_harmonic(1);
THDi = sqrt(sum(Irms_harmonic(2:end).^2)) / I1_rms;

fprintf('✓ Fundamental Current (h=1): %.4f A\n', I1_rms);
fprintf('✓ THDi (Phase A): %.2f %%\n', THDi*100);

[harmonics_sorted, harm_idx] = sort(Irms_harmonic(2:end), 'descend');
fprintf('\nTop 5 Harmonics:\n');
for i = 1:min(5, sum(harmonics_sorted > I1_rms*0.01))
    h = harm_idx(i) + 1;
    fprintf('  h=%2d: %.4f A (%.2f%% of fundamental)\n', h, harmonics_sorted(i), 100*harmonics_sorted(i)/I1_rms);
end

fprintf('\n');
fprintf('========================================\n');
fprintf('      POWER QUALITY SUMMARY\n');
fprintf('========================================\n');
fprintf('Analysis Window: %.4f s (%d cycles)\n', duration_win, actual_cycles);
fprintf('Sampling Frequency: %.1f Hz\n', Fs);
fprintf('----------------------------------------\n');
fprintf('Real Power (P):         %8.2f kW\n', P*1e-3);
fprintf('Apparent Power (S):     %8.2f kVA\n', S*1e-3);
fprintf('Reactive Power (Q1):    %8.2f kvar\n', Q1*1e-3);
fprintf('----------------------------------------\n');
fprintf('True Power Factor (PF): %8.4f\n', PF_true);
fprintf('Displacement PF (DPF):  %8.4f\n', DPF);
fprintf('THDi (Phase A):         %8.2f %%\n', THDi*100);
fprintf('----------------------------------------\n');
fprintf('Phase A: V=%6.1f V, I=%6.2f A\n', Vrms_a, Irms_a);
fprintf('Phase B: V=%6.1f V, I=%6.2f A\n', Vrms_b, Irms_b);
fprintf('Phase C: V=%6.1f V, I=%6.2f A\n', Vrms_c, Irms_c);
fprintf('========================================\n');

create_plots = true;

if create_plots
    figure('Name', 'Power Quality Analysis', 'Position', [100 100 1200 800]);

    subplot(3,2,1);
    plot(t_win*1000, va_win, 'r', t_win*1000, vb_win, 'g', t_win*1000, vc_win, 'b', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Three-Phase Voltages');
    legend('Phase A', 'Phase B', 'Phase C');
    xlim([t_win(1) t_win(1)+0.06]*1000);

    subplot(3,2,2);
    plot(t_win*1000, ia_win, 'r', t_win*1000, ib_win, 'g', t_win*1000, ic_win, 'b', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (ms)');
    ylabel('Current (A)');
    title('Three-Phase Currents');
    legend('Phase A', 'Phase B', 'Phase C');
    xlim([t_win(1) t_win(1)+0.06]*1000);

    subplot(3,2,3);
    plot(t_win*1000, p_inst, 'k', t_win*1000, P*ones(size(t_win)), 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (ms)');
    ylabel('Power (W)');
    title(sprintf('Instantaneous Power (P_{avg} = %.2f kW)', P*1e-3));
    legend('p(t)', 'Average P');
    xlim([t_win(1) t_win(1)+0.06]*1000);

    subplot(3,2,4);
    bar(1:max_harmonic, Irms_harmonic, 'FaceColor', [0.2 0.6 0.8]);
    grid on;
    xlabel('Harmonic Order');
    ylabel('RMS Current (A)');
    title(sprintf('Current Harmonic Spectrum (THDi = %.2f%%)', THDi*100));
    xlim([0 max_harmonic+1]);

    subplot(3,2,5);
    f_plot = f_axis(f_axis <= f1*max_harmonic);
    Ia_plot = Ia_mag(f_axis <= f1*max_harmonic);
    plot(f_plot, Ia_plot, 'b', 'LineWidth', 1);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('RMS Current (A)');
    title('FFT Spectrum (Phase A Current)');
    xlim([0 f1*max_harmonic]);

    subplot(3,2,6);
    axis off;
    summary_text = {
        '\bf Power Quality Summary'
        ' '
        sprintf('P = %.2f kW', P*1e-3)
        sprintf('S = %.2f kVA', S*1e-3)
        sprintf('Q1 = %.2f kvar', Q1*1e-3)
        ' '
        sprintf('PF = %.4f', PF_true)
        sprintf('DPF = %.4f', DPF)
        sprintf('THDi = %.2f %%', THDi*100)
        ' '
        sprintf('Duration: %.3f s', duration_win)
        sprintf('Cycles: %d', actual_cycles)
    };
    text(0.1, 0.9, summary_text, 'FontSize', 11, 'VerticalAlignment', 'top');
end

fprintf('\n✓ Analysis complete!\n');
