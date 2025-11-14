% Assuming Vab_Y and Vab_D are timeseries from PLECS
v1 = out.Vab_lv_1.Data;
v2 = out.Vab_lv_2.Data;
t  = out.Vab_lv_1.time;

% Extract steady-state section
f0 = 50;               % fundamental
T  = 1/f0;
t1 = t(end) - 10*T;    % last 10 cycles
idx = (t >= t1);
v1 = v1(idx);
v2 = v2(idx);
t  = t(idx);

% FFT calculation
N = length(t);
dt = mean(diff(t));
fs = 1/dt;
Y1 = fft(v1);
Y2 = fft(v2);

% Frequency index for 50 Hz
k = round(f0*N/fs);

phiY = angle(Y1(k));
phiD = angle(Y2(k));

% Phase difference in degrees
phase_shift = rad2deg( angle(exp(1j*(phiD - phiY))) );

fprintf('Phase shift (Vab_D relative to Vab_Y) = %.2f degrees\n', phase_shift);
