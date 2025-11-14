% Electrolyzer System Parameters
VLL_MV = 23000;    % MV grid voltage (line-line RMS)
f      = 50;       % Grid frequency (Hz)

% Electrolyzer rated power and DC bus voltage
P_DC   = 5e6;      % 5 MW load
V_DC   = 600;      % 600 V DC bus

% Short-circuit strength & derived grid data
SCR    = 10;                 % Short Circuit Ratio
S_sc   = SCR * P_DC;         % Short-circuit power (VA)

% Grid equivalent impedance (Thevenin)
% Here we choose X/R ≈ 10 (common for MV networks)
R_th   = 1;                  % ohm
L_th   = 10/(2*pi*f);        % H

% Transformer secondary LV voltage
VLL_LV = 200;       % LV side rated line-line voltage (RMS)
S_tr   = 6e6;       % Transformer rating (VA)
Z_sc   = 0.07;      % Transformer short-circuit impedance (pu)
trr1=VLL_MV/sqrt(3)/VLL_LV
trr2=VLL_MV/VLL_LV