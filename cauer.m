wp = 2; % [kHz]
ws = 3.2; % [kHz]

Amax = 0.8; % [dB]
Amin = 50; % [dB]

n_wp = wp / wp; % wp normalized
n_ws = ws / wp; % ws normalized

k = 1 / n_ws;
k_ = (1 - k ^ 2)^(1 / 2);
q0 = (1 / 2) * ((1 - k_^(1 / 2)) / (1 + k_^(1 / 2)));
q = q0 + 2 * q0 ^ 5 + 15 * q0 ^ 9 + 150 * q0 ^ 13;
d = (10 ^ (0.1 * Amin) - 1) / (10 ^ (0.1 * Amax) - 1);
n = log10(16 * d) / log10(1 / q);
n_ceil = ceil(n);


L = (1 / 2 * n) * log((10 ^ (0.05 * Amax) + 1) / (10 ^ (0.05 * Amax) - 1));


o0_N = 0;
M = 5;
for m = 0 : M
    o0_N = o0_N + (-1) ^ m * q ^ (m * (m + 1)) * sinh((2 * m + 1) * L);
end

o0_N = o0_N * 2 * q ^ (1 / 4);




o0_D = 0;
M = 5;
for m = 0 : M
    o0_D = o0_D + (-1) ^ m * q ^ (m ^ 2) * cosh(2 * m * L);
end

o0_D = 1 + 2 * o0_D;


o0 = o0_N / o0_D;


s = tf('s');

D0 = 1;
r = 0;
if (-1) ^ n == -1
    r = (n - 1) / 2;
    D0 = w_s ^ (-1/2) * s + o0;
else
    r = n / 2;
end



disp(n_ceil);