f_p = 2; % [kHz]
f_s = 3.2; % [kHz]

omega_p = 2 * pi * f_p; % [kHz]
omega_s = 2 * pi * f_s; % [kHz]

Amax = 0.8; % [dB]
Amin = 50; % [dB]

n_wp = omega_p / omega_p; % wp normalized
n_ws = omega_s / omega_p; % ws normalized

epsilon = (10 ^ (Amax / 10) - 1) ^ (1 / 2);

a = 10 ^ (0.1 * Amin) - 1;
b = 10 ^ (0.1 * Amax) - 1;
c = log10(a / b);
d = 2 * log10(n_ws);

n = ceil(c / d);

sigma_bar = zeros(1, n);
omega_bar = zeros(1, n);
s_bar = zeros(1, n);

for k = 1 : n
   factor = 1;

   sigma_bar(k) = epsilon ^ (-1 / n) * sin(((2 * k - 1) * pi) / (2 * n));
   omega_bar(k) = epsilon ^ (-1 / n) * cos(((2 * k - 1) * pi) / (2 * n));

   if sigma_bar(k) > 0
       factor = -1;
   end

   s_bar(k) = factor * (sigma_bar(k) + 1i * omega_bar(k));
end

G0 = epsilon ^ (-1);

T_bar = tf(G0, s_bar);






% Função de transferência normalizada
[num, den] = zp2tf([], s_bar, G0);

% Desnormalização da frequência
[num, den] = lp2lp(num, den, omega_p);

% Resposta em Frequência
[H, w] = freqs(num, den, 4096);

plot(w/(2*pi), ...
    mag2db(abs(H)))
axis([0 5 -60 5])
grid
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")