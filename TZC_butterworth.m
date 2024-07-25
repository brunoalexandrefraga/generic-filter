%------------- FILTRO ANALÓGICO -------------%
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

% Função de transferência normalizada
[num, den] = zp2tf([], s_bar, G0);

% Desnormalização da frequência
[num, den] = lp2lp(num, den, omega_p);

% Função de transferência desnormalizada
Ha = tf(num, den);







%------------- ANALÓGICO PARA DIGITAL -------------%
% Calcula Ga = Ha / s
den = [den, 0]; % den * s
num = [0, num]; % num

Ga = tf(num, den);

[r,p,k] = residue(num, den);


% Inicialize a função simbólica s
syms s t z n
g_t = 0;

% Transformada de Laplace inversa para cada termo de fração parcial
for i = 1:length(r)
    term = r(i) / (s - p(i));
    g_t = g_t + ilaplace(term, s, t);
end

%g_t = g_t * heaviside(t);

Delta_t = 1 / Fs;

t = n * Delta_t;
g_t = subs(g_t);

G_z = ztrans(g_t, n, z);

H_z = G_z * (z - 1) / z;


[H_z_num, H_z_den] = numden(H_z);
num_coeffs = sym2poly(H_z_num);
den_coeffs = sym2poly(H_z_den);











%------------- COEFICIENTES DIGITAIS -------------%
num_coeffs = real(num_coeffs);
den_coeffs = real(den_coeffs);

[ss,gn] = tf2sos(num_coeffs, den_coeffs);

ss = ss / 2 * 32678





















%------------- GRÁFICO DIGITAL -------------%
% Calculando a resposta em frequência
[Hz, Freq] = freqz(num_coeffs, den_coeffs, 'half', 4096);

plot(Freq, mag2db(abs(Hz)))
axis([0 pi -60 5])
grid
xlabel("Angular Frequency (rad/s)")
ylabel("Magnitude (dB)")

% Ajustando os ticks e labels do eixo x
xticks([0, pi/6, pi/3, pi/2, pi]);
xticklabels({'0', '\pi/6', '\pi/3', '\pi/2', '\pi'});