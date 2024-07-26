%------------- FILTRO ANALÓGICO -------------%
M = 5;

Fs = 48; % [kHz]

f_p = 2; % [kHz]
omega_p = 2*pi*f_p; % [rad/s]


f_s = 3.2; % [kHz]
omega_s = 2*pi*f_s; % [rad/s]


Amax = 0.8; % [dB]
Amin = 50; % [dB]

omega_p_bar = omega_p / omega_p;
omega_s_bar= omega_s / omega_p;

k = 1 / omega_s_bar;
k_line = (1 - k ^ 2)^(1 / 2);
q0 = (1 / 2) * ((1 - k_line^(1 / 2)) / (1 + k_line^(1 / 2)));
q = q0 + 2 * q0 ^ 5 + 15 * q0 ^ 9 + 150 * q0 ^ 13;
d = (10 ^ (0.1 * Amin) - 1) / (10 ^ (0.1 * Amax) - 1);
n = ceil(log10(16 * d) / log10(1 / q));

n_is_odd = mod(n, 2) == 1;

Lambda = (1 / (2 * n)) * log((10 ^ (0.05 * Amax) + 1) / (10 ^ (0.05 * Amax) - 1));

sigma0_bar_num = 0;
for m = 0 : M
    sigma0_bar_num = sigma0_bar_num + (-1) ^ m * q ^ (m * (m + 1)) * sinh((2 * m + 1) * Lambda);
end
sigma0_bar_num = sigma0_bar_num * 2 * q ^ (1 / 4);

sigma0_bar_den = 0;
for m = 1 : M
    sigma0_bar_den = sigma0_bar_den + (-1) ^ m * q ^ (m ^ 2) * cosh(2 * m * Lambda);
end
sigma0_bar_den = 1 + 2 * sigma0_bar_den;


sigma0_bar = sigma0_bar_num / sigma0_bar_den;


s_bar = tf('s');

D0 = 1;
r_den = 0;
if n_is_odd
    r_den = (n - 1) / 2;
    D0 = omega_s_bar ^ (-1/2) * s_bar + sigma0_bar;
else
    r_den = n / 2;
end

W = ((1 + k * sigma0_bar ^ 2) * (1 + (sigma0_bar ^ 2 / k))) ^ (1 / 2);


V = zeros(1, r_den);
Omega = zeros(1, r_den);
A0 = zeros(1, r_den);
B0 = zeros(1, r_den);
B1 = zeros(1, r_den);
G0 = 1;
T = 1;

for i = 1 : r_den
    mu = i;

    if ~n_is_odd
        mu = i - 1 / 2;
    end

    Omega_num = 0;
    for m = 0 : M
        Omega_num = Omega_num + (-1) ^ m * q ^ (m * (m + 1)) * sin(((2 * m + 1) * pi * mu) / n);
    end
    Omega_num = 2 * q ^ (1 / 4) * Omega_num;

    Omega_den = 0;
    for m = 1 : M
        Omega_den = Omega_den + (-1) ^ m * q ^ (m ^ 2) * cos((2 * m * pi * mu) / n);
    end
    Omega_den = 1 + 2 * Omega_den;

    Omega(i) = Omega_num / Omega_den;

    V(i) = ((1 - k * Omega(i) ^ 2) * (1 - Omega(i) ^ 2 / k)) ^ (1 / 2);

    A0(i) = 1 / Omega(i) ^ 2;

    B0(i) = ((sigma0_bar * V(i)) ^ 2 + (Omega(i) * W) ^ 2) / (1 + sigma0_bar ^ 2 * Omega(i) ^ 2) ^ 2;

    B1(i) = 2 * sigma0_bar * V(i) / (1 + sigma0_bar ^ 2 * Omega(i) ^ 2);

    G0 = G0 * B0(i) / A0(i) * omega_s_bar ^ (1 / 2);

    T = T * ((s_bar ^ 2) + omega_s_bar * A0(i)) / (s_bar ^ 2 + omega_s_bar ^ (1 / 2) * B1(i) * s_bar + omega_s_bar * B0(i));
end

if n_is_odd
    G0 = sigma0_bar * G0;
else
    G0 = 10 ^ (-0.05 * Amax) * G0;
end

T_bar = G0 / D0 * T;

% Resposta em Frequência
[num, den] = tfdata(T_bar, 'v');

% Desnormalização da frequência
[num, den] = lp2lp(num, den, omega_p);

% Função de transferência desnormalizada
H_s = tf(num, den);

























%------------- ANALÓGICO PARA DIGITAL -------------%
Ts = 1/Fs; % Período de amostragem

% Criar a função de transferência no domínio s
sys_s = tf(num, den);

% Obter polos e zeros no domínio s
[zeros_s, poles_s, gain_s] = tf2zp(num, den);

% Inicializar arrays para os novos numeradores e denominadores no domínio z
num_z = 1;
den_z = 1;

% Calcular os valores para cada par de polos e zeros complexos conjugados
for k = 1:2:length(poles_s)
    a = real(poles_s(k));
    b = imag(poles_s(k));
    if b ~= 0 % Considerar apenas polos/zeros complexos conjugados
        % Polinômio no domínio z para o par de polos/zeros
        poly_z = [1 -2*exp(-a*Ts)*cos(b*Ts) exp(-2*a*Ts)];
        
        % Atualizar denominador
        den_z = conv(den_z, poly_z);
    else
        poly_z = [1 exp(-a*Ts)];
        den_z = conv(den_z, poly_z);
    end
end

% Repetir o processo para zeros
for k = 1:2:length(zeros_s)
    a = real(zeros_s(k));
    b = imag(zeros_s(k));
    if b ~= 0 % Considerar apenas polos/zeros complexos conjugados
        % Polinômio no domínio z para o par de polos/zeros
        poly_z = [1 -2*exp(-a*Ts)*cos(b*Ts) exp(-2*a*Ts)];
        
        % Atualizar numerador
        num_z = conv(num_z, poly_z);
    else
        poly_z = [1 exp(-a*Ts)];
        num_z = conv(num_z, poly_z);
    end
end

% Ajustar o ganho para corresponder o ganho original
gain_z = gain_s * (abs(den_z(1)) / abs(num_z(1)));
num_z = num_z * gain_z;

% Inverter os coeficientes do numerador e do denominador para obter expoentes negativos
num_z = fliplr(num_z);
den_z = fliplr(den_z);

% Normalizando coeficientes
num_z = num_z / num_z(1);
den_z = den_z / den_z(1);

% Criar a função de transferência no domínio z com expoentes negativos
H_z = tf(num_z, den_z, Ts, 'Variable', 'z^-1');

















%------------- COEFICIENTES DIGITAIS -------------%
[ss,gn] = tf2sos(num_z, den_z);

%ss(1, 1) = 0;
%ss(1, 2) = 1;

ss = ss / 2 * 32678;

disp(ss);












%------------- GRÁFICO DIGITAL -------------%
% Calculando a resposta em frequência
[Hz, Freq] = freqz(num_z, den_z, 'half', 4096);

plot(Freq, mag2db(abs(Hz)))
axis([0 pi -60 5])
grid
xlabel("Angular Frequency (rad/s)")
ylabel("Magnitude (dB)")

% Ajustando os ticks e labels do eixo x
xticks([0, pi/6, pi/3, pi/2, pi]);
xticklabels({'0', '\pi/6', '\pi/3', '\pi/2', '\pi'});



% Plotar os polos e zeros
figure;
zplane(num_z, den_z);
title('Diagrama de Polos e Zeros');



% Plotar os polos e zeros
figure;
h1 = zplane(num_z, den_z);

num_z_truncated = num_z / 2 * 32678;
den_z_truncated = den_z / 2 * 32678;

hold on;
h2 = zplane(num_z_truncated, den_z_truncated);
set(h2(1), 'Marker', 's', 'MarkerEdgeColor', 'r'); % Quadrados vermelhos para zeros

title('Diagrama de Polos e Zeros');

hold off