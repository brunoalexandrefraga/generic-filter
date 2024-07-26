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