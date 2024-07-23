M = 5;

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
Ha = tf(num, den);



r_den = roots(den);

% Identifique raízes complexas conjugadas
r_complex = r_den(imag(r_den) ~= 0);

% Inicialize a lista de polinômios
den_biquad_polynomials = [];

% Forme os polinômios de segunda ordem
while ~isempty(r_complex)
    % Pegue o par de raízes complexas conjugadas
    root1 = r_complex(1);
    root2 = conj(root1);
    
    % Calcule os coeficientes do polinômio de segunda ordem
    a = 1;
    b = -(root1 + root2);
    c = root1 * root2;
    
    % Crie o polinômio de segunda ordem
    poly = [a -b c];
    
    % Adicione o polinômio à lista
    den_biquad_polynomials = [den_biquad_polynomials; poly];
    
    % Remova as raízes já processadas
    r_complex(1) = [];
end

% Identifique e adicione os polinômios com raízes reais
r_real = r_den(imag(r_den) == 0);

% Forme os polinômios de segunda ordem a partir das raízes reais
while length(r_real) >= 2
    % Pegue o par de raízes reais
    root1 = r_real(1);
    root2 = r_real(2);
    
    % Calcule os coeficientes do polinômio de segunda ordem
    a = 1;
    b = -(root1 + root2);
    c = root1 * root2;
    
    % Crie o polinômio de segunda ordem
    poly = [a -b c];
    
    % Adicione o polinômio à lista
    den_biquad_polynomials = [den_biquad_polynomials; poly];
    
    % Remova as raízes já processadas
    r_real(1:2) = [];
end

% Se houver uma raiz real isolada, adicione o polinômio correspondente
if ~isempty(r_real)
    root1 = r_real(1);
    poly = [0 1 -root1];
    den_biquad_polynomials = [den_biquad_polynomials; poly];
end

den_biquad_polynomials = unique(den_biquad_polynomials, 'rows', 'stable');






r_num = roots(num);

% Identifique raízes complexas conjugadas
r_complex = r_num(imag(r_num) ~= 0);

% Inicialize a lista de polinômios
num_biquad_polynomials = [];

% Forme os polinômios de segunda ordem
while ~isempty(r_complex)
    % Pegue o par de raízes complexas conjugadas
    root1 = r_complex(1);
    root2 = conj(root1);
    
    % Calcule os coeficientes do polinômio de segunda ordem
    a = 1;
    b = -(root1 + root2);
    c = root1 * root2;
    
    % Crie o polinômio de segunda ordem
    poly = [a -b c];
    
    % Adicione o polinômio à lista
    num_biquad_polynomials = [num_biquad_polynomials; poly];
    
    % Remova as raízes já processadas
    r_complex(1) = [];
end

% Identifique e adicione os polinômios com raízes reais
r_real = r_num(imag(r_num) == 0);

% Forme os polinômios de segunda ordem a partir das raízes reais
while length(r_real) >= 2
    % Pegue o par de raízes reais
    root1 = r_real(1);
    root2 = r_real(2);
    
    % Calcule os coeficientes do polinômio de segunda ordem
    a = 1;
    b = -(root1 + root2);
    c = root1 * root2;
    
    % Crie o polinômio de segunda ordem
    poly = [a -b c];
    
    % Adicione o polinômio à lista
    num_biquad_polynomials = [num_biquad_polynomials; poly];
    
    % Remova as raízes já processadas
    r_real(1:2) = [];
end

% Se houver uma raiz real isolada, adicione o polinômio correspondente
if ~isempty(r_real)
    root1 = r_real(1);
    poly = [0 1 -root1];
    num_biquad_polynomials = [num_biquad_polynomials; poly];
end

num_biquad_polynomials = unique(num_biquad_polynomials, 'rows', 'stable');



























Ga = Ha/s_bar;