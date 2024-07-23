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










r = roots(den);

% Identifique raízes complexas conjugadas
r_complex = r(imag(r) ~= 0);

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
r_real = r(imag(r) == 0);

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




r = roots(num);

% Identifique raízes complexas conjugadas
r_complex = r(imag(r) ~= 0);

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
r_real = r(imag(r) == 0);

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


























s = tf('s');

Ga = Ha / s;