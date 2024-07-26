% Dados da função de transferência no domínio s
num = [1 2 26]; % Numerador
den = [1 2 2]; % Denominador
Ts = 1/2; % Período de amostragem

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
        den_z = den_z .* poly_z;
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
        num_z = num_z .* poly_z;
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
sys_z = tf(num_z, den_z, Ts, 'Variable', 'z^-1');

% Exibir a função de transferência no domínio z
disp('Função de transferência no domínio z:');
sys_z