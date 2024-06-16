wp = 2; % [kHz]
ws = 3.2; % [kHz]

Amax = 0.8; % [dB]
Amin = 50; % [dB]

n_wp = wp / wp; % wp normalized
n_ws = ws / wp; % ws normalized

E = (10 ^ (Amax / 10) - 1) ^ (1 / 2);

a = 10 ^ (0.1 * Amin) - 1;
b = 10 ^ (0.1 * Amax) - 1;
c = log10(a / b);
d = 2 * log10(n_ws);

n = ceil(c / d);

n_o = zeros(n, 1);
n_w = zeros(n, 1);
n_s = zeros(n, 1);

for k = 1 : n
   factor = 1;

   n_o(k) = E ^ (-1 / n) * sin(((2 * k - 1) * pi) / (2 * n));
   n_w(k) = E ^ (-1 / n) * cos(((2 * k - 1) * pi) / (2 * n));

   if n_o(k) > 0
       factor = -1;
   end

   n_s(k) = factor * (n_o(k) + 1i * n_w(k));
end

plot(n_s, 'x');