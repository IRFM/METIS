% test correction Ti/te
tiote = logspace(-1,1,1001);
ne = 1e19;
ni = 1 .* ne;
Ai = 184;
Zi = 74;
% physical constantes
phys = cphys;

term = Zi .^ 2 .*  (ni ./ ne) .* (phys.me ./ phys.ua ./ Ai) .* (tiote);

rep = exp(-(3 .* sqrt(pi) ./ 4 .* term)  .^ (2/3));

figure;
subplot(2,2,1)
plot(tiote,term);
subplot(2,2,2)
plot(tiote,rep);

