% calcul  de l'increment de la section effcace d'arret des neutre rapides du au ions rapides
% pour la phase a basse densite
function sv = z0nbistopfast(A,E,pa,ne,te,pnbi,tause)

% debit 
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
S = pnbi ./ (E .* phys.e);

% ref : K. Okano, 28 EPS 2001
% changement d'unite
E   = E ./ A ./ 1e6;
te  = te ./ 1e3./ 10;
eb  = E ./ te;
pa  = min(1,abs(pa));
% coefficient
Eeb = 0.57658 + 1.6432 .* eb -4.8886 .* eb .^ 2 + 5.3932 .* eb .^ 3 - 2.9226 .* eb .^ 4 + 0.77156 .* eb .^ 5 -0.079095 .* eb .^ 6;
Eeb(eb > pi) = 0;
gp  = 0.30812 .* pa + 11.992 .* pa .^ 2 -68.971 .* pa .^ 3 +165.94 .* pa .^4 -177.33 .* pa .^ 5 + 69.833 .* pa .^ 6;
% force 
fbp = te .^ -0.9027 .* Eeb .* gp .* 1e-19;
% section
sv  = S .* tause ./ ne .* fbp;
