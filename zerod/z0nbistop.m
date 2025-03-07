% calcul la section effcace d'arret des neutre rapides
function sv = z0nbistop(A,E,ne,te,nhe,nlow,nhigh,zlow,zhigh)

% ref : R.K. Janev NF  29 , 1989 p 2125-
% changement d'unite
E   = log(E ./ A ./ 1e3);
n0 = 1e13;
nel = ne;
ne  = log(ne ./ 1e6 ./ n0);
te  = log(te ./ 1e3);
ne2 = ne .^ 2; 
E2  = E .^ 2; 

% coefficient
A111 = 4.4;
A112 = -2.49e-2;
A121 = 7.46e-2;
A122 = 2.27e-3;
A131 = 3.16e-3;
A132 = -2.78e-5;
A211 = 2.3e-1;
A212 = -1.15e-2;
A221 = -2.55e-3;
A222 = -6.20e-4;
A231 = 1.32e-3;
A232 = 3.38e-5;
% S1
s1  = A111      + A121 .* ne      + A112 .* te      + A122 .* ne .* te + ...
      A131 .* ne2    +  A132 .* ne2 .* te + A231 .* E .* ne2 + A232 .* E .* ne2 .* te + ...
      A211 .* E + A212 .* E .* te + A221 .* E .* ne + A222 .* E  .* te .* ne;
      
% coeff He
B111 = -2.36;
B112 = 1.85e-1;
B121 = -2.5e-1;
B122 = -3.81e-2;
B211 = 8.49e-1;
B212 = -4.78e-2;
B221 = 6.77e-2;
B222 = 1.05e-2;
B311 = -5.88e-2;
B312 = 4.34e-3;
B321 = -4.48e-3;
B322 = -6.76e-4;   
% She
she = B111  + B112 .* te + B121 .* ne + B122 .* ne .* te + B211 .* E + B212 .* E .* te + B221 .* E .* ne + ...
      B222 .* E .* ne .* te + B311 .* E2 + B312 .* E2 .* te + B321 .* E2 .* ne + B322 .* E2 .* ne .* te;

% Coef C
B111 =-1.49;
B112 =-1.54e-2;
B121 =-1.19e-1; 
B122 =-1.5e-2;
B211 =5.18e-1;
B212 =7.18e-3;
B221 =2.92e-2;
B222 =3.66e-3;
B311 =-3.36e-2;
B312 =3.41e-4;
B321 =-1.79e-3;
B322 =-2.04e-4;   
% Sc
sc  = B111  + B112 .* te + B121 .* ne + B122 .* ne .* te + B211 .* E + B212 .* E .* te + B221 .* E .* ne + ...
      B222 .* E .* ne .* te + B311 .* E2 + B312 .* E2 .* te + B321 .* E2 .* ne + B322 .* E2 .* ne .* te;
      
% coef O  
B111 =-1.41e0;
B112 =-4.08e-4;
B121 =-1.08e-1; 
B122 =-1.38e-2;
B211 =4.77e-1;
B212 =1.57e-3;
B221 =2.59e-2;
B222 =3.33e-3;
B311 =-3.05e-2;
B312 =7.35e-4;
B321 =-1.57e-3;
B322 =-1.86e-4;   

% SO    
so  = B111  + B112 .* te + B121 .* ne + B122 .* ne .* te + B211 .* E + B212 .* E .* te + B221 .* E .* ne + ...
      B222 .* E .* ne .* te + B311 .* E2 + B312 .* E2 .* te + B321 .* E2 .* ne + B322 .* E2 .* ne .* te;

% Coef Fe
B111 =-1.03e0;
B112 =1.06e-1;
B121 =-5.58e-2; 
B122 =-3.72e-3;
B211 =3.22e-1;
B212 =-3.75e-2;
B221 =1.24e-2;
B222 =8.61e-4;
B311 =-1.87e-2;
B312 =3.53e-3;
B321 =-7.43e-4;
B322 =-5.12e-5;   
% S Fe
sfe = B111  + B112 .* te + B121 .* ne + B122 .* ne .* te + B211 .* E + B212 .* E .* te + B221 .* E .* ne + ...
      B222 .* E .* ne .* te + B311 .* E2 + B312 .* E2 .* te + B321 .* E2 .* ne + B322 .* E2 .* ne .* te;
      
% calcul final
if zlow < 6
	vx   = max(0,(zlow - 2) ./ 4);
	slow = sc .* vx + (1-vx) .* she;
elseif (zlow >= 6) && (zlow <= 8)
	vx   = max(0,(zlow - 6) ./ 2);
	slow = so .* vx + (1-vx) .* sc;
else
	vx   = min(1,max(0,(zlow - 8) ./ 18));
	slow = sfe .* vx + (1-vx) .* so;

end      
if zhigh < 6
	vx   = max(0,(zhigh - 2) ./ 4);
	shigh = sc .* vx + (1-vx) .* she;
elseif (zhigh >= 6) && (zhigh <= 8)
	vx   = max(0,(zhigh - 6) ./ 2);
	shigh = so .* vx + (1-vx) .* sc;
else
	vx   = min(1,max(0,(zhigh - 8) ./ 18));
	shigh = sfe .* vx + (1-vx) .* so;
end      
% 
sv = 1e-4 .* 1e-16 .* exp(s1) ./ exp(E) .* (1 + (nhe .* 2 .* she + nlow .* zlow .* (zlow -1) .* slow +  ...
                                            nhigh .* zhigh .* (zhigh -1) .* shigh) ./ nel);
					    					    
      
