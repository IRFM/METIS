function ya0 = z0yield(gaz,wall,ein)

% parametres physiques
% The Plasma Boundary of magnetic fusion devices
% P.C Stangeby 
% chapitre 3
% No data for boron or He3, return an error if used
switch gaz
case 'H'
	switch wall
	case 'C'
		Eth  = 35;
		Etf  = 415;
		Q    = 0.035;
	case 'W'
		Eth  = 443;
		Etf  = 9870;
		Q    =  0.007;
	case 'Be'
		Eth  = 20;
		Etf  = 256;
		Q    =  0.1;
	otherwise
		error('Z0YIELD : unknow wall material');
	end
case 'D'
	switch wall
	case 'C'
		Eth  = 30;
		Etf  = 447;
		Q    =  0.1;
	case 'W'
		Eth  = 220;
		Etf  = 9923;
		Q    =  0.019;
	case 'Be'
		Eth  = 9;
		Etf  = 282;
		Q    = 0.3;
	otherwise
		error('Z0YIELD : unknow wall material');
	end
case 'T'
	switch wall
	case 'C'
		Eth  = 30;
		Etf  = 479;
		Q    =  0.2;
	case 'W'
		Eth  = 140;
		Etf  = 9977;
		Q    =  0.038;
	case 'Be'
		Eth  = 21;
		Etf  = 308;
		Q    = 0.24;
	otherwise
		error('Z0YIELD : unknow wall material');
	end
case 'He'
	switch wall
	case 'C'
		Eth  = 29;
		Etf  = 1087;
		Q    = 0.32;
	case 'W'
		Eth  = 110;
		Etf  = 20373;
		Q    = 0.106;
	case 'Be'
		Eth  = 30;
		Etf  = 780;
		Q    = 0.59;
	otherwise
		error('Z0YIELD : unknow wall material');
	end
otherwise
	error('Z0YIELD : unknow gaz');
end

e0    = logspace(-1,7,101);
epsi  = e0 ./ Etf;
sn    = 3.441 .* sqrt(epsi) .* log(epsi+exp(1)) ./ (1 + 6.355 .* sqrt(epsi) + epsi .* ( 6.882 .* sqrt(epsi) - 1.708));
d     = Eth ./ max(eps,e0);
g     = (1 - d .^ (2/3)) .* (1 - d .^ 2);
y     = Q .* sn .* g .* (e0>=Eth); 


% convolution par une gaussienne thermique
ein   = ein(:);
vt    = ones(size(ein));
ve    = ones(size(e0));
pe    = exp(- (vt*e0) ./ max(eps,ein * ve));

ya0   = trapz(e0,pe .* (vt * y),2) ./ max(eps,trapz(e0,pe,2));  


