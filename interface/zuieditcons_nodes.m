function nodes = zuieditcons_nodes(x,y,code_retour,spline_flag)

% si length(x) < 2
	modessc =0;
if length(x) <2
	x=[0,1];
	y=[0,0];
elseif length(x) ==2
	modessc = 1;
end

% orientation
if size(x,1) > 1
	spline_flag =0;
else
	spline_flag =1;
end
x=x(:);
y=y(:);

% selon le code retour
ymod   = abs(y);
yangle = angle(y);
yreal  = real(y);
yimag  = imag(y);
[ypuiss,ytoro,ypolo] = zdecodefce(y);

switch code_retour
case 'angle'
	y     = yangle;
	multi = 1;
case 'degres'
	multi = 1;
	y     = yangle / pi * 180;
case 'polo'
	multi = 1;
	y     = ypolo;
case 'toro'
	multi = 1;
	y     = ytoro;
case {'abs','pfce'}
	y  = ymod;
	yy = ymod;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end

case 'idn'
	y  = yreal;
	yy = yreal;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end

case 'didndt'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end

case 'idn2'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end

case 'nbar'
	y  = yreal;
	yy = yreal;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end

case 'gaspuff'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end

case 'fisrt'
	y  = yreal;
	yy = yreal;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
case 'second'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
case 'nhnd'
	y  = yreal;
	yy = yreal;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
case 'ntnd'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
otherwise
	%y  = y;
	yy = ymod;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
end

if modessc == 1
	multi =1 ;
end


% extraction des points
x_mem = x;
y_mem = y;
curve_init =0;
if all(~isfinite(y))
    x = linspace(min(x),max(x),11);
    y = zeros(size(x));
    multi = 1;
elseif length(x) <31
    % on garde tout les points
    % mode lineaire par defaut
    x = x';
    y = y';
else

   % choix des noeux
   [x,y]= zextraitnoeud(x,y,11 + 20 .* (1 - (spline_flag~=0)),spline_flag);

   curve_init = 1;
end

% retour
nodes = x;