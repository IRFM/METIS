function yy = tpchip_interp1(x,y,xx,method,extrap)

% gestiopn des arguments
if nargin < 4
	method = 'pchip';
end
if nargin < 5 
	extrap = NaN;
end

switch method
case 'pchip'
	% mise en forme (pas de matrice pour le moment)
	x = x(:);
	y = y(:);
	if x(1) > x(end)
		x = flipud(x);
		y = flipud(y);
	end
	
	% appel mexfile
	if ~ischar(extrap)
		yy    = extrap .* ones(size(xx));
		indin = find((xx <= max(x)) & (xx >= min(x)));
		[yy(indin),d,fail] = tpchip(x,y,xx(indin));
	else
		[yy,d,fail] = tpchip(x,y,xx);
	end
	
	% controle extrapolation
	if ~ischar(extrap)
		indout = find((xx > max(x)) | (xx < min(x)));
		if ~isempty(indout) 
			yy(indout) = extrap;
		end
	end
	if size(xx,1) ~= size(yy,1)
		yy =yy';
	end

%  	figure(123);clf
%   	plot(x,y,'.r',xx,yy,'ob',xx,interp1(x,y,xx,method,extrap),'xk');
%  	drawnow

otherwise
	% methode matlab
	yy = interp1(x,y,xx,method,extrap);	

end