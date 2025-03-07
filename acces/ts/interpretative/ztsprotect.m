% fonction de retrait des temps invalides
function [pout,tout,xout]  = ztsprotect(pin,tin,xin,leg)

% recherche des profils a 0
switch leg
case 'nefit'
   lim = 1e18;
case 'nifit'
   lim = 1e18;
case 'tefit'
   lim = 30;
otherwise
   lim = 0;
end

switch leg
case {'nefit','nifit','tefit'}
	ind = find(any(pin ~= 0,2) & all(pin >= -lim,2));
otherwise
	ind = find(any(pin ~= 0,2));
end
if (length(tin) - length(ind)) > (0.2 * size(pin,1))
   fprintf('Nombreux profils invalides :  correction impossible de %s !\n',leg);
	pout = pin;
	tout = tin;
	xout = xin;
else
	pout   = pin(ind,:);
	tout   = tin(ind);
	if size(xin,1) == size(pin,1)
	    xout  = xin(ind,:);
	else
		xout = xin;
	end
end

if length(tout) ~= length(tin)
    fprintf(' %d profils suprimer dans %s.\n',length(tin) - length(tout),leg);
end
