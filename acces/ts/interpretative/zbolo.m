% fit gaussien du profil de puissance rayonnee
% syntaxe : [centre,largeur,cprad] = bolo_width(choc,plotonoff);
%           x = linspace(0,1);   % coordonnees r/a 
%           fit = exp(-(x-centre) .^ 2 ./ largeur .^2);
% cprad = certification : 0 = ok
% choc  = numero du choc
% si plotonoff = 1 -> trace le graphe du resultat
function [centre,largeur,cprad] = zbolo(choc,plotonoff)

% test des entrees
if nargin < 2
	plotonoff = 0;
end

% valeur par defaut
centre  = NaN;
largeur = NaN;

% lecture des puissances 
[prad,tprad,iprad,cprad] = tsbase(choc,'gpbolo1');
if isempty(prad)
	return
end
[cert,vers,date,heure,unix,uniy,uniz,nomdon,type]=tsbase_cert(cprad);
if any(cert(1:3) == -2)
	return
end

% svd de filtrage
[u,s,v] = svd(prad,0);
% on garde le profil principal
forme_h = abs(v(1:8,1)');
forme_l = abs(v(end:-1:9,1)');

% lecture des angles
theta   = tsmat(choc,'DBOLO;CONFIG_CAM1;THETA');
if isempty(theta)
	cprad(2:(1+cprad(1))) = -2;
	return
end
pimp    = tsmat(choc,'DBOLO;CONFIG_CAM1;IMPACT');
ssol    = tsmat(choc,'DBOLO;CONFIG_CAM1;SOLIDE');

lz      = pimp ./ cos( pi/2 - abs(theta));
lz      = max(lz) -lz;
zsa     = lz ./ max(lz);
dforme_h  = gradient(forme_h);
if any(dforme_h > 0)
		indm = min(find((dforme_h>0) & (zsa < max(zsa))));
		forme_h(indm:end) = forme_h(indm);
end 
forme_h   = forme_h - min(forme_h);
forme_h   = forme_h .* abs(cos(theta));
forme_h   = forme_h ./ max(forme_h);
 
dforme_l  = gradient(forme_l);
if any(dforme_l > 0)
		indm = min(find((dforme_l>0) & (zsa < max(zsa))));
        if indm <= 2
            indm = 3;
        end
		forme_l(indm:end) = forme_l(indm);
end 
forme_l   = forme_l - min(forme_l);
forme_l   = forme_l .* abs(cos(theta));
forme_l   = forme_l ./ max(forme_l);
 

x       = linspace(0,1);
ff      = pchip(zsa,forme_h,x); 

% centre 
indc    = min(find(ff == max(ff)));
centre  = x(indc);
xc      = x(1:indc);
ffc     = ff(1:indc);
d       = abs(ffc - max(ffc)./exp(1));
lc      = xc(min(find(d == min(d))));
largeur = centre -lc;

if largeur > 0.5
   	disp('largeur non valide')
 	ff      = pchip(zsa,forme_l,x); 

	% centre 
	indc    = min(find(ff == max(ff)));
	centre  = x(indc);
	xc      = x(1:indc);
	ffc     = ff(1:indc);
	d       = abs(ffc - max(ffc)./exp(1));
	lc      = xc(min(find(d == min(d))));
	largeur = centre -lc;
	forme   = forme_l;
	if largeur > 0.5
  		disp('largeur non valide')
  		largeur = NaN;
		cprad(2:(1+cprad(1))) = -2;
	end	
	cprad(2:(1+cprad(1)))  = -1;

else
	forme   = forme_h;
end


if plotonoff
	h=figure;
	set(h,'defaultaxesfontsize',12,'defaultaxesfontweight', ...
	'bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
	fit = exp(-(x-centre) .^ 2 ./ largeur .^2);
	fit = max(ff) .* fit ./ max(fit);  
	plot(zsa,forme,'or',x,ff,'.b',x,fit,'k');
	title(int2str(choc))
	xlabel('r/a');
	ylabel('Prad (su)');
	legend('measurement','interpolation','fit','Location','NorthWest');
end
