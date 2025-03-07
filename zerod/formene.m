% calcul des profils de densite en mode H
function [nout,fbest] = formene(x,vpr,nem,ane,nebord,pep)

% a calculer une seule fois
persistent shapene3
persistent shapene1
x=linspace(0,1,21);
if isempty(shapene3) || any(~isfinite(shapene3(:))) || any(~isfinite(shapene1(:)))
    nbp3 = length(3:-0.01:0);
    shapene3 = NaN * ones(nbp3,length(x));
    ind = 1;
    for k=3:-0.01:0
        % calcul des profil normalises
        xx = [1 0.95 0 -0.05];
        yy = [0 k 1 1];
        shapene3(ind,:) = pchip(xx,yy,x);
        ind = ind +1 ;
    end
    nbp1 = length(1:-0.01:0);
    shapene1 = NaN * ones(nbp1,length(x));
    ind = 1;
    for k=1:-0.01:0
        % calcul des profil normalises
        xx = [1 0.95 0 -0.05];
        yy = [0 k 1 1];
        shapene1(ind,:) = pchip(xx,yy,x);
        ind = ind +1 ;
    end
end

% vecteur utiles
vt  = ones(size(nem,1),1);
ve  = ones(1,length(x));
vol = trapz(x,vpr,2);

% ne0 = (1 + ane) .* nem;
% calcul de la valeur attendue au centre
ne0  = nem .* (1 + ane);
% echelle des profiles
nes  = ne0 - nebord;
% rapport cible
rapc = ne0 ./ nem;
nbc  = nebord ./ nes;

% meilleur valeur du point a 0.95
fbest = NaN .* vt;
dbest = NaN .* vt;
nout  = NaN .* (vt * ve);
first = 1;

% profil creux
if any(ane < 0)
  kmax = 3;
  nbpk = size(shapene3,1);
else
  kmax = 1;
  nbpk = size(shapene1,1);
end

% boucle sur les piquage
% volume element
dvol = diff(cumtrapz(x,vpr,2),1,2);
%dxm = diff(x);
%dvol = (vpr(:,1:end-1) + vpr(:,2:end)) ./ 2 .* dxm(ones(size(vt)),:);

%for k=kmax:-0.01:0
for k=1:nbpk
	% calcul des profil normalises
	%xx = [1 0.95 0 -0.05];
	%yy = [0 k 1 1];
	%x=linspace(0,1,21);
	%y=pchip(xx,yy,x);
	switch kmax
    case 1
                y = shapene1(k,:);
    otherwise
                y = shapene3(k,:);
    end
	% calcul de la valeur moyenne 
	%ym =  trapz(x,(vt * y) .* vpr,2) ./ vol + nbc;
	%ym =  trapz(x,y(ones(size(vt)),:) .* vpr,2) ./ vol + nbc;
  dym    = (y(1:end-1) + y(2:end)) ./ 2;
  ym     = sum(dym(ones(size(vt)),:) .* dvol,2)./ vol + nbc;
	rap = (y(1) + nbc) ./ max(eps,ym);
	
	% premier passage
	if first == 1
		first = 0;
		fbest(:) = k;
		dbest = abs(rap - rapc);
		%nout  = nes * y + nebord * ve;
		nout  = nes * y + nebord(:,ones(size(ve)));
	elseif nargin > 5 
		% securite sur le gradient de Te
		%nep  = nes * y + nebord * ve;
		nep  = nes * y + nebord(:,ones(size(ve)));
		teout = pep ./ max(1e13,nep);
                dte   = max(diff(teout,1,2),[],2);
		dtest = abs(rap - rapc);
		indok = find((dtest < dbest) & (dte < 0));
		fbest(indok) = k;
		dbest(indok) = dtest(indok);
		if ~isempty(indok)
			%nout(indok,:)  = nes(indok) * y  + nebord(indok) * ve;
			nout(indok,:)  = nes(indok) * y  + nebord(indok,ones(size(ve)));
		end

        else
		dtest = abs(rap - rapc);
		indok = find(dtest < dbest);
		fbest(indok) = k;
		dbest(indok) = dtest(indok);
		if ~isempty(indok)
			%nout(indok,:)  = nes(indok) * y  + nebord(indok) * ve;
			nout(indok,:)  = nes(indok) * y  + nebord(indok,ones(size(ve)));
		end
	end
end

% cas pathologique
indbad = find(nem < nebord);
if ~isempty(indbad)
	nout(indbad,:) = nebord(indbad) * ve;
	fbest(indbad) = NaN;
end

%  figure(21)
%  clf
%  plot(fbest)
%  drawnow

%  figure(151);clf
%  subplot(2,2,1)
%  plot((nout(:,end) - nebord)./nebord);
%  subplot(2,2,2)
%  plot((nout(:,1) - ne0) ./ ne0);
%  subplot(2,2,3)
%  plot((trapz(x,nout .* vpr,2) ./trapz(x,vpr,2) - nem) ./ nem);
%  subplot(2,2,4)
%  plot((nout(:,1) ./(trapz(x,nout .* vpr,2) ./trapz(x,vpr,2)) - ane - 1) ./ (ane + 1));
%  drawnow



