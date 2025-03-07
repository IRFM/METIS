% calcul nd/ne + flag choc D
function [ndne,main,ndnev,badti,corti,certn,certi] = zndne(numchoc,zmain,temps,xfn,rsa,ne,ti,rm,a,puiss)

badti = 0;
corti = ones(size(temps));
corok = 1;
oldtemps = temps;
%if zmain == 1 
%     [fluxn,tn,cert] = tsbase(fix(numchoc),'gfluntn%1');
%else
%     [fluxn,tn,cert] = tsbase(fix(numchoc),'gfluntn%2');
%end
[fluxn,tn,void,certn] = tsbase(fix(numchoc),'gfluntn');

if isempty(fluxn)
	main = zmain;
	if zmain == 1
		ndne =  0.7 ;
		ndnev = ndne;
		certi = [];
	else
		ndne =  0.05 ;
		ndnev = ndne;
		certi = [];
	end
	ndnev = ndne;
	certi = certn;
	return
end

tn      = tn(:,1);
indbad = find(diff(tn) <=0);
while(~isempty(indbad) && ~isempty(tn))
   fluxn(indbad,:) = [];
   tn(indbad,:)    = [];
   indbad = find(diff(tn) <=0);   
end    


% correction idn
[iht,tiht,certi] = tsbase(fix(numchoc),'siht');
if ~isempty(iht)
   [tiht,iht] = antir(tiht,iht);
	fsh   = 1 ./ 2 ./ mean(diff(tiht));
	[bf,af] = butter(7,3/fsh);
	iht = medfilt1(abs(filtfilt(bf,af,iht)),5);
	sihtp = (iht > 0.5);
	siht = interp1(tiht,sihtp,tn(:,1),'nearest');
	ind = find(~isfinite(siht));
	if ~isempty(ind)
		siht(ind) = 0;
	end
	ind   = find(~siht);
	fluxn = fluxn(ind,:);
	tnn   = tn(ind);
	fluxn = interp1(tnn,fluxn,tn,'nearest');
end

fluxn   = 10 .^ fluxn;

l  = max(3,round(mean(diff(temps))./mean(diff(tn))));
if rem(l,2) == 0
   l = l +1;
end
[y1,ymin1,ymax1,ymean1,ystd1] = meanfilt1(fluxn(:,1),l);
[y2,ymin2,ymax2,ymean2,ystd2] = meanfilt1(fluxn(:,2),l);
fluxmin  = min(ymin1,ymin2);
fluxmax  = max(ymax1,ymax2);
fluxmean = 0.5 .* (ymean1 + ymean2);

fluxmin   = interp1(tn,fluxmin,temps,'nearest');
fluxmax   = interp1(tn,fluxmax,temps,'nearest');
fluxmean  = interp1(tn,fluxmean,temps,'nearest');

xxfn    = ones(length(temps),1) * xfn;
ind = find(sum(isnan(rsa),2)>0);
if ~isempty(ind)
  ne(ind,:)=[];;
  ti(ind,:)=[];;
  rsa(ind,:)=[];;
  rm(ind,:)=[];;
  xxfn(ind,:)=[];;
  a(ind) = [];
  fluxmax(ind) = [];
  fluxmin(ind) = [];
  fluxmean(ind) = [];
  puiss(ind)=[];
  corti(ind) = [];
  oldtemps(ind) = [];
end

ne      = tsplinet(rsa,ne,xxfn);
rm      = tsplinet(rsa,rm,xxfn);
ti      = tsplinet(rsa,ti,xxfn);
fn      = fluxneutron(ti / 1e3,ne,xfn,rm,a);

ndnevmax   = sqrt(fluxmax ./ fn);
ndnevmin   = sqrt(fluxmin ./ fn);
ndnevmean  = sqrt(fluxmean ./ fn);
if zmain   == 1 
	maskmax  = (fluxmax > 1e8) & (ndnevmax < 0.95) & (ndnevmax > 0.3) & (puiss < 0.5) & (imag(ndnevmax) == 0);
	maskmin  = (fluxmin > 1e8) & (ndnevmin < 0.95) & (ndnevmin > 0.3) & (puiss < 0.5) & (imag(ndnevmin) == 0);
else
	maskmax  = (fluxmax > 1e8) & (ndnevmax < 0.95)& (puiss < 0.5) & (imag(ndnevmax) == 0);
	maskmin  = (fluxmin > 1e8) & (ndnevmin < 0.95)& (puiss < 0.5) & (imag(ndnevmin) == 0);
end

indmax      = find(maskmax & ~maskmin);
indmin      = find(~maskmax & maskmin);     
ind2        = find(maskmax & maskmin);
ind0        = find(~maskmax & ~maskmin);

ndnev         = ndnevmean;
ndnev(ind2)   = 0.5 .* (ndnevmax(ind2) + ndnevmin(ind2));
ndnev(indmax) = ndnevmax(indmax);
ndnev(indmin) = ndnevmin(indmin);
fluxn         = fluxmean;
fluxn(ind2)   = 0.5 .* (fluxmax(ind2) + fluxmin(ind2));
fluxn(indmax) = fluxmax(indmax);
fluxn(indmin) = fluxmin(indmin);




indok   = find((fluxn > 1e8) & (ndnev < 0.95) & (puiss < 0.5) & (imag(ndnev) == 0));
indnok  = find((fluxn <= 1e8) | (ndnev >= 0.95) | (puiss >= 0.5) | (imag(ndnev) ~= 0));
if isempty(indok)
	indok   = find((fluxn > 1e8) & (ndnev < 0.95) & (imag(ndnev) == 0));
   indnok  = find((fluxn <= 1e8) | (ndnev >= 0.95) | (imag(ndnev) ~= 0));
   corok   = 0;
	if length(indok) < (0.5 .* (length(temps)))
		badti = 1;
	end
		
elseif length(indok) < (0.5 .* (length(temps) - length(find(puiss >= 0.5))))
	badti = 1;
end

if ~isempty(indok)
	ndne    = mean(ndnev(indok));
	ndnev(indnok) = ndne;
	
	% securite 
	ndnev    = min(2.*ndne,max(0.5.*ndne,ndnev));
	
   if corok == 1;
		indp     = find(puiss > 0.5);
		[fnc,inte] = fluxneutron(ti / 1e3,ndne .* ne,xfn,rm,a);
		tiv      = linspace(0.1,10,1001);
		svv      = sigvddn(tiv);
		tinte    = trapz(inte .*ti /1e3,2);
		svin     = sigvddn(tinte);
		svok     = fluxn ./ fnc .* svin;
		tiok     = spline(svv,tiv,real(svok));
		cor      = min(2,max(0.1,tiok ./tinte));
		fsh      = 1 ./ 2 ./ mean(diff(temps));
        order    = fix(min(7,length(cor)/3));
		[bf,af]  = butter(order,min(0.9,5/fsh));
		cor      = filtfilt(bf,af,cor);
		%corti(indp) =  max(0.92 .* cor(indp) ./ mean(cor(indok)),1);
      % seule une partie des neutrons porvient de l'augmentation de Ti 
		corti(indp) =  max((1./sqrt(2)) .* cor(indp) ./ mean(cor(indok)),1); 
		%corti(indp) =  cor(indp) ./ mean(cor(indok)); 
   end
else
	badti = 1;
	if zmain == 1
		ndne =  0.7 ;
	else
		ndne =  0.05 ;
	end
end
if ndne > 0.5 
   main = 1;
else
   main = zmain;
end

if length(oldtemps) ~= length(temps)
    
  ndnev = interp1(oldtemps,ndnev,temps,'nearest','extrap');
  corti = interp1(oldtemps,corti,temps,'nearest','extrap');
end

function [t,s] = antir(t,s)

if isempty(t)
	return
end

dt = gradient(t);
t0 = t(1);
indok = find(dt > 0);
dtm   = mean(dt(indok));
indnok  = find(dt <=0);
dt(indnok) =dtm;
t = cumsum(dt)+ t0;

