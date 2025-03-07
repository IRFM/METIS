function [zmain,nhnd,nonc,data,cert]=zsmasse(nchoc,zin);


data = [];
nhnd = [];
nonc = [];
zmain = [];



%RECUPERATION DES DONNEES BRUTES
% 15 canaux maxi pour Gmasse
[sp,tsp,void,cert] = tsbase(nchoc,'GMASSP3Q6B');
if isempty(sp)
 		[sp,tsp,void,cert] = tsbase(nchoc,'GMASSP3Q6A');
end
if isempty(sp)
  disp('pas de donnees GMASSE pour ce choc')
   return
end

[ip,tip] = tsbase(nchoc,'sipmes');

warning off

%MISE A 1 BIT DES VALEURS NEGATIVES
condi = sp < 1;
sp(condi) = ones (size(sp(condi)));
 
% FullRange
fr=1E-5;

%CALCUL DES RESULTATS
condi = sp < 55;

sp(condi)  = (fr*4.88e-3/220) .* sp(condi);
sp(~condi) = (fr/1000)*10 .^[(4.88e-3/3.3).*sp(~condi)];


tsp  = mean(tsp,2);
indr = find(diff(tsp) <= 0);
while ~isempty(indr)
	indr = find(diff(tsp) <= 0);
	if ~isempty(indr)
		sp(indr,:) = [];
		tsp(indr,:) = [];
	else 
   	break
	end
end

temps   = tsp;
riso    = NaN .* temps;
ri      = NaN .* temps;
rimp    = NaN .* temps;
nbe     = NaN .* temps;


%RECUPERATION DES PARAMETRES SUR SERVEUR 
ss=tsmat(nchoc,'DSTD2;PILOTAGE;BALCHALPM');
%ss=tsmat(nchoc,'DSTD2;PILOTAGE;BALCHA');
nli=size(ss,1);
inm=3:4:nli;
ira=inm+1;
nmasse=str2num(ss(inm,4:5));
range=str2num(ss(ira,4));

%  REMISE EN FORME DES SIGNAUX DES MASSES ACQUISES 2 FOIS
masse = [2 3 4 6 18 22 28 40 44];
for im = 1:9
  ind = find (nmasse==masse(im));
  if isempty(ind)
    eval(['m' int2str(masse(im)) '= [];'])
  elseif size(ind,1) > 1
    [mx,imax] = max(range(ind));
    if imax==2, ind = flipud(ind); end
    eval(['m' int2str(masse(im)) '= fusion(sp(:,ind(1)),sp(:,ind(2)),range(ind(1)),range(ind(2)));']) 
  else
    eval(['m' int2str(masse(im)) '= sp(:,ind)*10^(-range(ind));'])
  end
end

if  isempty(m6)
   M6 = 1e-15 * ones(size(m4));
else
   M6 = m6;
   condi  = m6 < 5e-13;
   M6(condi) = 1e-15 * ones(size(m6(condi)));
end

% Lectures des coefficients d'etalonnage
% C5c16576, C6c16576, coef
racine = fileparts(which('zsmasse'));

[nomfic]=ficpar(nchoc,'etalon',racine);
eval(sprintf('load %s/parametres/%s',racine, nomfic));

rop = polyval(C6c16576,log10(M6));

%   Etalonnage du 23/08/93
% 	COEFFICIENT D'ETALONNAGE
%Masse=[2=H2 ,3=HD ,4=D2  ,4=He ,18=H2O ,22=Ne ,28=CO ,40=Ar ,44=CO2]
%coef = [25e3 ,35e3 ,40625 ,76e3 ,96e3   ,253e4 ,11e4  ,19e4  ,2e5 ];

%	PRESSION EN PASCAL
H2=medfilt1(m2*coef(1),5);
HD=medfilt1(m3*coef(2),5);
D2=medfilt1(m4*coef(3),5);
D3=medfilt1((1000 * M6 ./ rop)*coef(3),5);
He=medfilt1((m4-1000 * M6 ./ rop)*coef(4),5);
% Mise a zero des valeurs de He negatives
condi = He < 0;
He(condi) = zeros(size(He(condi)));
H2O =medfilt1(m18* coef(5),5);
Ne =medfilt1(m22* coef(6),5);
CO =medfilt1(m28* coef(7),5);
Ar =medfilt1(m40* coef(8),5);
CO2 =medfilt1(m44* coef(9),5);

%Correction pour éviter ri=[]
% R.Masset le 27/07/99
nH2=0;nD2=0;nHD=0;
if ~isempty(H2), nH2=2*H2; end
if ~isempty(HD), nHD=HD; end
if ~isempty(D2), nD2=2*D2; end
nH=nH2+nHD;
pres_d =0;
pres_he =0;
% Critere de decision
cd2d3 = corrcoef(D2,D3);
cd2he = corrcoef(D2,He);
if ~all(isfinite(cd2he(:)))
  if zin == 1
    disp('Choc Deuterium')
    pres_d = 1;
    He=[];
    nD=nD2+nHD;
    ri=nH./(nH+nD); % RAPPORT ISOTOPIQUE  
  else
    disp('choc Helium')
	 pres_he = 1;
    D2=[];
    He=m4*coef(4);
  end
  zmain = zin;
elseif cd2d3(2) > 0.95 &  cd2he(2) <= 0.9  % D2=D2;
   disp('Choc Deuterium')
   pres_d = 1;
   He=[];
   nD=nD2+nHD;
   ri=nH./(nH+nD); % RAPPORT ISOTOPIQUE
	zmain = 1;
elseif cd2he(2) > 0.95 &  cd2d3(2) <= 0.9  
   disp('choc Helium')
	pres_he = 1;
   D2=[];
   He=m4*coef(4);
	zmain = 2;
else                 % He=He
   disp('choc Helium + Deuterium')
   pres_d  = 1;
	pres_he = 1;
   D2=D3;
   if ~isempty(D2), nD2=2 * D2; end
   nD=nD2+nHD;
   ri=nH./(nH+nD); % RAPPORT ISOTOPIQUE
	zmain = 1 + i * 2;
end

if sum(isfinite(ri)) > 10
	riso  = ri;
	indok = find((ri > 0) & (ri < 1) & isfinite(ri));
	if length(indok) > 10
		indnok = find((ri <= 0) |(ri >= 1) | ~isfinite(ri));
		riso(indnok) = mean(ri(indok));
	else
		riso = NaN .* temps;
	end
end

rimp  = (CO + CO2) ./ (H2O + CO + 2.* CO2); 
indok = find((rimp >= 0.01) &(rimp < 1) & isfinite(rimp));
if length(indok) > 10
		indnok = find((rimp < 0.01) | (rimp >= 1) | ~isfinite(rimp));
		rimp(indnok) = mean(rimp(indok));
else
		rimp = NaN .* temps;
end

% estimation 'zeff'
nbe   = 0;
nbz2  = 0;
nbd   = 0;
nbh   =0;
if ~isempty(H2) 
	nbe  = nbe + H2 * 2;
	nbz2 = nbz2 + H2 * 2;
	nbh  = nbh + H2  * 2;
end	
if ~isempty(HD) 
	nbe  = nbe  + HD * 2;
	nbz2 = nbz2 + HD * 2;
	nbd  = nbd  + HD;
	nbh  = nbh  + HD;
end	
if ~isempty(D2) 
	nbe  = nbe + D2 * 2;
	nbz2 = nbz2 + D2 * 2;
	nbd  = nbd + D2 *2;
end	
if ~isempty(D3) 
	nbe  = nbe + D3 * 2;
	nbz2 = nbz2 + D3 * 2;
	nbd  = nbd + D3 *2;
end	
if ~isempty(He) 
	nbe  = nbe + He * 2;
	nbz2 = nbz2 + He * 4;
end	
if ~isempty(H2O) 
	nbe  = nbe +  (1 + 1./riso) .* H2O * 10;
	nbz2 = nbz2 + (1 + 1./riso) .* H2O * 66;
	nbh  = nbh + H2O * 2;
end	
if ~isempty(Ne) 
	nbe  = nbe + Ne * 10;
	nbz2 = nbz2 + Ne * 100;
end	
if ~isempty(CO) 
	nbe  = nbe + CO * 14;
	nbz2 = nbz2 + CO * 100;
end	
if ~isempty(Ar) 
	nbe  = nbe + Ar * 18;
	nbz2 = nbz2 + Ar * 324;
end	
if ~isempty(CO2) 
	nbe  = nbe + CO2 * 22;
	nbz2 = nbz2 + CO2 * 164;
end	


zeff_lpt    =  nbz2 ./ nbe ;
ndne        =  nbd ./ nbe;

data.riso    = riso;
data.rimp    = rimp;
data.ndne    = ndne;
data.pres_d  = pres_d;
data.pres_he = pres_he;
data.zeff_lpt = zeff_lpt;
data.nbe      = nbe;
data.temps    = temps;
data.nhe      = m4*coef(4);

warning on

if ~all(isfinite(nbe))
	zmain = [];
	nhnd  = [];
	nonc     = [];
	return
end
if (max(nbe) - min(nbe)) < 0.05
	zmain = [];
	nhnd  = [];
	nonc     = [];
	return
end

tmin = tip(min(find(ip > 0.1)));
tmax = tip(max(find(ip > 0.1)));
indchoc = find((temps > tmin) & (temps < tmax));

if all(isfinite(data.riso(indchoc)))
	r = mean(data.riso(indchoc));
	nhnd = r ./ (1-r);
end
if all(isfinite(data.rimp(indchoc)))
	nonc = mean(data.rimp(indchoc));
end
