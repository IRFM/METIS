function [temps,consigne,remplissage,ispi,pompage,contenu,gaz,certout] = zchimie(choc,tintin)

if nargin < 2
	tintin =[];
end

% constantes de travail
zgaz   = [1 1 1 2 2 6 7 10 18];
agaz   = [1 2 3 3 4 NaN NaN NaN NaN];
fgaz   = [2 2 2 1 1 1 2 1 1];
strgaz = {'H','D','T','He3','He','C','N','Ne','Ar'}; 
gaz.zgaz  = zgaz;
gaz.agaz  = agaz;
gaz.fgaz  = fgaz;
gaz.strgaz  = strgaz;

temps       = [];
consigne    = [];
remplissage = [];
ispi        = [];
pompage     = [];
contenu     = [];

avo     =   6.02214199e23;   
npal    = avo / 22.4 / 1.013e5;
npam3   = 1000 * npal;

% lecture des donnees
% type de gaz injectes
H2_onoff  = tsmat(choc,'EXP=T=S;GAZ;IMPH2');
Ne_onoff  = tsmat(choc,'EXP=T=S;GAZ;IMPNE');
D2_onoff  = tsmat(choc,'EXP=T=S;GAZ;IMPD2');
Ar_onoff  = tsmat(choc,'EXP=T=S;GAZ;IMPAR');
N2_onoff  = tsmat(choc,'EXP=T=S;GAZ;IMPN2');
CD4_onoff = tsmat(choc,'EXP=T=S;GAZ;CD4');
HE3_onoff = tsmat(choc,'EXP=T=S;GAZ;HE3');

% debit injecte
[debit,tdebit,void,cert]  = tsbase(choc,'GDEBIT');
certout.gdebit = cert;
if isempty(debit)
	return
end
% pre traitement
[tdebit,debit] = antir(tdebit,debit);
debit(debit == -1)  = 1.1 .* max(debit(:));
debit(debit < 0) = 0;
for k = 1:size(debit,2)
	debit(:,k) = medfilt1(debit(:,k),11);
end


% ouverture des vannes
[vanne,tvanne,void,cert]   = tsbase(choc,'GTENSION');
certout.gtension = cert;

[tvanne,vanne] = antir(tvanne,vanne);
vanne          = vanne - ones(size(tvanne)) * mean(vanne(tvanne < -20,:),1);
maskv = (vanne > 10);
vanne_onoff = (sum(maskv,1) > 10);

% correction des debits
debit       = debit .* (ones(size(tdebit)) * vanne_onoff);
% recherche de la debut du choc
inddebgaz     = min(find(all(debit < 0.1,2)));
debit(1:inddebgaz,:) = 0;

% association debit%x -> gaz
deb.h       = zeros(size(tdebit));
deb.d       = zeros(size(tdebit));
deb.c       = zeros(size(tdebit));
deb.he      = zeros(size(tdebit));
deb.he3     = zeros(size(tdebit));
deb.ne      = zeros(size(tdebit));
deb.n       = zeros(size(tdebit));
deb.ar      = zeros(size(tdebit));
 
deb.he      = sum(debit(:,3:4),2);
deb.d       =  2 .* sum(debit(:,[5,6,7,10]),2) .* npam3;
if CD4_onoff 
 	deb.d      =  deb.d + 4 .* sum(debit(:,[8,9]),2) .* npam3;
 	deb.c      =  sum(debit(:,[8,9]),2) .* npam3;
else
 	deb.d      =  deb.d + 2 .* sum(debit(:,[8,9]),2) .* npam3;
end
if HE3_onoff
	deb.he3    = debit(:,13) .* npam3;
else
	deb.he     = debit(:,13) .* npam3;
end

if H2_onoff
   deb.h      =  2 .* sum(debit(:,[1,2,11,12]),2) .* npam3;
elseif Ne_onoff
   deb.ne     =  sum(debit(:,[1,2,11,12]),2) .* npam3;
elseif D2_onoff
   deb.d      =  deb.d +2 .* sum(debit(:,[1,2,11,12]),2) .* npam3;
elseif Ar_onoff
   deb.ar     =  sum(debit(:,[1,2,11,12]),2) .* npam3;
elseif N2_onoff
   deb.n     =  2 .* sum(debit(:,[1,2,11,12]),2) .* npam3;
end

% ici il faut ajouter l'ispi
try
	[inj,tinj] = zispi(choc);
catch
	inj =[];
end

if isempty(inj)
   % volume injecte (ispi)
%   [ mano,tmano,void,cert]  = tsbase(choc,'GMANO%13');
   [ mano,tmano,void,cert]  = tsbase(choc,'GISPI%5');
   certout.gmano = cert;
   [tmano,mano] = antir(tmano,mano);
   indp         = min(find((mano(1:(end-1)) > 1000) & (mano(2:end) < 500)));
   if isempty(indp)
	indp = max(find(tmano < 0));
   else
	indp = indp +1;
   end
   mano         = mano(indp:end) - mean(mano(indp:(indp+10)));
   tmano        = tmano(indp:end);
   %manozero     = mean(mano(tmano < 0));
   %mano         = mano - manozero;
   manoscale  = -mean(mano(mano <0));
   mano(mano < 0) = 0;
   mano         = medfilt1(mano,11);
   if (any(mano  > (3*manoscale))) & (std(mano) < manoscale) 
	disp('choc avec ispi')
	dmanodt      = pdederive(tmano,mano,2,2,1,1);
	seuil        = 3 .* max(std(dmanodt(dmanodt >= 0 )),std(dmanodt(dmanodt <=0)));
	mask         = (dmanodt > seuil);
	maskp        = cat(1,0,~mask);
	maskm        = cat(1,mask,0);
	indtir       = find((maskp.*maskm));
	indtir(end+1) = length(tmano);
        if indtir > 2
	    manotir      = zeros(length(indtir)-1,1);
	    tempstir     = zeros(length(indtir)-1,1);
	    for k = 1:(length(indtir)-1)
		manotir(k)  = mean(mano(indtir(k):indtir(k+1)));
		tempstir(k) = tmano(indtir(k));
	    end
	    % signal de sortie ispi
	    ispi.times = tempstir;
	    ispi.nbatomes = abs(gradient(manotir)) .* npal .* 2;
	    %correction des erreur de mesures
            [manotir,tempstir] = antir(manotir,tempstir);   
            if length(manotir) >= 2
 	       contenuispi        = interp1(tempstir,manotir,tdebit,'linear');
	       contenuispi(~isfinite(contenuispi)) = 0;
	       fsh   = 1./ 2 ./ mean(diff(tdebit));
	       [bf,af] = butter(7,2./fsh);
	       contenuispi = medfilt1(abs(filtfilt(bf,af,contenuispi)),11);
	       debispi            = pdederive(tdebit,contenuispi,2,2,1,1) .* npal .* 2;
	       debispi(~isfinite(debispi)) = 0;
	       debispi(debispi <=0) = 0;
	       debispi(tdebit < min(tempstir)) = 0;
	       [v1,v2,v3,debispi] = meanfilt1(debispi,11);
	       debispi            = min(5e21,debispi - min(debispi));
	       deb.d              = deb.d + debispi';
            end
      else
	ispi.times    = [];
	ispi.nbatomes = [];
      end
   else
	ispi.times    = [];
	ispi.nbatomes = [];
   end
else
	ispi.times    = tinj;
	ispi.nbatomes =  2 .* inj;
     
end      
% consigne
consigne = cat(2,deb.h,deb.d,0 .* tdebit,deb.he3,deb.he,deb.c,deb.n,deb.ne,deb.ar);

% remplissage
remplissage  = trapz(tdebit(tdebit<0),consigne(tdebit < 0,:),1);

% contenu
contenu      = cumtrapz(tdebit,consigne,1);

%  la vitesse de pompage
pomp_lpt       = tsmat(choc,'EXP=T=S;POMPAGE;LPT');
pomp_caisson   = tsmat(choc,'EXP=T=S;POMPAGE;CAISSON');
pomp_divers    = tsmat(choc,'EXP=T=S;POMPAGE;DIVERS');
if ~isempty(pomp_lpt) & ~isempty(pomp_caisson) & ~isempty(pomp_divers)
  vpomp          = (0.5 .* npam3 / 10) .* (sum(pomp_lpt) + ...
                  0.3 .* sum(pomp_caisson(1:2)) +  ...
						0.1 .* sum(pomp_divers));
else
   vpomp = 0;
end

% correction du contenu
pompage        = vpomp .* ((ones(size(contenu,1),1) * fgaz) .* contenu ./  ...
                          (1 + sum(contenu,2) * ones(1,size(contenu,2))));

% temps du choc
if isempty(tintin)
	[ip,tip] = tsbase(choc,'sipmes');
	imin     = min(find(ip > 0.1));
	imax     = max(find(ip > 0.1));
	temps    = tip(imin:imax);
else
	temps    = tintin;
end
consigne    = interp1(tdebit,consigne,temps,'nearest');
contenu     = interp1(tdebit,contenu,temps,'nearest');
pompage     = interp1(tdebit,pompage,temps,'nearest');



function [t,s] = antir(t,s)

if isempty(t)
	return
end

nb = 30;
t  = mean(t,2);
indr = find(diff(t) <= 0);
while ~isempty(indr) & (nb > 0)
	indr = find(diff(t) <= 0);
	if ~isempty(indr)
		s(indr,:) = [];
		t(indr,:) = [];
	else 
   	break
	end
	nb  = nb - 1 ;
end
