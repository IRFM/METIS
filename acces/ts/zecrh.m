function [cons,polar,phi_pol,phi_tor,temps,cert] = zecrh(choc,temps)


% selon le choc
cert.cert   = -1;
cert.vers   =  2;
cert.date   = datestr(now,2);
cert.heure  = datestr(now,13);
cert.unix   = 's';
cert.uniy   = 'MW';
cert.uniz   = ''  ;
cert.nomdon = 'ZPECRH';
cert.type   =  'Z';
cert.commentaire = 'Version provisoire du traitement tfce a venir';
if nargin < 2
	temps =[];
end
cons    = zeros(size(temps,1),3);
polar   = [NaN,NaN,NaN];
phi_pol = [NaN,NaN,NaN];
phi_tor = [NaN,NaN,NaN];


% patch pour choc a probleme
switch choc
case {30000,30767,31317,31314}
   disp('Probleme de base temps - donnees inutilisable');
   return 
end

% selon le choc
pecrh_mat = load(fullfile(fileparts(which('zecrh')),'Pecrh'));
pecrh_mat = pecrh_mat.Pecrh;
indch     = find(choc == fix(pecrh_mat(:,1)));
if isempty(indch)

   % lecture des puissances
   [pecrh1,tpecrh1]       = tsbase(fix(choc),'sika1');
   if isempty(pecrh1)
   	[pecrh1,tpecrh1]       = tsbase(fix(choc),'sicata1');
   end
   if isempty(pecrh1)
	   [pecrh1,tpecrh1]       = tsbase(fix(choc),'spinca1');
   end
   [pecrh2,tpecrh2]       = tsbase(fix(choc),'sika2');
   if isempty(pecrh2)
	   [pecrh2,tpecrh2]       = tsbase(fix(choc),'sicata2');
   end
   if isempty(pecrh2)
		[pecrh2,tpecrh2]       = tsbase(fix(choc),'spinca12');
   end
%
% patch pour base de temps mauvaise en fin de choc
% nouveau programme en cours de validation
% choc concerne, journee du 17 septembre 2003
% V. Basiuk 
%
   indval = min(find(gradient(tpecrh1)<0));
   if ~isempty(indval)
      tpecrh1(indval:end)=[];
      pecrh1(indval:end)=[];      
   end
   indval = min(find(gradient(tpecrh2)<0));
   if ~isempty(indval)
      tpecrh2(indval:end)=[];
      pecrh2(indval:end)=[];      
   end
   
   if ~isempty(pecrh1)
      indana=find(pecrh1 < max(pecrh1)*0.5);
%      if ~isempty(find(tpecrh1< 0.2))
%	      pecrh1 = pecrh1- mean(pecrh1(tpecrh1< 0.2));
%	      pecrh1(pecrh1 <= (10 * std(pecrh1(tpecrh1< 0.2)))) = 0;
%      end
      if ~isempty(indana)
	      pecrh1 = pecrh1- mean(pecrh1(indana));
	      pecrh1(pecrh1 <= (10 * std(pecrh1(indana)))) = 0;
      end

      pecrh1 = pecrh1 - min(pecrh1);
      pecrh1 = pecrh1 ./ max(pecrh1+eps) .* 0.300;
      [tpecrh1,pecrh1] = antir(tpecrh1,pecrh1);
      ind = min(find(gradient(tpecrh1(10:end))<0));
      if ~isempty(ind)
        tpecrh1(ind:end) = [];
        pecrh1(ind:end) = [];
      end
      if isempty(temps)
		   temps = tpecrh1;
	   else
		   pecrh1 = interp1(tpecrh1,pecrh1,temps,'nearest');
	   end
   end
   if ~isempty(pecrh2)
      indana=find(pecrh2 < max(pecrh2)*0.5);
      if ~isempty(indana)
	      pecrh2 = pecrh2- mean(pecrh2(indana));
	      pecrh2(pecrh2 <= (10 * std(pecrh2(indana)))) = 0;
      end
%      if ~isempty(find(tpecrh2< 0.5))
%	      pecrh2 = pecrh2- mean(pecrh2(tpecrh2< 0.5));
%	      pecrh2(pecrh2 <= (10 * std(pecrh2(tpecrh2< 0.5)))) = 0;
%      end
      pecrh2 = pecrh2 - min(pecrh2);
      pecrh2 = pecrh2 ./ max(pecrh2+eps) .* 0.450;
	   [tpecrh2,pecrh2] = antir(tpecrh2,pecrh2);
      ind = min(find(gradient(tpecrh2(10:end))<0));
      if ~isempty(ind)
        tpecrh2(ind:end) = [];
        pecrh2(ind:end) = [];
      end
      if isempty(temps)
		   temps = tpecrh2;
	   else
		   pecrh2 = interp1(tpecrh2,pecrh2,temps,'nearest');
	   end
   end

   % recherche des informations
   chemin   = fileparts(which('zecrh'));
   [s,r] = unix(sprintf('grep ^%d %s',choc,fullfile(chemin,'dataecrh.dat')));
   if s == 0
      [nbc,r]       = strtok(r);
      [polar1,r]    = strtok(r);
      [puiss1,r]    = strtok(r);
      [pol1,r]      = strtok(r);
      [tor1,r]      = strtok(r);
      [polar2,r]    = strtok(r);
      [puiss2,r]    = strtok(r);
      [pol2,r]      = strtok(r);
      [tor2,r]      = strtok(r);
	   polar         = [str2num(polar2),str2num(polar1),NaN];
	   phi_pol       = [str2num(pol2),str2num(pol1),NaN];
	   phi_tor       = [str2num(tor2),str2num(tor1),NaN];
	   puiss1        = str2num(puiss1) / 1000;
	   if ~isempty(puiss1) & isfinite(puiss1)
		   pecrh1 = pecrh1 - min(pecrh1);
		   pecrh1 = pecrh1 ./ max(pecrh1) .*  puiss1;
	   end
	   puiss2        = str2num(puiss2) / 1000;
	   if ~isempty(puiss2) & isfinite(puiss2)
		   pecrh2 = pecrh2 - min(pecrh2);
		   pecrh2 = pecrh2 ./ max(pecrh2) .*  puiss2;
	   end
	 
   else
   	disp('ECRH : pas d''information pour ce choc');
    cert.cert = -2;
   end
   if ~isempty(pecrh1) & ~isempty(pecrh2)
	   cons = cat(2,pecrh2,pecrh1,0 .* temps);
   elseif ~isempty(pecrh2)
	   cons = cat(2,pecrh2,zeros(size(pecrh2)),0 .* temps);
   elseif ~isempty(pecrh1)
	   cons = cat(2,zeros(size(pecrh1)),pecrh1,0 .* temps);
   elseif ~isempty(temps)
	   cons = zeros(size(temps,1),3);
   end
   cons(~isfinite(cons)) = 0;
   cert.size = size(cons);
else
   % parametres
   facteur = pecrh_mat(indch,2);
   numero  = pecrh_mat(indch,3);
   ang_pol = pecrh_mat(indch,4);
   ang_tor = pecrh_mat(indch,5);
   % lecture des puissances
   [pecrh,tpecrh]= tsbase(choc,sprintf('spia%d',numero));
   if isempty(temps)
		temps = tpecrh;
	else
		pecrh = interp1(tpecrh,pecrh,temps,'nearest');
	end
   % normaliation 
   pecrh = pecrh .* facteur;
   
   switch numero
   case 1
 	   cons = cat(2,zeros(size(pecrh)),pecrh,zeros(size(pecrh)));
      polar   = [NaN,2,NaN];
      phi_pol = [NaN,ang_pol,NaN];
      phi_tor = [NaN,ang_tor,NaN];
   case 2
 	   cons = cat(2,pecrh,zeros(size(pecrh)),zeros(size(pecrh)));
      polar   = [2,NaN,NaN];
      phi_pol = [ang_pol,NaN,NaN];
      phi_tor = [ang_tor,NaN,NaN];
   end

end


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
