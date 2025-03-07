% cree un profil de ti a partir des donnees de tbragg
% le profil est homotetique a ne
function [ti,sbrag] = ztibragg(numchoc,temps,x,ne,te,forme,pfci_ion,phys,ae0)

renorme  = 1; % facteur de renormalisation de ti 

% corrction de ae0

indnok = find(~isfinite(ae0) | (ae0 <= 0)| (ae0 > 1));
indok  = find(isfinite(ae0) & (ae0 > 0)& (ae0 <= 1));
if ~isempty(indnok)
	ae0(indnok) = mean(ae0(indok));
end
ae0 = ae0 * ones(1,length(x));

% lecture de donnees de tbragg
times = temps;
[tivb,ttivb]=tsbase(fix(numchoc),'stibrag');
if ~isempty(tivb)
  [brag,sbrag]   = cgcgettrait(fix(numchoc),'tbragg@');
else
  brag.ti = [];
  sbrag.ti = [];
end
% loi d'echelle pour le points manquants
tisl                   = 0.6234 .* (max(ne')'/1e19) .^0.2846 .* (max(te')'/1e3) .^ 0.6538;

% ti0 a partir de dbrag
if isempty(brag.ti)
   brag.ti =NaN;
end
if ~all(~isfinite(brag.ti))
 		ti0  	   = brag.ti .* 1e3;
		indnok   = find( (temps < min(sbrag.ti.temps)) |  ...
		               (temps > max(sbrag.ti.temps)) | ...
							~isfinite(brag.ti));
	   indok    = find( (temps >= min(sbrag.ti.temps)) &  ...
		                 (temps <= max(sbrag.ti.temps)) & ...
							  isfinite(brag.ti));
      if isempty(indok)
         % si hors intervalle
         indok = find(isfinite(brag.ti));
      end
		indokb   = find((sbrag.ti.temps >=	min(temps(indok))) &	 ...
				  			(sbrag.ti.temps <=	max(temps(indok))) & ...
							 isfinite(sbrag.ti.data));
		% renormalisation		
		if length(indok)> 10	
		   if forme   > 0
				%rappsl   = mean(tisl(indok) ./ (max(te(indok,:)')'/1e3));
				%rapp     = mean(brag.ti(indok) ./ brag.te(indok));
				renorme  = mean(max(te(indok,:),[],2)/1e3) ./ mean( sbrag.te.data(indokb));
		   end
		   tisl     =  mean(sbrag.ti.data(indokb)) ./ mean(tisl(indok)).* tisl;
			
			% securite sur tisl 
			% securite ti0 < te0 si pas de fci
			indi = find(mean(te(:,1:3)/1e3,2) < mean(sbrag.te.data(indokb)));
			if ~isempty(indi)
				tisl(indi) = max(tisl(indi), mean(sbrag.ti.data(indokb)));
			end
		end
		
		
		if ~isempty(indnok)
		    ti0(indnok) = tisl(indnok) .* 1e3;
		end
		% securite
		if any(~isfinite(ti0))
			ti0(~isfinite(ti0)) = tisl(~isfinite(ti0)) .* 1e3;
		end
else
		disp('Pas de donnees tbragg')
		ti0   = tisl .* 1e3;
end

% anti spike
if forme > 0
	ti0   = renorme .* smooth(ti0,3);
else
	ti0   = smooth(ti0,3);
end

% securite ti0 < te0 si pas de fci
indi = find(pfci_ion <= 0.3);
if ~isempty(indi)
	ti0(indi) = min(ti0(indi),te(indi,1));
end

% correction du bord pour le profil de ne 
pp          = [ -0.2867   27.7006 -622.8086];
nbar        = trapz(x,ne,2);
nebord      = max(3e17,exp(polyval(pp,log(nbar))));
ne          = ne - ne(:,end) * ones(1,size(ne,2)) + nebord * ones(1,size(ne,2));
nbarnew     = trapz(x,ne,2);
ne          = ne  .* ((nbar ./ nbarnew) * ones(1,size(ne,2)));
tia         = max(min(200,ti0 .* nebord ./ ne(:,1)),30);
ti0         = max(ti0,1.5 .* tia);

% calcul du profil de ti
if (forme == 0) | (forme == 1)
	tif         = 0.5 .* ( ne ./ (max(ne')'* ones(1,length(x))) + ...
 											 te ./ (max(te')'* ones(1,length(x))));
elseif forme  ==  -2
	tif         = te ./ (max(te')'* ones(1,length(x)));
elseif forme == -1
	tif         =  sqrt(ne ./ (max(ne')'* ones(1,length(x))) .* te ./ (max(te')'* ones(1,length(x))));
else 
	tif         = ne ./ (max(ne')'* ones(1,length(x)));
end 	
								
tif         = tif - ( tif(:,end) * 	ones(1,length(x)));
tif         = tif ./ (tif(:,1) * 	ones(1,length(x)));

for k =1:size(tif,1)
			tif(k,:)=zmonotone(x,tif(k,:));
end

ti  			= ((ti0 - tia)* ones(1,length(x))) .* tif + tia * ones(1,length(x));


% dans le cas ou tbragg est absent -> correction complete de ti 
if  isempty(brag.ti) & (forme >= 0)
	xx          = ones(size(ti,1),1) * x;
	tauex       = 6.14e14 ./ ne .* phys.mp ./ 2 ./ phys.me .* (te ./1e3).^1.5;
	factei      = 1.5 .* ae0 .* ne ./ tauex .* phys.e;
	qi          =   trapz(x,factei .* ti .* xx ,2);
	qe          =   trapz(x,factei .* te .* xx,2); 
	indnok      =   find(qe <= 0);
	qe(indnok)  =   1;
	g           =   (qe+ pfci_ion .* 1e6) ./ qi;
	indok       =   find((qe > 0) & (g > 0.3));
	indnok      =   find((qe <= 0) | (g <= 0.3));
	gm          =   mean(g(indok));
	if ~isempty(indnok)
			g(indnok)   = gm;
	end
	gmax        = te(:,1) ./ (ti(:,1) + eps);
	indnok      = find(ti(:,1) == 0);
	if ~isempty(indnok)
			gmax(indnok)   = 1;
	end
	g           = min(g,gmax);
	g           = smooth(g,7);
	ti          =   ti  .* (g * ones(1,size(ti,2)));


end
										  
% securite sur Ti/Te dans la partie externe du profil
xx          = ones(size(ti,1),1) * x;
tec         = te + (max(tia,te(:,end))- te(:,end)) * ones(1,length(x));
ind0        = ti <= tec;
tauex       = 6.14e14 ./ ne .* phys.mp ./ 2 ./ phys.me .* (tec ./1e3).^1.5;
factei      = 1.5 .* ae0 .* ne ./ tauex .* phys.e;
qis         = - trapz(x,factei .* ti .* xx .* ind0,2);
qes         =   trapz(x,factei .* tec .* xx .* ind0,2); 
qe          =   trapz(x,factei .* tec .* xx .* (~ ind0),2); 
indnok      =   find(qe < qes);
indok       =   find(qe >= qes);
qe(indnok)  =   1;
g           =   (qe + qes + qis + pfci_ion .* 1e6) ./ qe;
g(indnok)   =   1;
g(g < 1)    =   1;
g           =   smooth(g,7);
g(g < 1)    =   1;
g           =   g * ones(1,size(ti,2));
timax       =   ti .* ind0 + g .* tec .* (~ind0);
ti(indok,:) =   min(timax(indok,:),ti(indok,:));

