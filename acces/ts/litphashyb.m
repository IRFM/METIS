function [npar01,npar02,fi1,fi2,cert]=litphashyb(numchoc,tphyb,phyb)

% Fonction : lit dans la base le dephasage et les indices paralleles principaux des antennes hybrides
%
% Appel : [npar01,npar02,fi1,fi2]=litphashyb(numchoc);
%
% Entree : 
% numchoc : numero du choc
%
% Sorties :
% npar01, npar02 : n// central du pic principal de l'antenne pour coupleurs 1 et 2
% fi1, fi2 : dephasage moyen sur le creneau hybride (en degres)
% 
%
% F. Imbeaux CEA-DRFC - juillet 2002
%
% Liste des modifications :
% * 25/02/2003 : leger changement des formules donnant le n// suite au mail d'Annika

% le choc 43541 correspond au premier choc avec C3 et C4
% note :ch/ntt-2009.018 J. Hillairet

if isempty(phyb)
	npar01 = 1.8;
	npar02 = 1.8;
	fi1    = 0 ;
	fi2    = 0;
	return
end


numchoc     = fix(numchoc);
[fi,t,void,cert]      = tsbase(numchoc,'GPHASHYB');
if isempty(t)
	fi1=tsbase(numchoc,'RDPHI1');
	if isempty(fi1)
	   fi1 = 0;
	end
	fi2=tsbase(numchoc,'RDPHI2');
	if isempty(fi2)
		fi2 = 0;
	end
	if  numchoc >= 43541
		% pour C3
		npar01 = 2.01 + (fi1 + 90) ./ 312.5;
		if (npar01> 2.3) | (npar01 <=1.72)
			npar01 = 2.02; 
			disp('Attention, npar01 force a 2.02')
		end
		% pour C4
		npar02 = 1.72 + (fi2 -180) ./ 313.56;
		if (npar02 > 2.0) | (npar02 <=1.43)
			npar02 = 1.72; 
			disp('Attention, npar02 force a 1.72')
		end

	else
		npar01 = 1.84+fi1/200; % pour coupleur ancienne generation
		if (npar01> 3) | (npar01 <=1)
			npar01 = 1.84; %valeur par defaut ancien coupleur
			disp('Attention, npar01 force a 1.84')
		end
		if numchoc > 28000
			npar02 = 2.03+(fi2+90)/310; % pour coupleur nouvelle generation
		else
			npar02 = 1.84+fi2/200; % pour coupleur ancienne generation
		end
		if (npar02 > 3) | (npar02 <=1)
			if numchoc > 28000
				npar02 = 2.03; %valeur par defaut nouveau coupleur
				disp('Attention, npar02 force a 2.03')
			else
				npar02 = 1.84; %valeur par defaut ancien coupleur
				disp('Attention, npar02 force a 1.84')
			end
		end
	end
	fi1      = fi1 .* ones(size(tphyb(:,1)));
	fi2      = fi2 .* ones(size(tphyb(:,1)));
	npar01   = npar01 .* ones(size(tphyb(:,1)));
	npar02   = npar02 .* ones(size(tphyb(:,1)));
	
else
	dt          = gradient(t(:,1));
	indok       = find(dt > 0);
	indbad      = find(dt <=0);
	dtm         = mean(dt(indok));
	dt(indbad)  = dtm;
	t           = cumsum(dt) + t(1,1);

      
	phyb   = interp1(tphyb(:,1),phyb,t,'nearest');
	fi1=fi(:,1);
	fi2=fi(:,2);

	% traduction en n//
	if  numchoc >= 43541
		% pour C3
		npar01 = 2.01 + (fi1 + 90) ./ 312.5;
		ind1 = find(phyb(:,1)>0.3);
		ind11 = find(npar01(ind1)>2.3 | npar01(ind1)<1.72);
		if ~isempty(ind11)
   			npar01(ind1(ind11)) = 2.02 * ones(size(ind11)); %valeur par defaut ancien coupleur
   			disp(['Attention, ',int2str(length(ind11)),' valeurs aberrantes de la phase ont ete remplacees par la valeur par defaut (coupleur 1)'])
		end
	else
		npar01 = 1.84+fi1/200; % pour coupleur ancienne generation
		ind1 = find(phyb(:,1)>0.3);
		ind11 = find(npar01(ind1)>3 | npar01(ind1)<=1);
		if ~isempty(ind11)
   			npar01(ind1(ind11)) = 1.84 * ones(size(ind11)); %valeur par defaut ancien coupleur
   			disp(['Attention, ',int2str(length(ind11)),' valeurs aberrantes de la phase ont ete remplacees par la valeur par defaut (coupleur 1)'])
		end
	end
	%
	if  numchoc >= 43541
		% pour C4
		npar02 = 1.72 + (fi2 -180) ./ 313.56;
		ind2 = find(phyb(:,2)>0.3);
		ind22 = find(npar02(ind2)>2.0 | npar02(ind2)<1.43);
		if ~isempty(ind22)
			npar02(ind2(ind22)) = 1.72 * ones(size(ind22)); %valeur par defaut nouveau coupleur
			disp(['Attention, ',int2str(length(ind22)),' valeurs aberrantes de la phase ont ete remplacees par la valeur par defaut (coupleur 2)'])
		end
	elseif numchoc > 28000
 	  	npar02 = 2.03+(fi2+90)/310; % pour coupleur nouvelle generation
		ind2 = find(phyb(:,2)>0.3);
		ind22 = find(npar02(ind2)>3 | npar02(ind2)<=1);
		if ~isempty(ind22)
			npar02(ind2(ind22)) = 2.03 * ones(size(ind22)); %valeur par defaut nouveau coupleur
			disp(['Attention, ',int2str(length(ind22)),' valeurs aberrantes de la phase ont ete remplacees par la valeur par defaut (coupleur 2)'])
		end
	else
  	 	npar02 = 1.84+fi2/200; % pour coupleur ancienne generation
		ind2 = find(phyb(:,2)>0.3);
		ind22 = find(npar02(ind2)>3 | npar02(ind2)<=1);
		if ~isempty(ind22)
			npar02(ind2(ind22)) = 1.84 * ones(size(ind22)); %valeur par defaut ancien coupleur
			disp(['Attention, ',int2str(length(ind22)),' valeurs aberrantes de la phase ont ete remplacees par la valeur par defaut (coupleur 2)'])
		end
	end


	npar01   = interp1(t,npar01,tphyb(:,1),'nearest');
	npar02   = interp1(t,npar02,tphyb(:,1),'nearest');
	fi1      = interp1(t,fi1,tphyb(:,1),'nearest');
	fi2      = interp1(t,fi2,tphyb(:,1),'nearest');
	%[fi1]=tsbase(numchoc,'RDPHI1');   % dephasage moyen sur creneau hybride en degres
	%[fi2]=tsbase(numchoc,'RDPHI2');   % dephasage moyen sur creneau hybride en degres
end
