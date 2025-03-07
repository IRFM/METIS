% ZSAMPLE reechantillone un signal sur le temps et l'espace
%----------------------------------------------------------
% fichier zsample.m ->  zsample
%
%
% fonction Matlab 5 :
%
% Cette fonction gere le filtrage et le reechantillonage d'un signal 
% ou d'un groupe de signaux sur le temps et/ou l'espace.
% Elle permet de debruite le signal via un filtrage par ondelette.
% Elle peut lisser les groupe via un filtrage par svd.
% Elle reechantillone sur le temps uniquement si les vecteurs de temps
% sont different. Il en va de meme pour le reechatillonage spatiale. 
% Les valeurs hors intervalle sont mise a une valeur par defaut.
%
% syntaxe  :
%  
%  1 - pour un signal simple :
%  
%       signal = zsample(signal,temps_signal,temps,{param})
%    
%  2 - pour un groupe :
%      
%       [signal,kmax] = zsample(signal,temps_signal,rho_signal,temps,rho,{param})
%
%
% entrees :
%
%       signal         =   signal a reechantilloner   {size(signal) =[N,M]}
%       temps_signal   =   vecteur temps du signal a reechantilloner  {size(temps_signal) = [N,1]}
%       rho_signal     =   vecteur ou matrice espace du signal a reechantilloner {size(rho_signal) = [1,M] 
%                                                                               ou size(rho_signal) =[N,m]}
%       temps          =   vecteur temps du signal en sortie {size(temps) = [NN,1]}
%       rho            =   vecteur ou matrice espace du signal en sortie {size(rho) = [1,MM] ou size(rho) = [NN,MM]}
%       param          =   structure parametre optionnel de la fonction
%                       
% sorties :
% 
%     signal           =   signal reechantillone {size(signal) =[NN,MM]}
%     kmax             =   nombre de coefficient non nuls garder apres la svd
%     
% 
% parametre de la fonction :
% 
%     param.ondelette  =  active (1) ou non (0) le debruitage par ondelette
%     param.energie    =  gere le filtrage par svd :
%                          *  si <1, regle la coupure sur la variation denrgie de voie a voie.
%                             cette variation doit etre < param.energie a la coupure.
%                          
%                          *  si == 1, pas de svd
%                          
%                          * si >1, fixe le nombre de coeficients conserves
%                          
%     param.defaut.temps     =  valeur par defaut mise hors de l'intervalle de temps apres reechantillonage temporel
%     param.defaut.espace    =  valeur par defaut mise hors de l'intervalle d'espace apres reechantillonage spatiale
%     param.defaut.inf       =  valeur par defaut mise a la place des NaN et Inf  (si non vide)
%     param.plus             =  si 1 le signal est defini positif .
%     
%  valeur par defaut :
%    1 - pour un signal simple :
%           param.ondelette        = 1;
%           param.defaut.temps     = NaN;
%           param.defaut.espace    = 0;
%           param.defaut.inf       = [];
%           param.plus             = 0;
%           
%    2 - pour un groupe de signaux :
%           param.ondelette         = 0;
%           param.energie           = 0.01;
%           param.defaut.temps     = NaN;
%           param.defaut.espace    = 0;
%           param.defaut.inf       = [];
%           param.plus             = 0;
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 08/10/2004.
% 
% 
% liste des modifications : 
%
% * 08/10/2004 -> ajout securite sur dt =0
%
%--------------------------------------------------------------
%
function [signal,kmax] = zsample(signal,temps_signal,rho_signal,temps,rho,param)


% sortie par defaut
kmax=[];

% test des entrees
if nargin <3
	error('Nombre d''arguments incorrect !');
end
% squeeze de la matrice
signal =squeeze(signal);

% securite dt =0 
indnok = find(diff(temps_signal) <=0);
if ~isempty(indnok)
	while ~isempty(indnok) 
		if size(rho_signal,1) == size(temps_signal,1)
			rho_signal(indnok,:) =[];
		end
   		signal(indnok,:) = [];
   		temps_signal(indnok,:) =[];
		indnok = find(diff(temps_signal) <=0);
	end
end


% choix du mode 1 D ou 2D
ss=size(signal);
if length(ss) >2
	error('Signal doit etre un vecteur ou une matrice 2D de type signal(temps,espace)');
end
if all(ss ==1)
	return
end

if any(ss==1)
	dimmode =1;
	if nargin == 4
		param = temps;
	elseif nargin ~= 6
		param=[];
	end
	temps = rho_signal;
	clear rho_signal
	
else
	dimmode =2;
	if nargin < 5
		error('Il faut les 5 arguments pour les signaux 2d');
	end
	if nargin <6
		param=[];
	end
end

% parametres automatique si param est vide
if isempty(param)
	if dimmode == 1
		param.ondelette =1;      % debruitage
	else
		param.ondelette =0;
		param.energie =0.01;
	end

  param.defaut.temps     = NaN;
  param.defaut.espace    = 0;
  param.defaut.inf       = [];
  param.plus             = 0;
end

if (param.plus  == 1)&(~isempty(signal))
	ind = find(signal<0);
	if ~isempty(ind)
		signal(ind) = zeros(1,length(ind));
	end
end

% cas du signal 1D
if dimmode ==1
	% on redesse le signal si necessaire 
	signal = signal(:);
	
	if param.ondelette == 1
		% debruitage par ondelette
		lev=min(7,floor(log(ss(1))/log(2))-1);   % profondeur de la decompostion en ondelette
		signal=zwfilter(signal',lev)'; % fliltrage par ondelette
	end
		
	% reechantillonage
	if ~isequal(temps_signal,temps)
		signal = tsample(signal,temps_signal,temps,'fen');
		indice = find( ~((temps>=min(temps_signal))&(temps<=max(temps_signal))));
		if ~isempty(indice)
  		     signal(indice) = param.defaut.temps .* ones(1,length(indice));
		end
	end
	
	% retrait des infinis
	indice =find(~isfinite(signal));
	if ~isempty(indice) & ~isempty(param.defaut.inf)
		signal(indice) = param.defaut.inf .* ones(1,length(indice));
	end
	
	% fin du cas 1D
else
	% cas du signal 2D
	if param.ondelette == 1
		% debruitage par ondelette
		lev=min(7,floor(log(ss(1))/log(2))-1);   % profondeur de la decompostion en ondelette
		for k =1:ss(2)
			signal(:,k)=zwfilter(signal(:,k)',lev)'; % fliltrage par ondelette
		end
	end


	% mise en forme des rho
	if any(size(rho_signal)==1)
		rho_signal = ones(size(temps_signal,1),1) * rho_signal;
	end
	if any(size(rho)==1)
		rho = ones(size(temps,1),1) * rho;
	end
	
	
	% test sur l'energie a conserver pour la svd
	if param.energie <1
		svdmode = 1;    % svd critere energie
	else
		param.energie=fix(param.energie);
		if param.energie >1
			svdmode = 2; % nombre de voie gardee
		else
			svdmode = 0; % pas de svd
		end
	end
	
	% filtrage svd
	if svdmode  ~= 0
		% avec svd
		% decomposition svd 
		[u,s,v]=svd(signal,0);
		
		% calcul de la coupure
		if svdmode == 1 
			% calcul de l'energie
			d =diag(s);
			energie = sqrt(cumsum(d) .^ 2);
  			energie = energie ./max(energie);
  			de=abs(pdederive(1:length(d),energie,2,2,1,2));
  			kmax = min(find( de <= param.energie)); 
		else
			kmax = fix(param.energie);
		end
		
		% coupure de la svd
		u=u(:,1:kmax);
		v=v(:,1:kmax);
		s=s(1:kmax,1:kmax);
		
		% reconstitution svd
	   signal=u*s*v';		
	        
	end
	  
	% reechantillonage 
	% sur le temps
	if ~isequal(temps_signal,temps)
		signal = tsample(signal,temps_signal,temps,'fen');
		rho_signal = tsample(rho_signal,temps_signal,temps,'fen');
	end
	% sur l'espace
	if ~isequal(rho_signal,rho)
		signal = tsplinet(rho_signal,signal,rho);
	end
	
   % valeur hors intervalle     
	indice = find( ~((rho>=(min(rho_signal')'*ones(1,size(signal,2))))&(rho<=(max(rho_signal')'*ones(1,size(signal,2))))));
	if ~isempty(indice) 
  		 signal(indice) = param.defaut.espace .* ones(1,length(indice));
	end
	indice = find( ~((temps>=min(temps_signal))&(temps<=max(temps_signal))));
	if ~isempty(indice)
  		  signal(indice,:) = param.defaut.temps .* ones(length(indice),size(signal,2));
   end
  		 
	% retrait des infinis
	indice =find(~isfinite(signal));
	if ~isempty(indice) & ~isempty(param.defaut.inf)
		signal(indice) = param.defaut.inf .* ones(1,length(indice));
	end
end

% signal positif en sortie aussi
if (param.plus  == 1)&(~isempty(signal))
	ind = find(signal<0);
	if ~isempty(ind)
		signal(ind) = zeros(1,length(ind));
	end
end
