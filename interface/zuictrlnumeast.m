%  ZUICONTROLNUMEAST  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zuictrlnumeast.m  ->  zuicontrolnumeast 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   zuicontrolnumeast(action) 
%  
% entrees :  
%  action = 
%  
% sorties :  
%  
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  4.1  du  08/04/2008  
%  
%*@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function zuicontrolnumeast(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesEAST');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.etat,'string','')

% selon ation
switch lower(action)

	case 'numchoc'
		numchoc = zuidata(h.numchoc);
		if isempty(numchoc)
			set(h.etat,'string','enter a shot number ...')
			zuidata(h.numchoc,1,[]);
		   zuireset(h.validation);
			return
		end
		if numchoc < 0
			set(h.etat,'string','numchoc < 0 ...')
			zuidata(h.numchoc,1,[]);
		   zuireset(h.validation);
			return
		end
		
	case 'chemin'
		chemin = zuidata(h.chemin);
		[s,t]  = unix(sprintf('ls %s',chemin));
		if s ~=0
			set(h.etat,'string','non valid path ...')
			zuidata(h.chemin,'','');
		end
		
	
	otherwise
		warning('action not taken into account')
	   
end

