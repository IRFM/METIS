% ZUIFAITNUM  courte description
%------------------------------------------------------------------------------- 
% fichier : nom_du_fichier -> nom_fonction_principale, nom_fonction_locale 1,... 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   contents = zcontents 
%  
% entrees :  
%	action : tag du uicontrol active
%  
% sorties :  
%   contents  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX
% version  1.7  du  29/09/2001  
%  
% #auto#   Help genere automatiquement
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function zuifaitnum(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesTS');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.etat,'string','')

% selon ation
switch lower(action)

	case {'annulation','close'}
		zuicloseone(hfig);	
	
	case {'init','raz'}
		zuiformvisible(hfig);
		zuiformreset(hfig);
		zuiuploadform(hfig);
		zuireset(h.raz);
	
	case 'validation'
		numchoc = zuidata(h.numchoc);
		try
		    [text,void1,void2,void3]    = tsbase(numchoc,'scprof');
		catch
		    text='';
		end
		if isempty(text)
			set(h.etat,'string','No data for this shot and occurence ...')
			zuidata(h.numchoc,[],[]);
		        zuireset(h.validation);
			helpdlg('You have to use Tprof to generate the profils', ...
			        'No data for this shot and occurence ...');
			return
		end
		chemin = zuidata(h.chemin);
		[s,t]  = unix(sprintf('ls %s',chemin));
		if s ~=0
			set(h.etat,'string','non valid path ...')
			zuidata(h.chemin,'','');
		   zuireset(h.validation);
			return
		end

		zuiformcache(hfig) ;
		zuidownloadform(hfig);
		zuicloseone(hfig);	
		
	otherwise
		warning('action not taken into account')
	   
end

