%  ZUIEDIT_PARAM_DEVICE_FCT  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zuiedit_param_device_fct.m  ->  zuiedit_param_device_fct 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   zuiedit_param_device_fct(action) 
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
function zuiedit_param_device_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_device');

% information pour l'assistant
zuicr(hfig,action) ;
try
   hoc = getfield(h,action) ;
catch
   hoc =[];
end
  
fonction = evalin('base','param.fonction') ;

liste_module =zlist_module;
liste_coef =liste_module.coef;

% selon ation
switch lower(action)


% machine
case {'chng_machine'}
	liste_fct =  liste_module.machine; % zlistdevice;
	hout = zuichgfct('machine',liste_fct,h.machine,h.fct_machine) ;
	zuireset(hoc) ;

case {'para_machine'}
	zuifuninterface('machine') ;
	zuireset(hoc) ;

% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_param_device;
		[hfig,h] = zuiformhandle('ed_param_device');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	zuidisable(h.edit_machine) ;
	zuireset(h.raz) ;
	
% Annulation
case {'annulation','close'}
		zuicloseone(hfig);	

% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig);
	zuisavenonok;
	zuireset(h.validation);
	
otherwise
	warning('action not taken into account')
	   
end

