function zuiedit_param_postrait_calc_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_postrait_calc');

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

% post
case {'chng_post'}
	liste_fct =  liste_module.post; %zlistpost;
	hout = zuichgfct('post',liste_fct,h.post,h.fct_post) ;
	zuireset(hoc) ;

case {'para_post'}
	zuifuninterface('post') ;
	zuireset(hoc) ;
		
case {'mode_post'}
	value = zuidata(h.mode_post) ;
	if value==2
		zuienable(h.edit_post) ;
	else
		zuidisable(h.edit_post) ;
	end
		
case {'edit_post'}
	zuiedit_mode('data.mode.post',{0,2,3},{'no calculation','calculated','copy previous value'}) ;
	zuireset(hoc) ;

% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_asser_calc;
		[hfig,h] = zuiformhandle('ed_asser_calc');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	zuidisable(h.edit_asser) ;
	zuireset(h.raz) ;
	
% Annulation
case {'annulation','close'}
		zuicloseone(hfig);	

% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig);
	var = evalin('base','data.mode.asser') ;
	val = zuidata(h.mode_asser) ;
	if val~= 3
		zassignin('base','data.mode.asser',ones(size(var))*val) ;
	end
	zuisavenonok;
	zuireset(h.validation);
	
otherwise
	warning('action not taken into account')
	   
end
