function zuiedit_param_postrait_param_fct(action)

%  
[hfig,h] = zuiformhandle('ed_postrait_param') ;
if ~ishandle(hfig)
	return
end
hoc = getfield(h,action) ;

switch lower(action)

% Annulation
	case {'annulation','close'}
		zuicloseone(hfig);	
	
% raz
	case {'init','raz'}
		zuiformvisible(hfig) ;
		zuiformreset(hfig) ;
		zuiuploadform(hfig) ;
		zuireset(h.raz) ;
		
% Validation
	case 'validation'
		zuidownloadform(hfig);
		zuiformcache(hfig) ;
		zuisavenonok;
		zuireset(h.validation);

	otherwise
		warning('ation non prise en compte')
	   
end

