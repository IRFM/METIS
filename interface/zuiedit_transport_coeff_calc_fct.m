function zuiedit_transport_coeff_calc_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_transport_coeff_calc');

% information pour l'assistant
zuicr(hfig,action) ;
try
   hoc = getfield(h,action) ;
catch
   hoc =[];
end
  
fonction = evalin('base','param.fonction') ;

% liste_coef
%liste_coef = {'zjetto     ', ...
%              'zbgbs      ', ...
%              'zcoefsimple', ...
%	      'zitg       ',...
%	      'zetg       ',...
%	      'zetg_hoang ',...
%	      'zetg_stable',...
%	      'zweiland   ',...
%	      'zself      ',...
%	      'zrlwcoef   ',...
%	      'other      '} ;
%
%liste_coef =  zlistcoef;
liste_module =zlist_module;
liste_coef =liste_module.coef;
% selon ation
switch lower(action)

% coefficients de transport
% coefa
case {'chng_coefa'}
	hout = zuichgfct('coefa',liste_coef,h.coefa,h.fct_coefa) ;
	zuireset(hoc) ;

case {'para_coefa'}
	% disp('changer coefa')
	zuifuninterface('coefa') ;
	zuireset(hoc) ;
		
case {'mode_coefa'}
	% disp('changer mode coefa')
 	value = zuidata(h.mode_coefa) ;
	if value==3
		zuienable(h.edit_coefa) ;
	else
		zuidisable(h.edit_coefa) ;
	end
		
case {'edit_coefa'}
	% disp('edition mode coefa')
	zuireset(hoc) ;
		
% coefb
case {'chng_coefb'}
	% disp('changer coefb')
	hout = zuichgfct('coefb',liste_coef,h.coefb,h.fct_coefb) ;
	zuireset(hoc) ;

case {'para_coefb'}
	% disp('changer coefb')
	zuifuninterface('coefb') ;
	zuireset(hoc) ;
		
case {'mode_coefb'}
	% disp('changer mode coefb')
 	value = zuidata(h.mode_coefb) ;
	if value==3
		zuienable(h.edit_coefb) ;
	else
		zuidisable(h.edit_coefb) ;
	end
		
case {'edit_coefb'}
	% disp('edition mode coefb')
	zuireset(hoc) ;
		
% coefc
case {'chng_coefc'}
	hout = zuichgfct('coefc',liste_coef,h.coefc,h.fct_coefc) ;
	zuireset(hoc) ;

case {'para_coefc'}
	% disp('changer coefc')
	zuifuninterface('coefc') ;
	zuireset(hoc) ;
		
case {'mode_coefc'}
	% disp('changer mode coefc')
 	value = zuidata(h.mode_coefc) ;
	if value==3
		zuienable(h.edit_coefc) ;
	else
		zuidisable(h.edit_coefc) ;
	end
		
case {'edit_coefc'}
	% disp('edition mode coefc')
	zuireset(hoc) ;
		
% coefd
case {'chng_coefd'}
	hout = zuichgfct('coefd',liste_coef,h.coefd,h.fct_coefd) ;
	zuireset(hoc) ;

case {'para_coefd'}
	% disp('changer coefd')
	zuifuninterface('coefd') ;
	zuireset(hoc) ;
		
case {'mode_coefd'}
	% disp('changer mode coefd')
 	value = zuidata(h.mode_coefd) ;
	if value==3
		zuienable(h.edit_coefd) ;
	else
		zuidisable(h.edit_coefd) ;
	end
		
case {'edit_coefd'}
	% disp('edition mode coefd')
	zuireset(hoc) ;
		
% coefe
case {'chng_coefe'}
	hout = zuichgfct('coefe',liste_coef,h.coefe,h.fct_coefe) ;
	zuireset(hoc) ;

case {'para_coefe'}
	% disp('changer coefe')
	zuifuninterface('coefe') ;
	zuireset(hoc) ;
		
case {'mode_coefe'}
	% disp('changer mode coefe')
 	value = zuidata(h.mode_coefe) ;
	if value==3
		zuienable(h.edit_coefe) ;
	else
		zuidisable(h.edit_coefe) ;
	end
		
case {'edit_coefe'}
	% disp('edition mode coefe')
	zuireset(hoc) ;
		
% coeff
case {'chng_coeff'}
	hout = zuichgfct('coeff',liste_coef,h.coeff,h.fct_coeff) ;
	zuireset(hoc) ;

case {'para_coeff'}
	% disp('changer coeff')
	zuifuninterface('coeff');
	zuireset(hoc);
		
case {'mode_coeff'}
	% disp('changer mode coeff')
 	value = zuidata(h.mode_coeff) ;
	if value==3
		zuienable(h.edit_coeff) ;
	else
		zuidisable(h.edit_coeff) ;
	end
		
case {'edit_coeff'}
	% disp('edition mode coeff')
	zuireset(hoc) ;
		
% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_param_modext;
		[hfig,h] = zuiformhandle('ed_param_modext');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
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

