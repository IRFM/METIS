function zuiedit_iondens_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_iondens');

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

% impur
case {'chng_impur'}
	% disp('changer impur')
	liste_fct =liste_module.impur; %{'zinebcompo','zitc','zinebcompovb','other'} ;
	hout = zuichgfct('impur',liste_fct,h.impur,h.fct_impur) ;
	zuireset(hoc) ;

case {'para_impur'}
	% disp('changer impur')
	zuifuninterface('impur') ;
	zuireset(hoc) ;
		
case {'mode_impur'}
	% disp('changer mode impur')
 	value = zuidata(h.mode_impur) ;
        zuiedit_source_enable(h,'impur',value,1);
		
case {'edit_impur'}
	% disp('edition mode impur')
	zuiedit_mode('data.mode.impur',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

case {'refval_impur'}
        zuiedit_data_cons_inj;
        zuireset(hoc) ;

case {'presprof_impur'}
        zuiedit_iondens_presprof;
        zuireset(hoc) ;

% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_iondens;
		[hfig,h] = zuiformhandle('ed_iodens');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	zuidisable(h.edit_impur) ;
	zuireset(h.raz) ;
	
% Annulation
case {'annulation','close'}
		zuicloseone(hfig);	

% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig);
	var = evalin('base','data.mode.impur') ;
	val = zuidata(h.mode_impur) ;
	if val~= 3
		zassignin('base','data.mode.impur',ones(size(var))*val) ;
	end			
		
	zuisavenonok;
	zuireset(h.validation);
	
otherwise
	warning('action not taken into account')
	   
end

