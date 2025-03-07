function zuiedit_graphique_calc_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_graphique_calc');

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


% plot
case {'chng_plot'}
	% disp('changer plot')
	liste_fct = liste_module.plot; %{'zplot','other'} ;
	hout = zuichgfct('plot',liste_fct,h.plot,h.fct_plot) ;
	zuireset(hoc) ;

case {'para_plot'}
	% disp('changer plot')
	zuifuninterface('plot');
	zuireset(hoc);
		
case {'mode_plot'}
	% disp('changer mode plot')
 	value = zuidata(h.mode_plot) ;
	if value==3
		zuienable(h.edit_plot) ;
	else
		zuidisable(h.edit_plot) ;
	end
		
case {'edit_plot'}
	% disp('edition mode plot')
	zuiedit_mode('data.mode.plot',{0,1},...
                    {'off','on'}) ;
	zuireset(hoc);

% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_graphique_calc;
		[hfig,h] = zuiformhandle('ed_graphique_calc');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	zuidisable(h.edit_plot) ;
	zuireset(h.raz) ;
	
% Annulation
case {'annulation','close'}
		zuicloseone(hfig);	

% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig);
	var = evalin('base','data.mode.plot') ;
	val = zuidata(h.mode_plot) ;
	if val~= 3
		zassignin('base','data.mode.plot',ones(size(var))*val) ;
	end			
	zuisavenonok;
	zuireset(h.validation);
	
otherwise
	warning('action not taken into account')
	   
end

