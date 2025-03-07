function zuiedit_calcequi_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_calcequi');

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

% fonctions equilibre
% equil
case {'chng_equi'}
	% disp('changer equi')
	liste_fct = liste_module.equi;    %{'zequi_helena','other'} ;
	hout = zuichgfct('equi',liste_fct,h.equi,h.fct_equi) ;
	zuireset(hoc) ;

case {'para_equi'}
	% disp('parametres equi') ;
	zuifuninterface('equi') ;
	zuireset(hoc) ;
		
case {'mode_equi'}
	% disp('changer mode equi')
 	value = zuidata(h.mode_equi) ;
        zuiedit_source_enable(h,'equi',value,3);	
	
case {'edit_equi'}
	% disp('edition mode equi')
	zuiedit_mode('data.mode.equi',{0,1,2,3}, ...
	            {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

case {'refval_equi'}
        % disp('parametres equi') ;
        zuiedit_data_geom ; 
        zuireset(hoc) ;




% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_calcequi;
		[hfig,h] = zuiformhandle('ed_calcequi');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	zuidisable(h.edit_equi) ;
	zuireset(h.raz) ;
	
% Annulation
case {'annulation','close'}
		zuicloseone(hfig);	

% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig);
	var = evalin('base','data.mode.equi') ;
	val = zuidata(h.mode_equi) ;
	if val~= 3
		zassignin('base','data.mode.equi',ones(size(var))*val) ;
	end			
	zuisavenonok;
	zuireset(h.validation);
	
otherwise
	warning('action not taken into account')
	   
end

