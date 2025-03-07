% ZUIDEDIT_ASSER_CALC_FCT gestion callbacks du formulaire ZUIDEDIT_ASSER_CALC
%--------------------------------------------------------------
function zuiedit_asser_calc_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_asser_calc');

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

% asser
case {'chng_asser'}
        % disp('changer asser')
        liste_fct = {'zasserip','other'} ;
        hout = zuichgfct('asser',liste_fct,h.asser,h.fct_asser) ;
        zuireset(hoc) ;

case {'para_asser'}
        % disp('changer plot')
        zuifuninterface('asser');
        zuireset(hoc);

case {'mode_asser'}
        % disp('changer mode plot')
        value = zuidata(h.mode_asser) ;
        if value==3
                zuienable(h.edit_asser) ;
        else
                zuidisable(h.edit_asser) ;
        end

case {'edit_asser'}
        % disp('edition mode plot')
        zuiedit_mode('data.mode.asser',{1,2},...
                    {'read from data','calculated'}) ;
        zuireset(hoc);

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
