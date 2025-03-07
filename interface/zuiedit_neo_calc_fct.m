function zuiedit_neo_calc_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_neo_calc');

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

% neo
case {'chng_neo'}
       liste_fct = liste_module.neo; %{'zneofile','zneo','other'} ;
       hout = zuichgfct('neo',liste_fct,h.neo,h.fct_neo) ;
       zuireset(hoc) ;

case {'para_neo'}
       zuifuninterface('neo') ;
       zuireset(hoc) ;

case {'mode_neo'}
        value = zuidata(h.mode_neo); 
        if value==1
             zuidisable(h.fct_neo);
             zuidisable(h.chng_neo);
             zuidisable(h.para_neo);
             zuidisable(h.edit_neo);
             zuidisable(h.refval_neo);
             zuidisable(h.presprof_neo);
        else
             zuienable(h.fct_neo);
             zuienable(h.chng_neo);
             zuienable(h.para_neo);
             zuienable(h.edit_neo);
             zuienable(h.refval_neo);
             zuienable(h.presprof_neo);
        end
case {'edit_neo'}
        zuiedit_mode('data.mode.neo',{1,2}, ...
                    {'read from data','calculated'}) ;
        zuireset(hoc) ;

case {'refval_neo'}
        zuiedit_param_gene_funf('zmulti','param.cons.neomulti','Multiplicator');
        zuireset(hoc) ;

case {'presprof_neo'}
        zuiedit_neo_presprof;
        zuireset(hoc) ;

% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_neo_calc;
		[hfig,h] = zuiformhandle('ed_neo_calc');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
        zuidisable(h.edit_neo) ;
	zuireset(h.raz) ;
	
% Annulation
case {'annulation','close'}
		zuicloseone(hfig);	

% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig);

        var = evalin('base','data.mode.neo') ;
        val = zuidata(h.mode_neo) ;
        if val~= 3
                 zassignin('base','data.mode.neo',ones(size(var))*val) ;
        end

	zuisavenonok;
	zuireset(h.validation);
	
otherwise
	warning('action not taken into account')
	   
end

