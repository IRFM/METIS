function zuiedit_mhd_calc_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_mhd_calc');

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

% fonctions de MHD
% dds
case {'chng_mhd_dds'}
	liste_fct = liste_module.mhd.dds; %{'zddscrash','other'} ;
	hout = zuichgfct('mhd.dds',liste_fct,h.mhd_dds,h.fct_mhd_dds) ;
	zuireset(hoc) ;

case {'para_mhd_dds'}
	% disp('parametres mhd_dds')
	zuifuninterface('mhd.dds') ;
	zuireset(hoc) ;
		
case {'mode_mhd_dds'}
	% disp('changer mode mhd_dds')
 	value = zuidata(h.mode_mhd_dds) ;
	if value==2
		zuienable(h.edit_mhd_dds) ;
	else
		zuidisable(h.edit_mhd_dds) ;
	end
		
case {'edit_mhd_dds'}
	% disp('edition mode mhd_dds')
	zuiedit_mode('data.mode.mhd.dds',{0,1},{'set to zero','calculated'}) ;
	zuireset(hoc) ;
		
% elm
case {'chng_mhd_elm'}
	% disp('changer elm')
	liste_fct = liste_module.mhd.elm; %{'','other'} ;
	hout = zuichgfct('mhd.elm',liste_fct,h.mhd_elm,h.fct_mhd_elm) ;
	zuireset(hoc) ;

case {'para_mhd_elm'}
	% disp('parametres mhd_elm')
	zuifuninterface('mhd.elm') ;
	zuireset(hoc) ;
		
case {'mode_mhd_elm'}
	% disp('changer mode mhd_elm')
 	value = zuidata(h.mode_mhd_elm) ;
	if value==2
		zuienable(h.edit_mhd_elm) ;
	else
		zuidisable(h.edit_mhd_elm) ;
	end
		
case {'edit_mhd_elm'}
	% disp('edition mode mhd_elm')
	zuiedit_mode('data.mode.mhd.elm',{0,1},{'set to zero','calculated'}) ;
	zuireset(hoc) ;
% limite
case {'chng_mhd_limite'}
	liste_fct = liste_module.mhd.limite; %{'zlimmhd','other'} ;
	hout = zuichgfct('mhd.limite',liste_fct,h.mhd_limite,h.fct_mhd_limite) ;
	zuireset(hoc) ;

case {'para_mhd_limite'}
	% disp('parametres mhd_limite')
	zuifuninterface('mhd.limite') ;
	zuireset(hoc) ;
		
case {'mode_mhd_limite'}
	% disp('changer mode mhd_limite')
 	value = zuidata(h.mode_mhd_limite) ;
	if value==2
		zuienable(h.edit_mhd_limite) ;
	else
		zuidisable(h.edit_mhd_limite) ;
	end
		
case {'edit_mhd_limite'}
	% disp('edition mode limite')
	zuiedit_mode('data.mode.mhd.limite',{0,1},{'set to zero','calculated'}) ;
	zuireset(hoc);

% stabilite
case {'chng_mhd_stab'}
	liste_fct = liste_module.mhd.stab; %{'zlimmhd','other'} ;
        liste_fct
	hout = zuichgfct('mhd.stab',liste_fct,h.mhd_stab,h.fct_mhd_stab) ;
	zuireset(hoc) ;

case {'para_mhd_stab'}
	% disp('parametres mhd_stab')
	zuifuninterface('mhd.stab') ;
	zuireset(hoc) ;
		
case {'mode_mhd_stab'}
	% disp('changer mode mhd_stab')
 	value = zuidata(h.mode_mhd_stab) ;
	if value==3
		zuienable(h.edit_mhd_stab) ;
	else
		zuidisable(h.edit_mhd_stab) ;
	end
		
case {'edit_mhd_stab'}
	% disp('edition mode stab')
	zuiedit_mode('data.mode.mhd.stab',{0,1,2,3,4}, ...
	            {'set to zero','read from data','calculated','copy previous value',...
					'calcule en post-processing'}) ;
	zuireset(hoc);

case {'presprof_mhd_stab'}
        zuiedit_mhdstab_presprof;
        zuireset(hoc);



% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_mhd_calc;
		[hfig,h] = zuiformhandle('ed_mhd_calc');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	zuidisable(h.edit_mhd_dds) ;
	zuidisable(h.edit_mhd_elm) ;
	zuidisable(h.edit_mhd_limite) ;
	zuidisable(h.edit_mhd_stab) ;
	zuireset(h.raz) ;
	
% Annulation
case {'annulation','close'}
		zuicloseone(hfig);	

% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig);

        var = evalin('base','data.mode.mhd.dds') ;
        val = zuidata(h.mode_mhd_dds) ;
        if val~= 3
                zassignin('base','data.mode.mhd.dds',ones(size(var))*val) ;
        end

        var = evalin('base','data.mode.mhd.elm') ;
        val = zuidata(h.mode_mhd_elm) ;
        if val~= 3
                zassignin('base','data.mode.mhd.elm',ones(size(var))*val) ;
        end

        var = evalin('base','data.mode.mhd.limite') ;
        val = zuidata(h.mode_mhd_limite) ;
        if val~= 3
                zassignin('base','data.mode.mhd.limite',ones(size(var))*val) ;
        end

        var = evalin('base','data.mode.mhd.stab') ;
        val = zuidata(h.mode_mhd_stab) ;
        if val~= 3
                zassignin('base','data.mode.mhd.stab',ones(size(var))*val) ;
        end

	zuisavenonok;
	zuireset(h.validation);
	
otherwise
	warning('action not taken into account')
	   
end

