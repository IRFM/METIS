function zuiedit_source_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_source');

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

% rip
case {'chng_rip'}
	liste_fct =liste_module.rip;  %{'zripple_therm','other'} ;
	hout = zuichgfct('rip',liste_fct,h.rip,h.fct_rip) ;
	zuireset(hoc) ;

case {'para_rip'}
	zuifuninterface('rip') ;
	zuireset(hoc) ;
		
case {'mode_rip'}
 	value = zuidata(h.mode_rip) ;
        zuiedit_source_enable(h,'rip',value,2);
		
case {'edit_rip'}
	zuiedit_mode('data.mode.rip',{0,1,2}, ...
	            {'set to zero','read from data','calculated'}) ;
	zuireset(hoc) ;

case {'refval_rip'}
        zuireset(hoc) ;

case {'presprof_rip'}
        zuiedit_data_srcrip;
        zuireset(hoc) ;
		
% fonctions de sources
% fci
case {'chng_fci'}
	% disp('changer fci')
	liste_fct = liste_module.fci; %{'zfcifile','zfcisimple','other'} ;
	hout = zuichgfct('fci',liste_fct,h.fci,h.fct_fci) ;
	zuireset(hoc) ;

case {'para_fci'}
	% disp('parametres fci')
	zuifuninterface('fci') ;
	zuireset(hoc) ;
		
case {'mode_fci'}
	% disp('changer mode fci')
 	value = zuidata(h.mode_fci) ;
        zuiedit_source_enable(h,'fci',value,1);
	%if value==3
	%	zuienable(h.edit_fci) ;
	%else
	%	zuidisable(h.edit_fci) ;
	%end
        %if (value==3) || (value==1)
        %        zuienable(h.presprof_fci) ;
        %else
	%        zuidisable(h.presprof_fci) ;
        %end

		
case {'edit_fci'}
	% disp('edition mode fci')
	zuiedit_mode('data.mode.fci',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

case {'refval_fci'}
        zuiedit_data_cons_chauf_fci;
        zuireset(hoc) ;

case {'presprof_fci'}
        zuiedit_data_srcfci ;
        zuireset(hoc) ;		
% fce
case {'chng_fce'}
	liste_fct = liste_module.fce; %{'zfcesimple','zremafile','other'} ;
	hout = zuichgfct('fce',liste_fct,h.fce,h.fct_fce) ;
	zuireset(hoc) ;

case {'para_fce'}
	% disp('parametres fce')
	zuifuninterface('fce') ;
	zuireset(hoc) ;
		
case {'mode_fce'}
	% disp('changer mode fce')
 	value = zuidata(h.mode_fce) ;
        zuiedit_source_enable(h,'fce',value,1);
	%if value==3
	%	zuienable(h.edit_fce) ;
	%else
	%	zuidisable(h.edit_fce) ;
	%end
        %if (value==3) || (value==1)
        %        zuienable(h.presprof_fce) ;
        %else
        %        zuidisable(h.presprof_fce) ;
        %end

		
case {'edit_fce'}
	% disp('edition mode fce')
	zuiedit_mode('data.mode.fce',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

case {'refval_fce'}
        zuiedit_data_cons_chauf_fce;
        zuireset(hoc) ;

case {'presprof_fce'}
        zuiedit_data_srcfce;
        zuireset(hoc) ;
		
% hyb
case {'chng_hyb'}
	liste_fct = liste_module.hyb; %{'zhybsimple','zdelphe','zdke','other'} ;
	hout = zuichgfct('hyb',liste_fct,h.hyb,h.fct_hyb) ;
	zuireset(hoc) ;

case {'para_hyb'}
	% disp('parametres hyb')
	zuifuninterface('hyb') ;
	zuireset(hoc) ;
		
case {'mode_hyb'}
	% disp('changer mode hyb')
 	value = zuidata(h.mode_hyb) ;
        zuiedit_source_enable(h,'hyb',value,1);
		
case {'edit_hyb'}
	% disp('edition mode hyb')
	zuiedit_mode('data.mode.hyb',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

case {'refval_hyb'}
        zuiedit_data_cons_chauf_hyb;
        zuireset(hoc) ;

case {'presprof_hyb'}
        zuiedit_data_srchyb ;
        zuireset(hoc) ;
		
% idn
case {'chng_idn'}
	liste_fct = liste_module.idn; %{'zidnsimple','zsinbad2temps','zsinbadfile','zidnfast','other'} ;
	hout = zuichgfct('idn',liste_fct,h.idn,h.fct_idn) ;
	zuireset(hoc) ;

case {'para_idn'}
	% disp('parametres idn')
	zuifuninterface('idn') ;
	zuireset(hoc) ;
		
case {'mode_idn'}
	% disp('changer mode idn')
 	value = zuidata(h.mode_idn) ;
        zuiedit_source_enable(h,'idn',value,1);
		
case {'edit_idn'}
	% disp('edition mode idn')
	zuiedit_mode('data.mode.idn',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

case {'refval_idn'}
        zuiedit_data_cons_chauf_idn;
        zuireset(hoc) ;

case {'presprof_idn'}
        zuiedit_data_srcidn ;
        zuireset(hoc) ;
		
% n0
case {'chng_n0'}
	liste_fct = liste_module.n0; %{'zneutres','zjonas','other'} ;
	hout = zuichgfct('n0',liste_fct,h.n0,h.fct_n0) ;
	zuireset(hoc) ;

case {'para_n0'}
	% disp('parametres n0')
	zuifuninterface('n0') ;
	zuireset(hoc) ;
		
case {'mode_n0'}
	% disp('changer mode n0')
 	value = zuidata(h.mode_n0) ;
        zuiedit_source_enable(h,'n0',value,2);		

case {'edit_n0'}
	% disp('edition mode n0')
	zuiedit_mode('data.mode.n0',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

case {'presprof_n0'}
        zuiedit_data_srcn0 ;
        zuireset(hoc) ;
		
% bord
case {'chng_bord'}
	liste_fct = liste_module.bord; %{'zrecycle','zmurgaz','other'} ;
	hout = zuichgfct('bord',liste_fct,h.bord,h.fct_bord) ;
	zuireset(hoc) ;

case {'para_bord'}
	% disp('changer bord')
	zuifuninterface('bord') ;
	zuireset(hoc) ;
		
case {'mode_bord'}
	% disp('changer mode bord')
 	value = zuidata(h.mode_bord) ;
        zuiedit_source_enable(h,'bord',value,3);	

case {'edit_bord'}
	% disp('edition mode bord')
	zuiedit_mode('data.mode.bord',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

case {'refval_bord'}
        zuiedit_source_bord;
        zuireset(hoc) ;

% glacon
case {'chng_glacon'}
	liste_fct = liste_module.glacon; %{'zglaquelc','other'} ;
	hout = zuichgfct('glacon',liste_fct,h.glacon,h.fct_glacon) ;
	zuireset(hoc) ;

case {'para_glacon'}
	% disp('changer glacon')
	zuifuninterface('glacon') ;
	zuireset(hoc) ;
		
case {'mode_glacon'}
	% disp('changer mode glacon')
 	value = zuidata(h.mode_glacon) ;
	if value==0
             zuidisable(h.fct_glacon);
             zuidisable(h.chng_glacon);
             zuidisable(h.para_glacon);
             zuidisable(h.edit_glacon);
             zuidisable(h.refval_glacon);
        else
             zuienable(h.fct_glacon);
             zuienable(h.chng_glacon);
             zuienable(h.para_glacon);
             zuienable(h.edit_glacon);
             zuienable(h.refval_glacon);
        end
	
case {'edit_glacon'}
	% disp('edition mode glacon')
	zuiedit_mode('data.mode.glacon',{0,1},...
                    {'off','on'}) ;
	zuireset(hoc) ;

case {'refval_glacon'}
        zuiedit_source_pellet;
        zuireset(hoc) ;

% fusion
case {'chng_fus'}
	% disp('changer fus')
	liste_fct = liste_module.fus; %{'zfusion','zfusion_spot','other'} ;
	hout = zuichgfct('fus',liste_fct,h.fus,h.fct_fus) ;
	zuireset(hoc) ;

case {'para_fus'}
	% disp('changer fus')
	zuifuninterface('fus') ;
	zuireset(hoc) ;
		
case {'mode_fus'}
	% disp('changer mode fus')
 	value = zuidata(h.mode_fus) ;
	zuiedit_source_enable(h,'fus',value,2);
	
case {'edit_fus'}
	% disp('edition mode fus')
	zuiedit_mode('data.mode.fus',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

case {'presprof_fus'}
        zuiedit_data_srcfus ;
        zuireset(hoc) ;

% cyclo
case {'chng_cyclo'}
	% disp('changer cyclo')
	liste_fct =liste_module.cyclo; %{'zalbajar','zcytran','autre'} ;
	hout = zuichgfct('cyclo',liste_fct,h.cyclo,h.fct_cyclo) ;
	zuireset(hoc) ;

case {'para_cyclo'}
	% disp('changer cyclo')
	zuifuninterface('cyclo') ;
	zuireset(hoc) ;
		
case {'mode_cyclo'}
	% disp('changer mode cyclo')
 	value = zuidata(h.mode_cyclo) ;
        zuiedit_source_enable(h,'cyclo',value,4);
		
case {'edit_cyclo'}
	% disp('edition mode cyclo')
	zuiedit_mode('data.mode.cyclo',{0,1,2,3},...
                    {'set to zero','read from data','calculated','k->k+1'}) ;
	zuireset(hoc) ;

case {'presprof_ext'}
        zuiedit_data_srcext ;
        zuireset(hoc) ;

% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_source;
		[hfig,h] = zuiformhandle('ed_source');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	zuidisable(h.edit_rip) ;
	zuidisable(h.edit_fci) ;
	zuidisable(h.edit_fce) ;
	zuidisable(h.edit_hyb) ;
	zuidisable(h.edit_n0) ;
	zuidisable(h.edit_bord) ;
	zuidisable(h.edit_glacon) ;
	zuidisable(h.edit_fus) ;
	zuireset(h.raz) ;
	
% Annulation
case {'annulation','close'}
		zuicloseone(hfig);	

% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig);
		
	var = evalin('base','data.mode.rip') ;
	val = zuidata(h.mode_rip) ;
	if val~= 3
		zassignin('base','data.mode.rip',ones(size(var))*val) ;
	end			

       var = evalin('base','data.mode.fci') ;
	val = zuidata(h.mode_fci) ;
	if val~= 3
		zassignin('base','data.mode.fci',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.fce') ;
	val = zuidata(h.mode_fce) ;
	if val~= 3
		zassignin('base','data.mode.fce',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.hyb') ;
	val = zuidata(h.mode_hyb) ;
	if val~= 3
		zassignin('base','data.mode.hyb',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.idn') ;
	val = zuidata(h.mode_idn) ;
	if val~= 3
		zassignin('base','data.mode.idn',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.n0') ;
	val = zuidata(h.mode_n0) ;
	if val~= 3
		zassignin('base','data.mode.n0',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.bord') ;
	val = zuidata(h.mode_bord) ;
	if val~= 3
		zassignin('base','data.mode.bord',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.glacon') ;
	val = zuidata(h.mode_glacon) ;
	if val <= 1
		zassignin('base','data.mode.glacon',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.fus') ;
	val = zuidata(h.mode_fus) ;
	if val~= 3
		zassignin('base','data.mode.fus',ones(size(var))*val) ;
	end			

	var = evalin('base','data.mode.cyclo') ;
	val = zuidata(h.mode_cyclo) ;
	if val~= 3
		zassignin('base','data.mode.cyclo',ones(size(var))*val) ;
	end			

	zuisavenonok;
	zuireset(h.validation);
	
otherwise
	warning('action not taken into account')
	   
end

