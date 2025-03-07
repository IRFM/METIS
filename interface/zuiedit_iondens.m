% fichier zuiedit_param_modext.m ->  
%		zuicreeform : creation du formulaire
%		zuiedit_param_modext_fct	: fonction de control des callbacks
%		zuiedit_param_modext_crtl	: fonction de control des donnees

function zuiedit_iondens

[hfig,h] = zuiformhandle('ed_iondens');
if ishandle(hfig)
	zuiformvisible(hfig) ;
	return
end

info=zinfo ;
popupmode =      {'    off mode                    ', ...
                  '    prescribed mode             ', ...
	          '    calculated mode             ', ...
	          '    complex mode                '} ;
valeurmode = {0,1,2,3} ;

% pour mhd.stabilite
popupmode_stab = {'    void mode                    ', ...
                  '    input mode                   ', ...
	          '    calculated mode              ', ...
	          '    complex mode                 ', ...
		  '    post-processing mode        '} ;
valeurmode_stab = {0,1,2,3,4} ;

% pour la mhd
popupmode_m = {' no calculation        ', ...
               ' calculated mode       ', ...
 	            ' complex mode          '} ;
valeurmode_m = {0,1,2} ;

% pour la neo
popupmode_neo = {' input mode             ', ...
                 ' calculated mode        ', ...
                 ' complex mode           '} ;
valeurmode_neo= {1,2,3} ;
 

% pour les glacons
popupmode_g = {' no pellet         ', ...
               ' pellet            ', ...
 	       ' complex mode      '} ;
valeurmode_g = {0,1,2} ;


fonction = evalin('base','param.fonction') ;
%--------------------------------------------------------------------
% le formulaire

nb1   = 10;  % longueur titre module
nb2   = 19;  % longueur nom de la fonction

form = {};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'libel_loadfile','text@full','Plasma composition',[],''};
form{length(form)+1} = {col1};

col1 = {'libel_loadfile','text@full','ions densities',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

%col1 = {'sources','text','',[],''} ;
%col2 = {'sources','text','calculation mode',1,''} ;
%col3 = {'sources','text','module name',[],'',''} ;
%col4 = {'sources','text','change module',[],''};
%col5 = {'sources','text','parameters',[],''} ;
%col6 = {'sources','text','edition mode',[],'',''} ;
%col7 = {'sources','text','reference values',[],'',''} ;
%col8 = {'sources','text','prescribed profiles',[],'',''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

% impur
% -----
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('impur','impurities',popupmode,valeurmode,nb1,nb2,fonction,1);
%col1 = {'impur','text','impurities',10,info.param.fonction.impur,''} ;
%fct = getfield(fonction,'impur') ; if isempty(fct) fct=' '; end
%col2 = {'fct_impur','text',fct,12,info.param.fonction.impur,''} ;
%col3 = {'chng_impur','radio','modify',0,' modify ...'} ;
%col4 = {'para_impur','radio','parameters',0,'parameter function'} ;
%val = evalin('base','data.mode.impur') ;
%ind = fctval(val) + 1 ;
%col5 = {'mode_impur','popup',popupmode,ind,'mode',valeurmode,''} ;
%if ind~=4
%       col6 = {'edit_impur','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%       col6 = {'edit_impur','radio','edition mode',0,'edition mode',[],''} ;
%end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',1,''};
col7 = {'void','frame','void',3,''};
col8 = {'void','frame','void',3,''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

hout=zuicreeform('Edition','ed_iondens','zuiedit_iondens_fct','',form);
