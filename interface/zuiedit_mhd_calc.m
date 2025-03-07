function zuiedit_neo_calc

[hfig,h] = zuiformhandle('ed_mhd_calc');
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
popupmode_neo = {' off mode             ', ...
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

col1 = {'libel_loadfile','text@full','MHD',[],''};
form{length(form)+1} = {col1};

col1 = {'libel_loadfile','text@full','calculation mode',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

%  fonctions de sources
% ---------------------
col1 = {'sources','text','',[],''} ;
col2 = {'sources','text','calculation mode',1,''} ;
col3 = {'sources','text','module name',[],'',''} ;
col4 = {'sources','text','change module',[],''};
col5 = {'sources','text','parameters',[],''} ;
col6 = {'sources','text','edition mode',[],'',''} ;
col7 = {'sources','text','prescribed profiles',[],'',''} ;
col8 = {'sources','text','',[],'',''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% mhd.dds
% -------
col1 = {'mhd_dds','text','sawteeth',nb1,info.param.fonction.mhd.dds,''} ;
fct = getfield(fonction.mhd,'dds') ; if isempty(fct) fct=' '; end
val = evalin('base','data.mode.mhd.dds') ;
ind = fctval(val) + 1 ;
col2 = {'mode_mhd_dds','popup',popupmode_m,ind,'mode',valeurmode_m,''} ;
col3 = {'fct_mhd_dds','text',fct,nb2,info.param.fonction.mhd.dds,''} ;
col4 = {'chng_mhd_dds','radio','change module',0,''} ;
col5 = {'para_mhd_dds','radio','parameters',0,'parameter function'} ;
if ind~=3
	col6 = {'edit_mhd_dds','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_mhd_dds','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
col7 = {'sources','text','',[],'',''} ;
col8 = {'sources','text','',[],'',''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;


% mhd.elm
% -------
col1 = {'mhd_elm','text','elm',nb1,info.param.fonction.mhd.elm,''} ;
fct = getfield(fonction.mhd,'elm') ; if isempty(fct) fct=' '; end
val = evalin('base','data.mode.mhd.elm') ;
ind = fctval(val) + 1 ;
col2 = {'mode_mhd_elm','popup',popupmode_m,ind,'mode',valeurmode_m,''} ;
col3 = {'fct_mhd_elm','text',fct,nb2,info.param.fonction.mhd.elm,''} ;
col4 = {'chng_mhd_elm','radio','change module',0,''} ;
col5 = {'para_mhd_elm','radio','parameters',0,'parameter function'} ;
if ind~=3
	col6 = {'edit_mhd_elm','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_mhd_elm','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
col7 = {'sources','text','',[],'',''} ;
col8 = {'sources','text','',[],'',''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;


% mhd.limite
% ----------
col1 = {'mhd_limite','text','stability threshold',nb1,info.param.fonction.mhd.limite,''} ;
fct = getfield(fonction.mhd,'limite') ; if isempty(fct) fct=' '; end
val = evalin('base','data.mode.mhd.limite') ;
ind = fctval(val) + 1 ;
col2 = {'mode_mhd_limite','popup',popupmode_m,ind,'mode',valeurmode_m,''} ;
col3 = {'fct_mhd_limite','text',fct,nb2,info.param.fonction.mhd.limite,''} ;
col4 = {'chng_mhd_limite','radio','change module',0,''} ;
col5 = {'para_mhd_limite','radio','parameters',0,'parameter function'} ;
if ind~=3
	col6 = {'edit_mhd_limite','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_mhd_limite','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
col7 = {'sources','text','',[],'',''} ;
col8 = {'sources','text','',[],'',''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% mhd.stab
% ----------
col1 = {'mhd_stab','text','detailed stability',nb1,info.param.fonction.mhd.stab,''} ;
fct = getfield(fonction.mhd,'stab') ; if isempty(fct) fct=' '; end
val = evalin('base','data.mode.mhd.stab') ;
ind = fctval4(val) + 1 ;
col2 = {'mode_mhd_stab','popup',popupmode_stab,ind,'mode',valeurmode_stab,''} ;
col4 = {'chng_mhd_stab','radio','change module',0,''} ;
col3 = {'fct_mhd_stab','text',fct,nb2,info.param.fonction.mhd.stab,''} ;
col5 = {'para_mhd_stab','radio','parameters',0,'parameter function'} ;
if ind~=4
	col6 = {'edit_mhd_stab','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_mhd_stab','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
col7 = {'presprof_mhd_stab','radio','prescribed profile',0,''} ;
col8 = {'sources','text','',[],'',''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',1,''};
col7 = {'void','frame','void',3,''};
col8 = {'void','frame','void',3,''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;


hout=zuicreeform('Edition','ed_mhd_calc','zuiedit_mhd_calc_fct','',form);
