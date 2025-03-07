function zuiedit_neo_calc

[hfig,h] = zuiformhandle('ed_neo_calc');
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

col1 = {'libel_loadfile','text@full','Neoclassical',[],''};
form{length(form)+1} = {col1};

col1 = {'libel_loadfile','text@full','calculation mode',[],''};
form{length(form)+1} = {col1};


sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

%  fonctions de sources
% ---------------------
%col1 = {'sources','text','source type',[],''} ;
%col2 = {'sources','text','calculation mode',1,''} ;
%col3 = {'sources','text','module name',[],'',''} ;
%col4 = {'sources','text','',[],''};
%col5 = {'sources','text','',[],''} ;
%col6 = {'sources','text','',[],'',''} ;
%col7 = {'sources','text','',[],'',''} ;
%col8 = {'sources','text','',[],'',''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

% neo
% ---
col1 = {'neo','text','neoclassical',10,info.param.fonction.neo,''} ;
fct = getfield(fonction,'neo') ; if isempty(fct) fct=' '; end
val = evalin('base','data.mode.neo') ;
ind = fctval(val) + 1 ;
col2 = {'mode_neo','popup',popupmode_neo,ind,'',valeurmode_neo,''} ;
if ind==1
      col3 = {'fct_neo' ,'text' ,fct            ,12,'',[],'','Enable','off'} ;
      col4 = {'chng_neo','radio','change module',0 ,'',[],'','Enable','off'} ;
      col5 = {'para_neo','radio','parameters'   ,0 ,'',[],'','Enable','off'} ;
      col6 = {'edit_neo','radio','edition mode' ,0 ,'',[],'','Enable','off'} ;
      col7 = {'refval_neo','radio','reference values',0,'',[],'','Enable','off'} ;
      col8 = {'presprof_neo','radio','prescribed profile',0,'',[],'','Enable','off'} ;
else
      col3 = {'fct_neo','text',fct,12,'',''} ;
      col4 = {'chng_neo','radio','change module',0,''} ;
      col5 = {'para_neo','radio','parameters',0,''} ;
      col6 = {'edit_neo','radio','edition mode',0,''} ;
      col7 = {'refval_neo','radio','reference values',0,''} ;
      col8 = {'presprof_neo','radio','prescribed profile',0,'',[],'','Enable','on'} ;
end

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

hout=zuicreeform('Edition','ed_neo_calc','zuiedit_neo_calc_fct','',form);
