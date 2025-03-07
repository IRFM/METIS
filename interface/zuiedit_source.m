% ZUIEDIT_PARAM_MODEXT formulaire 'modules externes'des  parametres  généraux
%--------------------------------------------------------------
% fichier zuiedit_param_modext.m ->  
%		zuicreeform : creation du formulaire
%		zuiedit_param_modext_fct	: fonction de control des callbacks
%		zuiedit_param_modext_crtl	: fonction de control des donnees
%
% fonction Matlab 5 :
%	formulaire 'modules externes' sous 'Edition'
%	fonction de creation de GUI pour le formulaire
%	des modules externes , parametres  généraux
% 	sous le mode edition du formulaire principal
% 
%
% syntaxe  :
%   	zuiedit_parma_modext ;
%
% entrees
%
% sorties :
%
% fonction ecrite par C. Passeron, poste 61 19
% version 3.0, du 10/01/2005.
% 
% liste des modifications : 
%
%   * 12/09/2001 -> changement de la liste des options pour les glacons  (J-F Artaud)
%   * 16/10/2001 -> bug dans la declaration du module equi (.impur a la place de .equi)
%   * 14/03/2002 -> ajout de la stabilite MHD
%   * 19/03/2002 -> correction bug mode mhd
%   * 10/12/2002 -> interface anglais
%   * 03/09/2003 -> ajout gestion module ripple
%   * 11/09/2003 -> correction detection mode complexe sur neo
%   * 10/11/2003 -> ajout des modules "machine" et "post"
%   * 10/01/2005 -> ajout du module cyclo
%--------------------------------------------------------------

function zuiedit_source

[hfig,h] = zuiformhandle('ed_source');
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

col1 = {'libel_loadfile','text@full','Sources',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

%  fonctions de sources
% ---------------------
col1 = {'sources','text','source type',[],''} ;
col2 = {'sources','text','calculation mode',1,''} ;
col3 = {'sources','text','module name',[],'',''} ;
col4 = {'sources','text','change module',[],''};
col5 = {'sources','text','parameters',[],''} ;
col6 = {'sources','text','edition mode',[],'',''} ;
col7 = {'sources','text','reference values',[],'',''} ;
col8 = {'sources','text','prescribed profiles',[],'',''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

% fci
% ---

[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('fci','ICRF',popupmode,valeurmode,nb1,nb2,fonction,1);
%col1 = {'fci','text','ICRF',nb1,info.param.fonction.fci,''} ;
%fct = getfield(fonction,'fci') ; if isempty(fct) fct=' '; end
%val = evalin('base','data.mode.fci') ;
%ind = fctval(val) + 1 ;
%col2 = {'mode_fci','popup',popupmode,ind,'',valeurmode,''} ;
%col3 = {'fct_fci','text',fct,nb2,'',''} ;
%col4 = {'chng_fci','radio','change module',0,''} ;
%col5 = {'para_fci','radio','parameters',0,'parameter function'} ;
%if ind~=4
%	col6 = {'edit_fci','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%	col6 = {'edit_fci','radio','edition mode',0,'edition mode',[]} ;
%end
%col7 = {'refval_fci','radio','reference values',0,''} ;
%if (ind~=4) && (ind~=2)
%        col8 = {'presprof_fci','radio','prescribed profiles',0,'',[],'','Enable','off'} ;
%else
%        col8 = {'presprof_fci','radio','prescribed profiles',0,'',[]} ;
%end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

% fce
% ---
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('fce','ECRF',popupmode,valeurmode,nb1,nb2,fonction,1);
%col1 = {'fce','text','ECRF',nb1,info.param.fonction.fce,''} ;
%fct = getfield(fonction,'fce') ; if isempty(fct) fct=' '; end
%val = evalin('base','data.mode.fce') ;
%ind = fctval(val) + 1 ;
%col2 = {'mode_fce','popup',popupmode,ind,'',valeurmode,''} ;
%col3 = {'fct_fce','text',fct,nb2,'',''} ;
%col4 = {'chng_fce','radio','change module',0,''} ;
%col5 = {'para_fce','radio','parameters',0,'parameter function'} ;
%if ind~=4
%	col6 = {'edit_fce','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%	col6 = {'edit_fce','radio','edition mode',0,'edition mode',[],''} ;
%end
%col7 = {'refval_fce','radio','reference values',0,''} ;
%if (ind~=4) && (ind~=2)
%        col8 = {'presprof_fce','radio','prescribed profiles',0,'',[],'','Enable','off'} ;
%else
%        col8 = {'presprof_fce','radio','prescribed profiles',0,'',[]} ;
%end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

% hyb
% ---
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('hyb','LH',popupmode,valeurmode,nb1,nb2,fonction,1);
%col1 = {'hyb','text','LH',nb1,info.param.fonction.hyb,''} ;
%fct = getfield(fonction,'hyb') ; if isempty(fct) fct=' '; end
%val = evalin('base','data.mode.hyb') ;
%ind = fctval(val) + 1 ;
%col2 = {'mode_hyb','popup',popupmode,ind,'',valeurmode,''} ;
%col3 = {'fct_hyb','text',fct,nb2,'',''} ;
%col4 = {'chng_hyb','radio','change module',0,''} ;
%col5 = {'para_hyb','radio','parameters',0,'parameter function'} ;
%if ind~=4
%	col6 = {'edit_hyb','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%	col6 = {'edit_hyb','radio','edition mode',0,'edition mode',[],''} ;
%end
%col7 = {'refval_hyb','radio','reference values',0,''} ;
%if (ind~=4) && (ind~=2)
%        col8 = {'presprof_hyb','radio','prescribed profiles',0,'',[],'','Enable','off'} ;
%else
%        col8 = {'presprof_hyb','radio','prescribed profiles',0,'',[]} ;
%end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;


% idn
% ---
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('idn','NBI',popupmode,valeurmode,nb1,nb2,fonction,1);
%col1 = {'idn','text','NBI',nb1,info.param.fonction.idn,''} ;
%fct = getfield(fonction,'idn') ; if isempty(fct) fct=' '; end
%val = evalin('base','data.mode.idn') ;
%ind = fctval(val) + 1 ;
%col2 = {'mode_idn','popup',popupmode,ind,'',valeurmode,''} ;
%col3 = {'fct_idn','text',fct,nb2,'',''} ;
%col4 = {'chng_idn','radio','change module',0,''} ;
%col5 = {'para_idn','radio','parameters',0,'parameter function'} ;
%if ind~=4
%	col6 = {'edit_idn','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%	col6 = {'edit_idn','radio','edition mode',0,'edition mode',[],''} ;
%end
%col7 = {'refval_idn','radio','reference values',0,''} ;
%if (ind~=4) && (ind~=2)
%        col8 = {'presprof_idn','radio','prescribed profiles',0,'',[],'','Enable','off'} ;
%else
%        col8 = {'presprof_idn','radio','prescribed profiles',0,'',[]} ;
%end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

% n0
% --
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('n0','neutrals',popupmode,valeurmode,nb1,nb2,fonction,2);

%col1 = {'n0','text','neutrals',nb1,info.param.fonction.n0,''} ;
%fct = getfield(fonction,'n0') ; if isempty(fct) fct=' '; end
%val = evalin('base','data.mode.n0') ;
%ind = fctval(val) + 1 ;
%col2 = {'mode_n0','popup',popupmode,ind,'',valeurmode,''} ;
%col3 = {'fct_n0','text',fct,nb2,'',''} ;
%col4 = {'chng_n0','radio','change module',0,''} ;
%col5 = {'para_n0','radio','parameters',0,'parameter function'} ;
%if ind~=4
%	col6 = {'edit_n0','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%	col6 = {'edit_n0','radio','edition mode',0,'edition mode',[],''} ;
%end
%col7 = {'sources','text','',[],'',''} ;
%if (ind~=4) && (ind~=2)
%        col8 = {'presprof_n0','radio','prescribed profiles',0,'',[],'','Enable','off'} ;
%else
%        col8 = {'presprof_n0','radio','prescribed profiles',0,'',[]} ;
%end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;


% rip
% ---
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('rip','ripple',popupmode,valeurmode,nb1,nb2,fonction,2);
%col1 = {'rip','text','ripple',nb1,info.param.fonction.rip,''} ;
%fct = getfield(fonction,'rip') ; if isempty(fct) fct=' '; end
%val = evalin('base','data.mode.rip') ;
%ind = fctval(val) + 1 ;
%col2 = {'mode_rip','popup',popupmode,ind,'',valeurmode,''} ;
%col3 = {'fct_rip','text',fct,nb2,'',''} ;
%col4 = {'chng_rip','radio','change module',0,''} ;
%col5 = {'para_rip','radio','parameter',0,'parameter function'} ;
%if ind~=4
%        col6 = {'edit_rip','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%        col6 = {'edit_rip','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
%end
%col7 = {'sources','text','',[],'',''} ;
%if (ind~=4) && (ind~=2)
%        col8 = {'presprof_rip','radio','prescribed profiles',0,'',[],'','Enable','off'} ;
%else
%        col8 = {'presprof_rip','radio','prescribed profiles',0,'',[]} ;
%end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;


% bord
% ----
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('bord','edge',popupmode,valeurmode,nb1,nb2,fonction,3);
%col1 = {'bord','text','edge',10,info.param.fonction.bord,''} ;
%fct = getfield(fonction,'bord') ; if isempty(fct) fct=' '; end
%val = evalin('base','data.mode.bord') ;
%ind = fctval(val) + 1 ;
%col2 = {'mode_bord','popup',popupmode,ind,'',valeurmode,''} ;
%col3 = {'fct_bord','text',fct,12,'',''} ;
%col4 = {'chng_bord','radio','change module',0,''} ;
%col5 = {'para_bord','radio','parameters',0,'parameter function'} ;
%if ind~=4
%	col6 = {'edit_bord','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%	col6 = {'edit_bord','radio','edition mode',0,'edition mode',[],''} ;
%end
%col7 = {'refval_bord','radio','reference values',0,''} ;
%col8 = {'sources','text','',[],'',''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

% glacon
% ------
%[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('glacon','pellet',popupmode,valeurmode,nb1,nb2,fonction,3);
col1 = {'glacon','text','pellet',10,info.param.fonction.glacon,''} ;
fct = getfield(fonction,'glacon') ; if isempty(fct) fct=' '; end
val = evalin('base','data.mode.glacon') ;
ind = fctval(val) + 1 ;
col2 = {'mode_glacon','popup',popupmode_g,ind,'',valeurmode_g,''} ;
col8 = {'sources','text','',[],'',''} ;
if ind==1
      col3 = {'fct_glacon' ,'text' ,fct            ,12,'',[],'','Enable','off'} ;
      col4 = {'chng_glacon','radio','change module',0 ,'',[],'','Enable','off'} ;
      col5 = {'para_glacon','radio','parameters'   ,0 ,'',[],'','Enable','off'} ;
      col6 = {'edit_glacon','radio','edition mode' ,0 ,'',[],'','Enable','off'} ;
      col7 = {'refval_glacon','radio','reference values',0,'',[],'','Enable','off'} ;
else
      col3 = {'fct_glacon','text',fct,12,'',''} ;
      col4 = {'chng_glacon','radio','change module',0,''} ;
      col5 = {'para_glacon','radio','parameters',0,''} ;
      col6 = {'edit_glacon','radio','edition mode',0,''} ;
      col7 = {'refval_glacon','radio','reference values',0,''} ;
end
%if ind~=4
%	col6 = {'edit_glacon','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%	col6 = {'edit_glacon','radio','edition mode',0,'edition mode',[],''} ;
%end
%col7 = {'refval_glacon','radio','reference values',0,''} ;
%col8 = {'sources','text','',[],'',''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

% fus
% ---
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('fus','fusion',popupmode,valeurmode,nb1,nb2,fonction,2);
%col1 = {'fus','text','fusion',10,info.param.fonction.fus,''} ;
%fct = getfield(fonction,'fus') ; if isempty(fct) fct=' '; end
%val = evalin('base','data.mode.fus') ;
%ind = fctval(val) + 1 ;
%col2 = {'mode_fus','popup',popupmode,ind,'',valeurmode,''} ;
%col3 = {'fct_fus','text',fct,12,'',''} ;
%col4 = {'chng_fus','radio','change module',0,''} ;
%col5 = {'para_fus','radio','parameters',0,'parameter function'} ;
%if ind~=4
%	col6 = {'edit_fus','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%	col6 = {'edit_fus','radio','edition mode',0,'edition mode',[],''} ;
%end
%col7 = {'sources','text','',[],'',''} ;
%if (ind~=4) && (ind~=2)
%        col8 = {'presprof_fus','radio','prescribed profiles',0,'',[],'','Enable','off'} ;
%else
%        col8 = {'presprof_fus','radio','prescribed profiles',0,'',[]} ;
%end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;


% cyclo
% -----
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('cyclo','cyclo',popupmode,valeurmode,nb1,nb2,fonction,4);
%col1 = {'cyclo','text','cyclo',10,info.param.fonction.cyclo,''} ;
%fct = getfield(fonction,'cyclo') ; if isempty(fct) fct=' '; end
%val = evalin('base','data.mode.cyclo') ;
%ind = fctval(val) + 1 ;
%col2 = {'mode_cyclo','popup',popupmode,ind,'',valeurmode,''} ;
%col3 = {'fct_cyclo','text',fct,12,'',''} ;
%col4 = {'chng_cyclo','radio','change module',0,''} ;
%col5 = {'para_cyclo','radio','parameters',0,'function parameters'} ;
%if ind~=4
%	col6 = {'edit_cyclo','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
%else
%	col6 = {'edit_cyclo','radio','edition mode',0,'edition mode',[],''} ;
%end
%col7 = {'sources','text','',[],'',''} ;
%col8 = {'sources','text','',[],'',''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;

% ext
% -----
col1 = {'sources','text','external source',[],'',''} ;
col2 = {'sources','text','',[],'',''} ;
col3 = {'sources','text','',[],'',''} ;
col4 = {'sources','text','',[],'',''} ;
col5 = {'sources','text','',[],'',''} ;
col6 = {'sources','text','',[],'',''} ;
col7 = {'sources','text','',[],'',''} ;
col8 = {'presprof_ext','radio','prescribed profiles',0,'',[]}; 
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

hout=zuicreeform('Edition','ed_source','zuiedit_source_fct','',form);
