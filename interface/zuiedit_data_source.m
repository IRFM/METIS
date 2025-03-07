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

function zuiedit_param_modext

[hfig,h] = zuiformhandle('ed_param_modext');
if ishandle(hfig)
	zuiformvisible(hfig) ;
	return
end

info=zinfo ;
popupmode =      {'    void mode                    ', ...
                  '    input mode                   ', ...
	               '    calculated mode              ', ...
	               '    complex mode                 '} ;
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
col1 = {'libel_loadfile','text@full','Parameter control - External modules',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

%  fonctions d'equilibre
% ----------------------
titre = {'equilibre','text@full','equilibrium',1,info.param.fonction.equi,''} ;
form{length(form)+1} = {titre};

% equi
% ----
col1 = {'equi','text','equilibrium',nb1,info.param.fonction.equi,''} ;
fct = getfield(fonction,'equi') ; if isempty(fct) fct=' '; end
col2 = {'fct_equi','text',fct,nb2,info.param.fonction.equi,''} ;
col3 = {'chng_equi','radio','modify',0,' modify ...'} ;
col4 = {'para_equi','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.equi') ;
ind = fctval(val) +1 ;
col5 = {'mode_equi','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_equi','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_equi','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% neo
% ---
col1 = {'neo','text','neoclassical',nb1,info.param.fonction.neo,''} ;
fct = getfield(fonction,'neo') ; if isempty(fct) fct=' '; end
col2 = {'fct_neo','text',fct,nb2,info.param.fonction.neo,''} ;
col3 = {'chng_neo','radio','modify',0,' modify ...'} ;
col4 = {'para_neo','radio','parameter',0,'parameter function'} ;
val = evalin('base','data.mode.neo') ;
ind = fctval(val) + 1 ;
col5 = {'mode_neo','popup',popupmode_neo,ind-1,'mode',valeurmode_neo,''} ;
if ind~=4
	col6 = {'edit_neo','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_neo','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% rip
% ---
col1 = {'rip','text','ripple',nb1,info.param.fonction.rip,''} ;
fct = getfield(fonction,'rip') ; if isempty(fct) fct=' '; end
col2 = {'fct_rip','text',fct,nb2,info.param.fonction.rip,''} ;
col3 = {'chng_rip','radio','modify',0,' modify ...'} ;
col4 = {'para_rip','radio','parameter',0,'parameter function'} ;
val = evalin('base','data.mode.rip') ;
ind = fctval(val) + 1 ;
col5 = {'mode_rip','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_rip','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_rip','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''} ;
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

%  fonctions de mhd
% -----------------
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'mhd','text@full','MHD',1,'MHD',''};
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% mhd.dds
% -------
col1 = {'mhd_dds','text','sawteeth',nb1,info.param.fonction.mhd.dds,''} ;
fct = getfield(fonction.mhd,'dds') ; if isempty(fct) fct=' '; end
col2 = {'fct_mhd_dds','text',fct,nb2,info.param.fonction.mhd.dds,''} ;
col3 = {'chng_mhd_dds','radio','modify',0,' modify ...'} ;
col4 = {'para_mhd_dds','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.mhd.dds') ;
ind = fctval(val) + 1 ;
col5 = {'mode_mhd_dds','popup',popupmode_m,ind,'mode',valeurmode_m,''} ;
if ind~=3
	col6 = {'edit_mhd_dds','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_mhd_dds','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% mhd.elm
% -------
col1 = {'mhd_elm','text','elm',nb1,info.param.fonction.mhd.elm,''} ;
fct = getfield(fonction.mhd,'elm') ; if isempty(fct) fct=' '; end
col2 = {'fct_mhd_elm','text',fct,nb2,info.param.fonction.mhd.elm,''} ;
col3 = {'chng_mhd_elm','radio','modify',0,' modify ...'} ;
col4 = {'para_mhd_elm','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.mhd.elm') ;
ind = fctval(val) + 1 ;
col5 = {'mode_mhd_elm','popup',popupmode_m,ind,'mode',valeurmode_m,''} ;
if ind~=3
	col6 = {'edit_mhd_elm','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_mhd_elm','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% mhd.limite
% ----------
col1 = {'mhd_limite','text','stability threshold',nb1,info.param.fonction.mhd.limite,''} ;
fct = getfield(fonction.mhd,'limite') ; if isempty(fct) fct=' '; end
col2 = {'fct_mhd_limite','text',fct,nb2,info.param.fonction.mhd.limite,''} ;
col3 = {'chng_mhd_limite','radio','modify',0,' modify ...'} ;
col4 = {'para_mhd_limite','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.mhd.limite') ;
ind = fctval(val) + 1 ;
col5 = {'mode_mhd_limite','popup',popupmode_m,ind,'mode',valeurmode_m,''} ;
if ind~=3
	col6 = {'edit_mhd_limite','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_mhd_limite','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% mhd.stab
% ----------
col1 = {'mhd_stab','text','detailed stability',nb1,info.param.fonction.mhd.stab,''} ;
fct = getfield(fonction.mhd,'stab') ; if isempty(fct) fct=' '; end
col2 = {'fct_mhd_stab','text',fct,nb2,info.param.fonction.mhd.stab,''} ;
col3 = {'chng_mhd_stab','radio','modify',0,' modify ...'} ;
col4 = {'para_mhd_stab','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.mhd.stab') ;
ind = fctval4(val) + 1 ;
col5 = {'mode_mhd_stab','popup',popupmode_stab,ind,'mode',valeurmode_stab,''} ;
if ind~=4
	col6 = {'edit_mhd_stab','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_mhd_stab','radio','edition mode',0,'edition mode',[],'','Enable','on'} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

%  fonctions de sources
% ---------------------
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'sources','text@full','sources',1,'source term functions',''};
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% fci
% ---
col1 = {'fci','text','ICRF',nb1,info.param.fonction.fci,''} ;
fct = getfield(fonction,'fci') ; if isempty(fct) fct=' '; end
col2 = {'fct_fci','text',fct,nb2,info.param.fonction.fci,''} ;
col3 = {'chng_fci','radio','modify',0,' modify ...'} ;
col4 = {'para_fci','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.fci') ;
ind = fctval(val) + 1 ;
col5 = {'mode_fci','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_fci','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_fci','radio','edition mode',0,'edition mode',[]} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% fce
% ---
col1 = {'fce','text','ECRF',nb1,info.param.fonction.fce,''} ;
fct = getfield(fonction,'fce') ; if isempty(fct) fct=' '; end
col2 = {'fct_fce','text',fct,nb2,info.param.fonction.fce,''} ;
col3 = {'chng_fce','radio','modify',0,' modify ...'} ;
col4 = {'para_fce','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.fce') ;
ind = fctval(val) + 1 ;
col5 = {'mode_fce','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_fce','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_fce','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% hyb
% ---
col1 = {'hyb','text','LH',nb1,info.param.fonction.hyb,''} ;
fct = getfield(fonction,'hyb') ; if isempty(fct) fct=' '; end
col2 = {'fct_hyb','text',fct,nb2,info.param.fonction.hyb,''} ;
col3 = {'chng_hyb','radio','modify',0,' modify ...'} ;
col4 = {'para_hyb','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.hyb') ;
ind = fctval(val) + 1 ;
col5 = {'mode_hyb','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_hyb','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_hyb','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% idn
% ---
col1 = {'idn','text','NBI',nb1,info.param.fonction.idn,''} ;
fct = getfield(fonction,'idn') ; if isempty(fct) fct=' '; end
col2 = {'fct_idn','text',fct,nb2,info.param.fonction.idn,''} ;
col3 = {'chng_idn','radio','modify',0,' modify ...'} ;
col4 = {'para_idn','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.idn') ;
ind = fctval(val) + 1 ;
col5 = {'mode_idn','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_idn','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_idn','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% n0
% --
col1 = {'n0','text','neutrals',nb1,info.param.fonction.n0,''} ;
fct = getfield(fonction,'n0') ; if isempty(fct) fct=' '; end
col2 = {'fct_n0','text',fct,nb2,info.param.fonction.n0,''} ;
col3 = {'chng_n0','radio','modify',0,' modify ...'} ;
col4 = {'para_n0','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.n0') ;
ind = fctval(val) + 1 ;
col5 = {'mode_n0','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_n0','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_n0','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% bord
% ----
col1 = {'bord','text','edge',10,info.param.fonction.bord,''} ;
fct = getfield(fonction,'bord') ; if isempty(fct) fct=' '; end
col2 = {'fct_bord','text',fct,12,info.param.fonction.bord,''} ;
col3 = {'chng_bord','radio','modify',0,' modify ...'} ;
col4 = {'para_bord','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.bord') ;
ind = fctval(val) + 1 ;
col5 = {'mode_bord','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_bord','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_bord','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% glacon
% ------
col1 = {'glacon','text','pellet',10,info.param.fonction.glacon,''} ;
fct = getfield(fonction,'glacon') ; if isempty(fct) fct=' '; end
col2 = {'fct_glacon','text',fct,12,info.param.fonction.glacon,''} ;
col3 = {'chng_glacon','radio','modify',0,' modify ...'} ;
col4 = {'para_glacon','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.glacon') ;
ind = fctval(val) + 1 ;
col5 = {'mode_glacon','popup',popupmode_g,ind,'mode',valeurmode_g,''} ;
if ind~=4
	col6 = {'edit_glacon','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_glacon','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% fus
% ---
col1 = {'fus','text','fusion',10,info.param.fonction.fus,''} ;
fct = getfield(fonction,'fus') ; if isempty(fct) fct=' '; end
col2 = {'fct_fus','text',fct,12,info.param.fonction.fus,''} ;
col3 = {'chng_fus','radio','modify',0,' modify ...'} ;
col4 = {'para_fus','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.fus') ;
ind = fctval(val) + 1 ;
col5 = {'mode_fus','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_fus','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_fus','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;


% cyclo
% -----
col1 = {'cyclo','text','cyclo',10,info.param.fonction.cyclo,''} ;
fct = getfield(fonction,'cyclo') ; if isempty(fct) fct=' '; end
col2 = {'fct_cyclo','text',fct,12,info.param.fonction.cyclo,''} ;
col3 = {'chng_cyclo','radio','modify',0,' modify ...'} ;
col4 = {'para_cyclo','radio','parameters',0,'function parameters'} ;
val = evalin('base','data.mode.cyclo') ;
ind = fctval(val) + 1 ;
col5 = {'mode_cyclo','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_cyclo','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_cyclo','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;





col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

%  coefficients de transport
% --------------------------
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'transport','text@full','transport',1,'transport equation coefficient',''};
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% coefa
% -----
col1 = {'coefa','text','coefa',10,info.param.fonction.coefa,''} ;
fct = getfield(fonction,'coefa') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefa','text',fct,12,info.param.fonction.coefa,''} ;
col3 = {'chng_coefa','radio','modify',0,' modify ...'} ;
col4 = {'para_coefa','radio','parameters',0,'parameter function'} ;
%val = evalin('base','data.mode.coefa','[]') ;
%if ~isempty(val)
%	ind = val(1)+1 ;
%	col5 = {'mode_coefa','popup',popupmode,ind,'mode',valeurmode,'data.mode.coefa'} ;
%else 
%	col5 = {'mode_coefa','popup',popupmode,1,'mode',valeurmode,''} ;
%end
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefa','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% coefb
% -----
col1 = {'coefb','text','coefb',10,info.param.fonction.coefb,''} ;
fct = getfield(fonction,'coefb') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefb','text',fct,12,info.param.fonction.coefb,''} ;
col3 = {'chng_coefb','radio','modify',0,' modify ...'} ;
col4 = {'para_coefb','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefb','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% coefc
% -----
col1 = {'coefc','text','coefc',10,info.param.fonction.coefc,''} ;
fct = getfield(fonction,'coefc') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefc','text',fct,12,info.param.fonction.coefc,''} ;
col3 = {'chng_coefc','radio','modify',0,' modify ...'} ;
col4 = {'para_coefc','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefc','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% coefd
% -----
col1 = {'coefd','text','coefd',10,info.param.fonction.coefd,''} ;
fct = getfield(fonction,'coefd') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefd','text',fct,12,info.param.fonction.coefd,''} ;
col3 = {'chng_coefd','radio','modify',0,' modify ...'} ;
col4 = {'para_coefd','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefd','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% coefe
% -----
col1 = {'coefe','text','coefe',10,info.param.fonction.coefe,''} ;
fct = getfield(fonction,'coefe') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefe','text',fct,12,info.param.fonction.coefe,''} ;
col3 = {'chng_coefe','radio','modify',0,' modify ...'} ;
col4 = {'para_coefe','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefe','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% coeff
% -----
col1 = {'coeff','text','coeff',10,info.param.fonction.coeff,''} ;
fct = getfield(fonction,'coeff') ; if isempty(fct) fct=' '; end
col2 = {'fct_coeff','text',fct,12,info.param.fonction.coeff,''} ;
col3 = {'chng_coeff','radio','modify',0,' modify ...'} ;
col4 = {'para_coeff','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coeff','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

%  autres fonctions
% --------------------------
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'autres','text@full','miscellaneous',1,'miscellaneous',''};
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% impur
% -----
col1 = {'impur','text','impurities',10,info.param.fonction.impur,''} ;
fct = getfield(fonction,'impur') ; if isempty(fct) fct=' '; end
col2 = {'fct_impur','text',fct,12,info.param.fonction.impur,''} ;
col3 = {'chng_impur','radio','modify',0,' modify ...'} ;
col4 = {'para_impur','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.impur') ;
ind = fctval(val) + 1 ;
col5 = {'mode_impur','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_impur','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_impur','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% plot
% ----
col1 = {'plot','text','plot',10,info.param.fonction.plot,''} ;
fct = getfield(fonction,'plot') ; if isempty(fct) fct=' '; end
col2 = {'fct_plot','text',fct,12,info.param.fonction.plot,''} ;
col3 = {'chng_plot','radio','modify',0,' modify ...'} ;
col4 = {'para_plot','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.plot') ;
ind = fctval(val) + 1 ;
col5 = {'mode_plot','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_plot','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_plot','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% asser
% ----
col1 = {'asser','text','feedback',10,info.param.fonction.asser,''} ;
fct = getfield(fonction,'asser') ; if isempty(fct) fct=' '; end
col2 = {'fct_asser','text',fct,12,info.param.fonction.asser,''} ;
col3 = {'chng_asser','radio','modify',0,' modify ...'} ;
col4 = {'para_asser','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.asser') ;
ind = fctval(val) + 1 ;
col5 = {'mode_asser','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_asser','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_asser','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% post
% ----
col1 = {'post','text','post',10,info.param.fonction.post,''} ;
fct = getfield(fonction,'post') ; if isempty(fct) fct=' '; end
col2 = {'fct_post','text',fct,12,info.param.fonction.post,''} ;
col3 = {'chng_post','radio','modify',0,' modify ...'} ;
col4 = {'para_post','radio','parameters',0,'edit function parameters'} ;
val = evalin('base','data.mode.post') ;
ind = fctval(val)+1;
col5 = {'mode_post','popup',popupmode_m,min(ind,3),'mode',valeurmode_m,''} ;
if ind~=4
	col6 = {'edit_post','radio','edit  mode',0,'edit compute mode',[],'','Enable','off'} ;
else
	col6 = {'edit_post','radio','edit mode',0,'edit compute mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

% machine
% ----
col1 = {'machine','text','device',10,info.param.fonction.machine,''} ;
fct = getfield(fonction,'machine') ; if isempty(fct) fct=' '; end
col2 = {'fct_machine','text',fct,12,info.param.fonction.machine,''} ;
col3 = {'chng_machine','radio','modify',0,' changer ...'} ;
col4 = {'para_machine','radio','parameters',0,'edit function parameters '} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

hout=zuicreeform('Edition','ed_param_modext','zuiedit_param_modext_fct','',form);

function ind = fctval(val)

if all(~isfinite(val))
        ind = 0 ;
elseif all(val==0)
	ind = 0 ;
elseif all(val==1)
	ind = 1 ;
elseif all(val==2)
	ind = 2 ;
else
	ind = 3 ;
end

function ind = fctval4(val)

if all(~isfinite(val))
        ind = 0 ;
elseif all(val==0)
	ind = 0 ;
elseif all(val==1)
	ind = 1 ;
elseif all(val==2)
	ind = 2 ;
elseif all(val==3)
	ind = 3 ;
else
	ind = 4 ;
end

