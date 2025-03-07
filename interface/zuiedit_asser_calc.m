% ZUIEDIT_ASSER_CALC
%--------------------------------------------------------------
% fichier zuiedit_asser_calc.m ->
%               zuicreeform : creation du formulaire
%               zuiedit_asser_calc_fct  : fonction de control des callbacks
%
% fonction Matlab 5 :
%       formulaire decrivant le calcul des asservissements
%
%
% syntaxe  :
%       zuiedit_asser_calc_modext ;
%
% entrees
%
% sorties :

function zuiedit_asser_assercalc

[hfig,h] = zuiformhandle('ed_asser_calc');
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
col1 = {'libel_loadfile','text@full','calculation mode',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

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

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',1,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;


hout=zuicreeform('Edition','ed_asser_calc','zuiedit_asser_calc_fct','',form);



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

