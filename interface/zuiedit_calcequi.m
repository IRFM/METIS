function hout=zuiedit_calcequi

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('ed_calcequi') ;
if ishandle(hout)
        zuiformvisible(hout) ;
        return
end

% infos pour les tooltips
info = zinfo ;

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

% formulaire avec la sous structure from
form={};

% Titre
% -----
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text@full','Equilibrium',[],''};
form{length(form)+1} = {col1};

form{1} = {sepa};

col1 = {'libel','text@full','calculation mode',[],''};
form{length(form)+1} = {col1};

colj=  {'jump1','jump','jump',[],''};

% Séparation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% impur
% -----
[col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active('equi','equilibrium',popupmode,valeurmode,nb1,nb2,fonction,3);
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;


% separation
% ----------
col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',1,''};
col7 = {'void','frame','void',3,''};
col8 = {'void','frame','void',3,''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;


% Formulaire
% ----------
hout=zuicreeform('Edition','ed_calcequi','zuiedit_calcequi_fct','',form) ;


