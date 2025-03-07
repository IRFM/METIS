t = evalin('base','data.gene.temps') ;
% formulaire avec la sous structure from
form={};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

colj = {'void','jump','void',[],''} ;
col1 = {'slide_tps','jump','pour reserver beaucoup beaucoup et encore beaucoup de place',[],''} ;
form{length(form)+1} = {colj,col1} ;
col1 = {'slide_tps','slider',[t(1) t(end)],1,'temps donné pour affichage du profil'} ;
form{length(form)+1} = {colj,col1} ;

hout=zuicreeform('test slider','test_slide','zuiedit_profcmplx_fct','',form,'') ;
