% ZUIEDIT_PROFCMPLX formulaire d'edition de profil complexes
%-------------------------------------------------------------------------------------------------
% fichier zuiedit_profcmplx.m.m ->  
%         zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	formulaire d'edition de profil complexes
%
% syntaxe  :
%	hout=zuiedit_profcmplx(nom_var,var_t,var_x)) ;
%
% entree :
%	nom_var : nom de la variable de retour dans le workspace pour les donnes du profil
%	var_t   : nom de la variable de temps
%	var_x   : nom de la variable d'espace
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 10/12/2002.
% 
% liste des modifications : 
%  * 11/10/2001 -> ajout du trace de la donne initiale avant recomposition
%  * 17/10/2002 -> ajout de l'edition de te,ti, pe et pion
%  * 10/12/2002 -> interface anglais
%--------------------------------------------------------------------------------------------------
function hout=zuiedit_profcmplx(nom_var,var_t,var_x)

if nargin <3
	disp(' Problem inside the editor of complex formular (zuiedit_profcmplx)')
	return
end

% cas specifique de te,pe ti et pion
switch nom_var
case 'data.prof.te'
	evalin('base','param.edit.tepe =''te'';');
case 'data.prof.pe'
	evalin('base','param.edit.tepe =''pe'';');
case 'data.prof.ti'
	evalin('base','param.edit.tipion =''ti'';');
case 'data.prof.pion'
	evalin('base','param.edit.tipion =''pion'';');
end	

% nom du tag
[tagc,reste] = strtok(nom_var,'.') ;
while ~isempty(reste)
	[tagc,reste] = strtok(reste,'.') ;
end

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle([tagc '_profcmplx']) ;
if ishandle(hout)
        zuiformvisible(hout) ;
	return
end

% deomposition de la donnee
% --------------------------
y = evalin('base',nom_var) ;
% suppression des NaN et inf
ind    = find(~isfinite(y));
if ~isempty(ind)
	y(ind) = zeros(1,length(ind));
end
[u,k,v] = svd(y,0) ;
k       = diag(k);


fprintf('k1 = %g, k2 = %g, k3 = %g, sum(k(4:end)) = %g\n',k(1),k(2),k(3),sum(k(4:end)));
e       = k.^2;
e       = sqrt(sum(e(1:3))/sum(e));
fprintf('Energy kept = %4g %%\n',e*100);

nom_var = strrep(nom_var,'(','') ;
nom_var = strrep(nom_var,':','') ;
nom_var = strrep(nom_var,',','') ;
nom_var = strrep(nom_var,')','') ;


% mise en forme
% valeur entre -1 et 1
for ik = 1:3
 	up = max(u(:,ik));
	um = min(u(:,ik));
	if (up > 0) & (um >= 0)
		mu  = 1 ./ up;
	elseif (up <= 0) & (um < 0)
		mu  = 1 ./ um;
	elseif abs(um) > up
		mu  = 1 ./ um;
	else
		mu  = 1 ./ up;
	end
 	vp = max(v(:,ik));
	vm = min(v(:,ik));
	if (vp > 0) & (vm >= 0)
		mv  = 1 ./ vp;
	elseif (vp <= 0) & (vm < 0)
		mv  = 1 ./ vm;
	elseif abs(vm) > vp
		mv  = 1 ./ vm;
	else
		mv  = 1 ./ vp;
	end

	k(ik)   = k(ik) ./ mv ./ mu;
	u(:,ik) = u(:,ik) .* mu;
	v(:,ik) = v(:,ik) .* mv;

end

% Creation de la structure
% ------------------------
t = evalin('base',var_t) ;
cmd = ['complexe.' nom_var '.t = t ;'] ;
eval(cmd) ;
x = evalin('base',var_x) ;
cmd = ['complexe.' nom_var '.x = x ;']  ;
eval(cmd) ;
cmd = ['complexe.' nom_var '.k1 = k(1,:) ;'] ;
eval(cmd) ;
cmd = ['complexe.' nom_var '.u1 = u(:,1) ;'] ;
eval(cmd) ;
cmd = ['complexe.' nom_var '.v1 = v(:,1) ;'] ;
eval(cmd) ;
cmd = ['complexe.' nom_var '.k2 = k(2,:) ;'] ;
eval(cmd) ;
cmd = ['complexe.' nom_var '.u2 = u(:,2) ;'] ;
eval(cmd) ;
cmd = ['complexe.' nom_var '.v2 = v(:,2) ;'] ;
eval(cmd) ;
cmd = ['complexe.' nom_var '.k3 = k(3,:) ;'] ;
eval(cmd) ;
cmd = ['complexe.' nom_var '.u3 = u(:,3) ;'] ;
eval(cmd) ;
cmd = ['complexe.' nom_var '.v3 = v(:,3) ;'] ;
eval(cmd) ;

% donn� reconstitu�
yref = y ;
y = u(:,1) * k(1,:) * v(:,1)' + ...
    u(:,2) * k(2,:) * v(:,2)' + ...
    u(:,3) * k(3,:) * v(:,3)' ;
cmd = ['complexe.' nom_var '.y = y ;'] ;
eval(cmd) ;
zassignin('base','complexe',complexe) ;
    

% valeurs par d�aut
k1 = k(1,:) ;
k2 = k(2,:) ;
k3 = k(3,:) ;

% infos pour les tooltips
info = zinfo ;
aide = eval(['info.' nom_var]) ;

% formulaire avec la sous structure from
form={};

% Titre
% -----
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'nom_var','text@full',['Edition of : ',nom_var],[],aide};
form{1} = {col1};

colj = {'jump','jump','     ',[],''};
coljg = {'jump','jump','          ',[],''};

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% d�endance saptiale du mode
% ---------------------------
colv1 = {'listv1_plot','list','',10,''};
colv2 = {'listv2_plot','list','',10,''};
colv3 = {'listv3_plot','list','',10,''};
coly  = {'listy_plot','list','',10,''};
form{length(form)+1} = {colj,colv1,colj,colv2,colj,colv3,coljg,coly,colj} ;

col1 = {'jump_void','jump','texte pour reserver de la place',[],''};
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% amplitude du mode
% -----------------
col1  = {'text_k1'    ,'text' ,'k1' ,5,'amplitude of the first mode'} ;
col1a = {'edit_k1'    ,'edit' ,k1,2,'amplitude of the first mode'} ;
col2  = {'radio1_plus','radio',' + ',1,'allow to take into account mode 2'} ;
col3  = {'text_k2'    ,'text' ,'k2' ,5,'amplitude of the mode 2'} ;
col3a = {'edit_k2'    ,'edit' ,k2,2,'amplitude of the mode 2'} ;
col4  = {'radio2_plus','radio',' + ',1,'allow to take into account mode 3'} ;
col5  = {'text_k3'    ,'text' ,'k3' ,5,'amplitude of the mode 3'} ;
col5a = {'edit_k3'    ,'edit' ,k3,2,'amplitude of the mode 3'} ;
col6  = {'text_egal'  ,'text' ,' = ' ,5 ,' '} ;
col7  = {'radio_zoom','radio@right','zoom',0,'activate/deactivate the zoom '} ;
form{length(form)+1} = {colj,col1,col1a,col2,col3,col3a,col4,col5,col5a,col6,colj,col7} ;


sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% d�endance temporelle du mode
% -----------------------------
colu1 = {'listu1_plot','list','',10,''};
colu2 = {'listu2_plot','list','',10,''};
colu3 = {'listu3_plot','list','',10,''};
colimag = {'listimag_plot','list','',10,''};
form{length(form)+1} = {colj,colu1,colj,colu2,colj,colu3,coljg,colimag,colj} ;

col1 = {'jump_void','jump','texte pour reserver de la place',[],''};
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;
form{length(form)+1} = {colj,col1,colj,col1,colj,col1,coljg,col1,colj} ;

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% contr�e du temps
% -----------------
colj = {'void','jump','void',[],''} ;
col1 = {'jump','jump','pour reserver beaucoup beaucoup et encore beaucoup de place',[],''} ;
form{length(form)+1} = {colj,colj,colj,colj,colj,col1} ;

col1 = {'text_tps','text','time',5,'profile time'} ;
col2 = {'edit_tps','edit',t(1),2,'profile time'} ;
col3 = {'text_tmin','text',t(1),2,'tmin'} ;
col4 = {'slider_tps','slider',[t(1) t(end)],t(1),'profile time'} ;
col5 = {'text_tmax','text',t(end),2,'tmax'} ;
form{length(form)+1} = {colj,col1,col2,colj,col3,col4,colj,col5} ;

sepa ={'separation_comm','frame','',1,''};
form{length(form)+1} = {sepa};

% Bouton Quit
% ----------
comm{1}={'btn_annul','radio@left' ,'Cancel',0,'close the window'};
comm{2}={'btn_valid','radio@right','Ok',0,'close the window and update variable in the workspace'};

% Formulaire
% ----------
hout=zuicreeform('Complex profile edition',[tagc '_profcmplx'],'zuiedit_profcmplx_fct','',form,comm) ;

% --------------------------------------------------------------------------------------------------------
%
% Creation de la Partie graphique
% -------------------------------

[hform,hui] = zuiformhandle([tagc '_profcmplx']) ;

% definition du slider
% --------------------
posc = get(hui.slider_tps,'Position') ;
set(hui.slider_tps,'Position',[posc(1) posc(2) 0.5 posc(4)]) ;
set(hui.slider_tps,'Min',t(1),'Max',t(end)) ;

% repositionnement du text tmax en bout de slider
post =  get(hui.text_tmax,'Position') ;
set(hui.text_tmax,'Position',[posc(1)+0.5 post(2) post(3) post(4)]) ;

% v1
% --
hui.axes_plotv1=zuiplotinM(hui.listv1_plot);
axes(hui.axes_plotv1) ;
plot(x,v(:,1))
xlabel('x')
ylabel('f_1')
set(hui.axes_plotv1,'Tag','v1') ;
set(hui.axes_plotv1,'ButtonDownFcn','zuiedit_profcmplx_fct(''edit_prof'')');

% v2
% --
hui.axes_plotv2=zuiplotinM(hui.listv2_plot);
axes(hui.axes_plotv2) ;
plot(x,v(:,2))
xlabel('x')
ylabel('f_2')
set(hui.axes_plotv2,'Tag','v2') ;
%hui.context_axes = uicontextmenu;
%set(hui.axes_plotv2,'uicontextmenu',hui.context_axes) ;
%hui.itemv2 = uimenu(hui.context_axes, 'Label', 'editer le mode', 'Callback', 'zuiedit_profcmplx_fct(''edit_prof'')');
set(hui.axes_plotv2,'ButtonDownFcn','zuiedit_profcmplx_fct(''edit_prof'')');

% v3
% --
hui.axes_plotv3=zuiplotinM(hui.listv3_plot);
axes(hui.axes_plotv3) ;
plot(x,v(:,3))
xlabel('x')
ylabel('f_3')
set(hui.axes_plotv3,'Tag','v3') ;
set(hui.axes_plotv3,'ButtonDownFcn','zuiedit_profcmplx_fct(''edit_prof'')');

% recomposition
% -------------
hui.axes_ploty=zuiplotinM(hui.listy_plot) ;
axes(hui.axes_ploty) ;
% trace de la donn� initiale
hz =  plot(x,yref(1,:),'r.') ;
set(hz,'tag','yref')
hold on;
hz = plot(x,y(1,:)) ;
set(hz,'tag','y')
title('u1*k1*v1+*u2*k2*v2+u3*k3*v3') ;

set(hui.axes_ploty,'Tag','y') ;
set(hui.axes_ploty,'ButtonDownFcn','zuiedit_profcmplx_fct(''zoom'')');
xlabel('x') ;

% u1
% --
hui.axes_plotu1=zuiplotinM(hui.listu1_plot);
axes(hui.axes_plotu1) ;
plot(t,u(:,1))
xlabel('t')
ylabel('g_1')
set(hui.axes_plotu1,'Tag','u1') ;
set(hui.axes_plotu1,'Xlim',[t(1) t(end)]) ;
set(hui.axes_plotu1,'ButtonDownFcn','zuiedit_profcmplx_fct(''edit_cons'')');

% u2
% --
hui.axes_plotu2=zuiplotinM(hui.listu2_plot);
axes(hui.axes_plotu2) ;
plot(t,u(:,2))
xlabel('t')
ylabel('g_2')
set(hui.axes_plotu2,'Tag','u2') ;
set(hui.axes_plotu2,'Xlim',[t(1) t(end)]) ;
set(hui.axes_plotu2,'ButtonDownFcn','zuiedit_profcmplx_fct(''edit_cons'')');

% u3
% --
hui.axes_plotu3=zuiplotinM(hui.listu3_plot);
axes(hui.axes_plotu3)
plot(t,u(:,3))
xlabel('t')
ylabel('g_3')
set(hui.axes_plotu3,'Tag','u3') ;
set(hui.axes_plotu3,'Xlim',[t(1) t(end)]) ;
set(hui.axes_plotu3,'ButtonDownFcn','zuiedit_profcmplx_fct(''edit_cons'')');

% image
% -----
hui.axes_img=zuiplotinM(hui.listimag_plot);
axes(hui.axes_img) ;
imagesc(x,t,y) 
xlabel('x')
ylabel('t')
set(gca,'ydir','normal')
colormap('default')
colorbar
title('y')

% memorisation des valeurs
setappdata(hform,'nom_var',nom_var) ;
setappdata(hform,'yref',yref) ;
setappdata(hform,'k1def',k1) ;
setappdata(hform,'k2def',k2) ;
setappdata(hform,'k3def',k3) ;

% memorisation des handles
setappdata(hform,'zhandle',hui) 



