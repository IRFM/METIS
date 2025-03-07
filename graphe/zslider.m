% hf  = handle de la figure
% temps   = vecteur temps du slider
% callback = fontion a appeler ma_fonction(..,temps_zslider,...,option_zslider...)

function zslider(hf,temps,liste_opt,opt_def,callback)

if nargin > 1
   % creation
	figure(hf)
	clf
   % creatioon des objets
	set(hf,'resizefcn','zplotprof(''resize'');');
	setappdata(hf,'zslider_callback',callback);
	%
	% creation des uicontrols
	%
   uicontrol(hf,'Style','frame','Units','normalized' , ...
    		 		 'Position',[0 0 0.999 0.05], ...
                'Backgroundcolor',[0.6 0.6 0.6]); 

   hoption = uicontrol(hf,'Style','popup','Units','normalized' , ...
    		 	             'Position',[0 0.003 0.1 0.04], ...
                         'String',liste_opt, ...
                         'Value',opt_def, ...
                         'HorizontalAlignment','left', ...
                         'tag','option', ...
                         'userdata',0, ...
                         'Callback','zslider(''option'');'); 
								  
   htemps = uicontrol(hf,'Style','edit','Units','normalized' , ...
    		 	             'Position',[0.2 0.003 0.1 0.04], ...
                         'String','0', ...
                         'HorizontalAlignment','left', ...
                         'tag','temps', ...
                         'userdata',0, ...
                         'Callback','zslider(''temps'');'); 

   uicontrol(hf,'Style','text','Units','normalized' , ...
    		 		 'Position',[0.1 0.003 0.1 0.04], ...
                'String','t = ','HorizontalAlignment','left'); 

   uicontrol(hf,'Style','text','Units','normalized' , ...
     		 		  'Position',[0.3 0.003 0.05 0.04], ...
                 'String',' s','HorizontalAlignment','left'); 

   hcurs = uicontrol(hf,'Style','slider','Units','normalized' , ...
    		 	            'Position',[0.46 0.003 0.439 0.04], ...
                        'Min',1e201,'Max',-1e201, ...
                        'tag','curseur', ...
                        'Value',0,'Callback','zslider(''curseur'');');

   hmin =  uicontrol(hf,'Style','text','Units','normalized' , ...
     	     		         'Position',[0.36 0.003 0.1 0.04], ...
                        'String',['0',' s'], ...
                        'tag','minimum', ...
                        'HorizontalAlignment','left');

   hmax = uicontrol(hf,'Style','text','Units','normalized' , ...
     		 	           'Position',[0.899 0.003 0.1 0.04], ...
                       'String',['0', ' s'], ...
                       'tag','maximum', ...
                       'HorizontalAlignment','right');
	
	ind=finite(temps);
	tmin=min(min(temps(ind)));
	tmax=max(max(temps(ind)));
	
	set(hcurs,'min',tmin,'value',tmin);
	set(hcurs,'max',tmax);
	set(htemps,'string',num2str(get(hcurs,'value')));
	set(hmin,'string',[num2str(get(hcurs,'min')),' s']);
	set(hmax,'string',[num2str(get(hcurs,'max')),' s']);
	drawnow
   zplotprof('resize');
	zslider('init');
else
   % callback
	action = hf;
	hf     = gcf;
	ho     = gco;
	%
	% recherche des uicontrols 
	%
	hcurs    = findobj(hf,'tag','curseur');
	hoption  = findobj(hf,'tag','option');
	htemps   = findobj(hf,'tag','temps');
	hmin     = findobj(hf,'tag','minimum');
	hmax     = findobj(hf,'tag','maximum');
	if strcmp(action,'temps')
		 %
		 %
		 % mise a jour des champs  de controle
		 %
		 tc=str2num(get(htemps,'string'));
		 if isempty(tc)
		    set(htemps,'string',num2str(get(hcurs,'value')));
		    return
		 elseif ~finite(tc)
		    set(htemps,'string',num2str(get(hcurs,'value')));
		    return
		 end
		 if get(hcurs,'min')>tc
		    tc=get(hcurs,'min');
		    set(htemps,'string',num2str(get(hcurs,'min')));
		 end
		 if get(hcurs,'max')<tc
		    tc=get(hcurs,'max');
		    set(htemps,'string',num2str(get(hcurs,'max')));
		 end
		 set(hcurs,'value',tc);
	 elseif strcmp(action,'curseur')
	    %
	    % mise a jour des champs de la fenetre de control
	    %
	    set(htemps,'string',num2str(get(hcurs,'value')));
	 end
	 %
	 % appel du calback
	 %
	 % callback = fontion a appeler ma_fonction(..,temps_zslider,...,option_zslider...)
    temps_zslider  = get(hcurs,'value');
	 option_zslider = popupstr(hoption);
	 assignin('base','temps_zslider',temps_zslider);
	 assignin('base','option_zslider',option_zslider);
	 callback       = getappdata(hf,'zslider_callback');
    delete(findobj(hf,'type','axes'));
	 evalin('base',callback); 
end
