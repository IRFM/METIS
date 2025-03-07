function zuiedit_graphique_fct(action)

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_graphique');
 
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)

case 'graphique_calc'
	if ishandle(hfig)
                zuiedit_graphique_calc;
		zuireset(h.graphique_calc)
	end

case 'graphique_plot'
        if ishandle(hfig)
                zuiedit_param_plot ;
                zuireset(h.graphique_plot)
        end

case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig);
		zuireset(h.btn_quit);
	end

otherwise

end

