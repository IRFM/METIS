function zuiedit_graphique_fct(action)

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_postrait');
 
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)

case 'postrait_param'
	if ishandle(hfig)
                zuiedit_param_postrait_param;
		zuireset(h.postrait_param)
	end

case 'postrait_calc'
        if ishandle(hfig)
                zuiedit_param_postrait_calc ;
                zuireset(h.postrait_calc)
        end

case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig);
		zuireset(h.btn_quit);
	end

otherwise

end

