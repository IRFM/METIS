function zuiedit_neo_fct(action)

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_neo');
 
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)

case 'neo_calcmode'
	if ishandle(hfig)
                zuiedit_neo_calc;
		zuireset(h.neo_calcmode)
	end

case 'radio_jboot'
        if ishandle(hfig)
                zuiedit_mode('data.mode.jboot',{1,2},{'read from data','calculated'}) ;
                zuireset(h.radio_jboot)
        end

case 'radio_eta'
        if ishandle(hfig)
                zuiedit_mode('data.mode.eta',{1,2},{'read from data','calculated'}) ;
                zuireset(h.radio_eta)
        end

case 'radio_qei'
        if ishandle(hfig)
                zuiedit_mode('data.mode.qei',{0,1,2},{'set to zero','read from data','calculated'}) ;
                zuireset(h.radio_qei)
        end

case 'radio_qneo'
        if ishandle(hfig)
                zuiedit_mode('data.mode.qneo',{0,1,2},{'set to zero','read from data','calculated'}) ;
                zuireset(h.radio_qneo)
        end

case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig);
		zuireset(h.btn_quit);
	end

otherwise

end

