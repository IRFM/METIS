function zuiedit_compo_fct(action)

%disp('callback : ')
%disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_plasma_composition');
 
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)

case 'radio_compo'
	if ishandle(hfig)
                zuiedit_param_comp ;
                zuireset(h.radio_compo);
	end

case 'inject'
        if ishandle(hfig)               
               zuiedit_iondens
               %zuiedit_data_cons_inj ;
               zuireset(h.inject) ;
        end

case 'radio_ae' 
        if ishandle(hfig)               
               zuiedit_mode('data.mode.ae',{1,2},{'read from data','calculated'}) ;
	       zuireset(h.radio_ae) ;
        end

case 'radio_zeff' 
        if ishandle(hfig)               
               zuiedit_mode('data.mode.zeff',{1,2},{'read from data','calculated'}) ;
	       zuireset(h.radio_zeff) ;
        end

case 'radio_prad' 
        if ishandle(hfig)               
               zuiedit_mode('data.mode.prad',{0,1,2},{'null','read from data','calculated'}) ;
	       zuireset(h.radio_prad) ;
        end

case 'radio_brem' 
        if ishandle(hfig)               
               zuiedit_mode('data.mode.brem',{0,1,2},{'null','read from data','calculated'}) ;
	       zuireset(h.radio_brem) ;
        end

case 'radio_scal'
        if ishandle(hfig)
                zupdatecompo(1);
                zuireset(h.radio_scal)
        end


case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig);
		zuireset(h.btn_quit);
	end

otherwise

end

