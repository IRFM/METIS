function [hdlg,value] = z0dpatience(mode)

if nargin == 0

    hdlg = findobj(0,'tag','METIS_waitbar');
    if ~isempty(hdlg) && isappdata(0,'METIS_waitbar_handle');
	hdlg = getappdata(0,'METIS_waitbar_handle');
    end

    if ~isempty(hdlg)
	value = max(get(findobj(hdlg,'type','patch'),'xdata'));
        norme = max(get(findobj(hdlg,'type','axes'),'xlim'));
	value = value ./ norme;
    else
	value = [];
    end

elseif ischar(mode)

    hdlg = findobj(0,'tag','METIS_waitbar');
    if ~isempty(hdlg)
        % adaptation matlab2015b
	for k=1:length(hdlg)
	  try
	    delete(hdlg(k));
	  catch
	    try
	      close(hdlg(k));
	    end
	  end
	end
	try
	    close(getappdata(0,'METIS_waitbar_handle'));
	    delete(getappdata(0,'METIS_waitbar_handle'));
	catch
	  % rien
	end
	rmappdata(0,'METIS_waitbar_handle');
    end
    langue      =  lower(getappdata(0,'langue_cronos'));
    
    switch mode
    case 'delete'
	return
    case 'full'

	    switch langue
	    case 'francais'
		message = 'Appel du simulateur Metis';
	    otherwise
		message = 'Call of Metis simulator' ;
	    end
	   scale = [0.25,0];

    case 'fast'

	    switch langue
	    case 'francais'
		message =  'Appel du simulateur Metis, mode rapide';
	    otherwise
		message = 'Call of Metis simulator, fast computing mode';
	    end
	    scale = [1,0];

    case 'fit'

	    switch langue
	    case 'francais'
		message = 'Appel du simulateur Metis, fit de l''efficaciteb LH';
	    otherwise
		message = 'Call of Metis simulator, fit of LH efficiency';
	  end
	scale = [1,0];
     case 'evolution'

	    switch langue
	    case 'francais'
		message = 'Appel du simulateur Metis en mode evolution';
	    otherwise
		message = 'Call of Metis simulator in evolution mode';
	  end
	scale = [1,0];
   end

    hdlg = waitbar(0,message,'CreateCancelBtn','setappdata(0,''STOP_METIS'',1);','userdata',scale);
    set(hdlg,'tag','METIS_waitbar');
    hb = findobj(hdlg,'type','uicontrol','style','push');
    set(hb,'style','radio','value',0);
    setappdata(0,'METIS_waitbar_handle',hdlg);

    value = 0;

else

    hdlg = findobj(0,'tag','METIS_waitbar');
    if ~isempty(hdlg)
        hdlg_mem = hdlg;
	hdlg = hdlg(1);
        scale = get(hdlg,'userdata');
	waitbar(mode.* scale(1) + scale(2),hdlg);
	hb = findobj(hdlg,'type','uicontrol','style','radio');
	if get(hb,'value') == 1
            % destruction des waitbar mal refermees
	    if length(hdlg_mem) > 1 
	      delete(waitbar(1,hdlg_mem(2:end)));
            end
            % mode working point
            if evalin('base','isfield(z0dinput,''working_point'')')
		evalin('base','z0dinput = jeux1.z0dinput;');
            end
            % ended with error
	    error('METIS was stopped at the request of the user'); 
	end 
    end

end
