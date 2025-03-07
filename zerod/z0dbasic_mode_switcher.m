% change GUI from basic mode to expert mode and vis versa.
function z0dbasic_mode_switcher(flag)

tags = commute_tag_liste;
hfig = zuiformhandle('zeroda');
if isempty(hfig)
  if isappdata(0,'GUI_GLOBAL_BASIC_MODE')
      return
  else
        flag = 0;
  end
end
switch flag
case 1
    % switch to basic mode
    for k=1:length(tags)
	hd = findobj(hfig,'tag',tags{k});
	if ~isempty(hd)
	    set(hd,'visible','off');
	end
    end
    setappdata(0,'GUI_GLOBAL_BASIC_MODE',1)
otherwise
    % switch to expert mode
    for k=1:length(tags)
	hd = findobj(hfig,'tag',tags{k});
	if ~isempty(hd)
	    set(hd,'visible','on');
	end
    end
    setappdata(0,'GUI_GLOBAL_BASIC_MODE',0)
end

% close parameter figures
closeandreopen('z0dsepanew2');
closeandreopen('metis4imas');
closeandreopen('metis4itm');
closeandreopen('zerod');
% redraw all
drawnow



function closeandreopen(fun)

% section menu close and reopen new version
switch fun
case 'z0dsepanew2'
    hf = findobj(0,'type','figure','tag',fun);
    if ~isempty(hf)
       close(hf);
       zuicreefunform('z0dsepanew2','sepa_option',1,0,'z0dinput=z0separatrix(z0dinput,sepa_option);');
    end
    return
otherwise
    hf = findobj(0,'type','figure','tag',sprintf('section_form_%s',fun));
    if ~isempty(hf)
       close(hf);
       zuicreefunform(fun,'z0dinput.option',1,[],'',1);
    end
end

% GUI for parameter edition: close and reopen new version
hf = findobj(0,'type','figure','tag',fun);
if isempty(hf)
  return
end
[s,r] = strtok(get(hf,'name'),'#');
if ~isempty(r)
  section = r(2:end-1);
  close(hf);
  zuicreefunform(sprintf('%s#%s',fun,section),'z0dinput.option',1,[],'',1);
else
  close(hf);
end


function tags = commute_tag_liste

tags ={'radio_audit','radio_pdf','radio_flux','radio_make_flux','radio_ftnbi','radio_temps','radio_noise','radio_external', ...
      'radio_luke','radio_remove','radio_external','radio_hyb','radio_evolution','radio_evolution_restart', ...
      'radio_qlkANN_k','radio_coher0d1d','radio_qlkANN_k_wp','radio_divertor','radio_breakdown','radio_2pts','radio_ramp', ...
      'radio_cost','radio_conv','radio_nbijet','radio_lhacc','radio_ddsts','radio_void_8','radio_void_9','radio_void_10',...
      'radio_export_cpos','radio_qlkANN_k','radio_qlkANN_k_wp','radio_qlkz_std_wp', ...
      'radio_coils_currents','radio_plot_coils_currents','radio_eqdsk','radio_equi_gs'};