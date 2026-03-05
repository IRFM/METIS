% This file is for the standalone simulation of METIS+TGLF
% set parameters
waitfor(zuicreefunform('gui_tglf','tglfparams',1,0,'','visible'));
if getappdata(0,'ZGUI_CANCEL')
    return;
end

z0plot_reference;
title('choose a time')
drawnow	
[twp,void] = ginput(1);

% Ask user to enter the time manually
prompt = {'Check the time value:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {num2str(twp)};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    return; 	% Exit if the user cancels the input dialog
end

twp = str2double(answer{1});
close(gcf);
drawnow

[tmp, itime] = min(abs(post.profil0d.temps - twp));
currentpath = pwd;

% Run tglf and get chis + some plots
switch tglfparams.SWITCH_MODEL
case 0
  [Chie,Chii,De,Ve] = run_tglf_3sp(tglfparams,post,itime,1);
end

cd(currentpath)
