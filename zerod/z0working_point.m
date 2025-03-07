% script to compute working point using reference at given time slice
% test is this is not already a working point computation 
if isfield(z0dinput,'working_point')
    % do not restore initial setup
    reset_setup = false;
    twp = z0dinput.cons.temps(end-1);
    temps = z0dinput.cons.temps;
    % real time base
    tempo_nova = cat(1,twp - 3000,twp - 2000,twp - 1000,twp);
    z0dinput.cons.temps = tempo_nova;
    if isfield(z0dinput.geo,'temps')
        z0dinput.geo.temps = tempo_nova;
    end
    if isfield(z0dinput.exp0d,'temps')
        z0dinput.exp0d.temps = tempo_nova;
    end
    z0dinput.working_point = true;

else
    if isappdata(0,'WORKING_POINT_TIME')
      twp = getappdata(0,'WORKING_POINT_TIME');
    else
      % selection of the time slice
      z0plot_reference;
      title('choose a time')
      drawnow	
      [twp,void] = ginput(1);
      close(gcf);
      drawnow
    end
    % save original setup
    jeux1.z0dinput = z0dinput;
    reset_setup = true;
   
    % creating new setup with 4 time slices
    %
    tempo_nova = twp * ones(4,1);
    % interpolation des donnees des structures
    temps = z0dinput.cons.temps;
    % consignes
    noms = fieldnames(z0dinput.cons);
    for k=1:length(noms)
	    switch noms{k}
	    case 'temps'
		    z0dinput.cons.temps = tempo_nova;
	    otherwise
		    z0dinput.cons.(noms{k}) = interp1(temps,z0dinput.cons.(noms{k}),tempo_nova,'linear','extrap');
	    end
    end
    % geometrie
    noms = fieldnames(z0dinput.geo);
    for k=1:length(noms)
	    switch noms{k}
	    case 'temps'
		    z0dinput.geo.temps = tempo_nova;
	    otherwise
	      if ~isempty(z0dinput.geo.(noms{k}))
		z0dinput.geo.(noms{k}) = interp1(temps,z0dinput.geo.(noms{k}),tempo_nova,'linear','extrap');
	      end
	    end
    end
    % experience
    noms = fieldnames(z0dinput.exp0d);
    for k=1:length(noms)
	    switch noms{k}
	    case 'temps'
		    z0dinput.exp0d.temps = tempo_nova;
	    otherwise
		    if size(z0dinput.exp0d.(noms{k}),1) == length(temps)
			    warning off
			    z0dinput.exp0d.(noms{k}) = interp1(temps,z0dinput.exp0d.(noms{k}),tempo_nova,'nearest','extrap');
			    warning on
		    end
	    end
    end

    % real time base
    tempo_nova = cat(1,twp - 3000,twp - 2000,twp - 1000,twp);
    z0dinput.cons.temps = tempo_nova;
    if isfield(z0dinput.geo,'temps')
	z0dinput.geo.temps = tempo_nova;   
    end
    if isfield(z0dinput.exp0d,'temps')
	z0dinput.exp0d.temps = tempo_nova;   
    end
     z0dinput.working_point = true; 
end



% run simulation        % pour eviter les incoherence sur le mode restart
% information sur les donnees externes
noms = fieldnames(getappdata(0));
text_warn = '';
for k = 1:length(noms)
    if findstr(noms{k},'_EXP')
        fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
        text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
    end
end
if isdeployed && ~isempty(text_warn)
    warndlg(text_warn,'Pay attention: external data used');
end
clear  z0dstruct

% overwrite parameters whith reference parameters if defined
z0dinput.option = z0doverwriteparam(z0dinput.option);
% modification of parameters:
z0dinput.option.dwdt_method    = 'working_point';
z0dinput.option.runaway        = 0;
z0dinput.option.mode_expo_inte = 1;
z0dinput.option.evolution      = 0;
z0dinput.option.berror         = 0;


% recopie machine dans option
zerod_machine_name;
if isempty(strfind(z0dinput.option.machine,'operation point'))
  z0dinput.option.machine = sprintf('%s [operation point]',z0dinput.option.machine);
end
if isempty(strfind(z0dinput.machine,'operation point'))
  z0dinput.machine        = sprintf('%s [operation point]',z0dinput.machine);
end

% computation
[post.zerod,void,post.profil0d] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
% set time base for graph
%dt = 0.5 .* (max(temps) - min(temps));
%t1 = min(temps);
%t2 = mean(temps);
t1 = twp - max(1e-6,post.zerod.taue(end-1));
t2 = twp - max(1e-6,post.zerod.taue(end-1)) ./ 2;
t3 = twp;
t4 = twp + max(1e-6,post.zerod.taue(end-1)) ./ 2;
%t4 = max(temps);
tempo_nova = unique(sort(cat(1,t1,t2,t3,t4)));
if length(tempo_nova) < 4
    tempo_nova = cat(1,min(tempo_nova) - dt,tempo_nova);
end
if length(tempo_nova) < 4
    tempo_nova = cat(1,min(tempo_nova) + dt,tempo_nova);
end
if length(tempo_nova) < 4
    tempo_nova = cat(1,min(tempo_nova) - dt,tempo_nova);
end
if length(tempo_nova) < 4
    tempo_nova = cat(1,min(tempo_nova) + dt,tempo_nova);
end
z0dinput.cons.temps = tempo_nova;
if isfield(z0dinput.geo,'temps')
	z0dinput.geo.temps = tempo_nova;   
end
if isfield(z0dinput.exp0d,'temps')
	z0dinput.exp0d.temps = tempo_nova;   
end
post.zerod.temps    = tempo_nova; 
post.profil0d.temps = tempo_nova; 
post.z0dinput = z0dinput;
z0plotsc;

% restore original set up
if reset_setup
  z0dinput = jeux1.z0dinput;
end

