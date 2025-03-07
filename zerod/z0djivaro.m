% script de prolongation de la base temps de metis
% dialogue
prompt={'initial time of simulation or refine interval (s)','final time of simulation or refine interval (s)','time step (s) : > 0 => resample, < 0 => refine,  = 0 => change start and end time and if is text, remove interval between tmin and tmax'};
name='METIS time vector edition';
numlines=1;
defaultanswer={sprintf('%g',min(z0dinput.cons.temps)),sprintf('%g',max(z0dinput.cons.temps)),sprintf('%g',mean(diff(z0dinput.cons.temps)))};
answer=inputdlg(prompt,name,numlines,defaultanswer);
%
% nouvelle base temps
%
if length(answer) < 3
 return
end

temps      = z0dinput.cons.temps;
cutting    = 0;
if isempty(str2num(answer{3}))
  % juste cut but no resampling
  tmin = str2num(answer{1});
  tmax = str2num(answer{2});
  tempo_nova = z0dinput.cons.temps((z0dinput.cons.temps < tmin) | (z0dinput.cons.temps > tmax));
  if isempty(tempo_nova)
     disp('No intersecting time; time edition will be aborted');
     return
  end
  % but time between timin and tmax must be removed at the end
  cutting    = 1;
elseif ~isfinite(str2num(answer{3})) || (str2num(answer{3}) == 0)
  % juste cut but no resampling
  tmin = str2num(answer{1});
  tmax = str2num(answer{2});
  tempo_nova = z0dinput.cons.temps((z0dinput.cons.temps >= tmin) & (z0dinput.cons.temps <= tmax));
elseif str2num(answer{3}) > 0
  tempo_nova = (str2num(answer{1}):str2num(answer{3}):str2num(answer{2}))';
else
  ind_out = find((z0dinput.cons.temps < str2num(answer{1})) | (z0dinput.cons.temps > str2num(answer{2})));
  tempo_nova = union(z0dinput.cons.temps(ind_out),(str2num(answer{1}):abs(str2num(answer{3})):str2num(answer{2}))');
  dt_nova = cat(1,Inf,diff(tempo_nova));
  ind_bad = find((dt_nova < (abs(str2num(answer{3})) - sqrt(eps))) & (tempo_nova >= str2num(answer{1})) & (tempo_nova <= str2num(answer{2})));
  tempo_nova(ind_bad) = [];
  dt_nova(ind_bad) = [];

  h = findobj(0,'type','figure','tag','z0djivaro');
  if isempty(h)
	h=figure('tag','z0djivaro');
  else
	figure(h);
  end
  clf
  set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	  'defaultlinelinewidth',1,'color',[1 1 1])

  plot(tempo_nova,dt_nova,'.b',temps,cat(1,Inf,diff(temps)),'or');
  xlabel('time (s)');
  ylabel('dt (s)');
  title('time base vector');
  legend('New','Old');
  drawnow

end
indplus = find(tempo_nova > max(temps));
  
%
% interpolation des donnees des structures
%
% index for cutting
if cutting == 1
    ind_cut = max(1,fix(interp1(z0dinput.cons.temps,1:length(z0dinput.cons.temps),tempo_nova,'nearest','extrap')));
end
% consignes
noms = fieldnames(z0dinput.cons);
for k=1:length(noms)
	switch noms{k}
	case 'temps'
		z0dinput.cons.temps = tempo_nova;
        if cutting == 1
            z0dinput.cons.temps(tempo_nova >= tmax) = tempo_nova(tempo_nova >= tmax) - tmax + tmin;
        end
	otherwise
        if cutting == 1
            z0dinput.cons.(noms{k}) = z0dinput.cons.(noms{k})(ind_cut,:);
        else
            vlimmin = max(eps,0.1 .* min(z0dinput.cons.(noms{k})));
            vout    = z0dinput.cons.(noms{k})(end);
            z0dinput.cons.(noms{k}) = interp1(temps,z0dinput.cons.(noms{k}),tempo_nova,'linear','extrap');
            switch noms{k}
                case 'flux'
                    % rien
                otherwise
                    z0dinput.cons.(noms{k}) = max(z0dinput.cons.(noms{k}),vlimmin);
            end
            z0dinput.cons.(noms{k})(indplus:end) = vout;
        end
	end
end
% geometrie
noms = fieldnames(z0dinput.geo);
for k=1:length(noms)
	switch noms{k}
	case 'temps'
        z0dinput.geo.temps = tempo_nova;
        if cutting == 1
            z0dinput.geo.temps(tempo_nova >= tmax) = tempo_nova(tempo_nova >= tmax) - tmax + tmin;
        end
   otherwise
       if ~isempty(z0dinput.geo.(noms{k}))
           if cutting == 1
               z0dinput.geo.(noms{k}) = z0dinput.geo.(noms{k})(ind_cut,:);
           else              
               vlimmin = max(eps,0.1 .* min(z0dinput.geo.(noms{k})));
               vout    = z0dinput.geo.(noms{k})(end);
               z0dinput.geo.(noms{k}) = interp1(temps,z0dinput.geo.(noms{k}),tempo_nova,'linear','extrap');
               z0dinput.geo.(noms{k}) = max(z0dinput.geo.(noms{k}),vlimmin);
               z0dinput.geo.(noms{k})(indplus:end) = vout;
           end
       end
	end
end
% experience
noms = fieldnames(z0dinput.exp0d);
for k=1:length(noms)
	switch noms{k}
	case 'temps'
        z0dinput.exp0d.temps = tempo_nova;
        if cutting == 1
            z0dinput.exp0d.temps(tempo_nova >= tmax) = tempo_nova(tempo_nova >= tmax) - tmax + tmin;
        end
	otherwise
        if size(z0dinput.exp0d.(noms{k}),1) == length(temps)
            if cutting == 1
                z0dinput.exp0d.(noms{k}) = z0dinput.exp0d.(noms{k})(ind_cut,:);
            else
                warning off
                z0dinput.exp0d.(noms{k}) = interp1(temps,z0dinput.exp0d.(noms{k}),tempo_nova,'nearest','extrap');
                warning on
            end
        end
	end
end
