function rnh_out = isotope_ratio_west(shot,time,picrh)

% default value
rnh_out = 0.07;
if nargin < 2
  time = [];
end
t_igni = tsbase(shot,'rignitron');
if isempty(t_igni)
    error('Ignitron time is mandatory');
end

% try to read nh/nd+nh in IMAS if available
try
  vs = imas_west_get(shot,'spectrometer_visible');
catch
  vs = [];
end

% test presence of data for isotope_ratios
if ~isempty(vs) && (length(vs.channel) >232) 
  rnh_out = zeros(size(time));
  nbch = 0;
  for k=[208,233,176]
      vsc = vs.channel{k}.isotope_ratios;
      if ~isempty(vsc.signal_to_noise)
	    triso = vsc.time;
	    if isempty(time)
	      time = triso - t_igni;
	      rnh_out = zeros(size(time));
	    end
	    sob   = vsc.signal_to_noise;
	    for l=1:length(vsc.isotope)
	      switch vsc.isotope{l}.label
	      case 'D'
		  rnd   = vsc.isotope{l}.density_ratio;
	      case 'H'
		  rnh   = vsc.isotope{l}.density_ratio;
	      otherwise
		error(sprintf('isotope not taken into account %s',vsc.isotope{l}.label));
	      end
	    end
	    indok = find((sob>0) & (rnd > 0) & (rnh > 0) & ((rnd+rnh) > 0.95) & ((rnd +rnh) < 1.05));
	    %figure(k);plot(triso(indok),rnh(indok));
	    if length(indok) >= 3
	      rnh_out = rnh_out + interp1(triso(indok) - t_igni,rnh(indok),time,'nearest',0);
	      nbch = nbch + 1;
	    end
      end
  end
  rnh_out = rnh_out ./ nbch;
  % security if out of range
  rnh_out(rnh_out > 0.95) = 0.07;
  rnh_out(rnh_out < 1e-3) = 0.07;      
  %figure(21);plot(time,rnh_out);
  if all(isfinite(rnh_out))
    if nargin > 2
	picrh = max(eps,picrh);
	rnh_out = min(1,max(0,trapz(time,rnh_out .* picrh) ./ max(eps,trapz(time, picrh))));
    else
	rnh_out = mean(rnh_out);
    end
    fprintf('Using n_H / (n_H + n_D) computed from visible spectrometer measurements (time averaged value = %g)\n',rnh_out);
  else
    % default value
    rnh_out = 0.07;	
  end
end
