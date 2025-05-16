function  ps = make_time_pulse_schedule(ps)


% original time
time = ps.time;

% loop on field
noms = fieldnames(ps);
for k=1:length(noms)
   switch noms{k}
       case {'ids_properties','code'}
            %skip
       otherwise
            time = get_ps_time(ps.(noms{k}),time);
   end
end
time = sort(unique(time));
% set hogeneous time 
ps.ids_properties.homogeneous_time = 1;
ps.time = time;
% resample references
% loop on field
noms = fieldnames(ps);
for k=1:length(noms)
   switch noms{k}
       case {'ids_properties','code'}
            %skip
       otherwise
            ps.(noms{k}) = resample_ps_time(ps.(noms{k}),time);
   end
end



function time = get_ps_time(data,time)

if iscell(data)
    for l=1:length(data)
        time = get_ps_time(data{l},time);
    end
elseif isstruct(data)
    noms = fieldnames(data);
    for l=1:length(noms)
       switch noms{l}
           case 'reference'
               if ~isempty(data.(noms{l}).time)
                   time = union(time,data.(noms{l}).time);
               end
           otherwise
               time = get_ps_time(data.(noms{l}),time);              
       end
    end
else
    % nothing here
end

function data = resample_ps_time(data,time)

if iscell(data)
    for l=1:length(data)
        data{l} = resample_ps_time(data{l},time);
    end
elseif isstruct(data)
    noms = fieldnames(data);
    for l=1:length(noms)
        switch noms{l}
            case 'reference'
                if ~isempty(data.(noms{l}).data)
                    %resample
                    [t,it] = unique(data.(noms{l}).time);
                    val    = data.(noms{l}).data(it);
                    data.(noms{l}).data = interp1_imas(t,val,time,'linear',NaN);
                    indbad = find(~isfinite(data.(noms{l}).data));
                    if ~isempty(indbad)
                        data.(noms{l}).data(indbad) = interp1_imas(t,val,time(indbad),'nearest','extrap');
                    end
                    data.(noms{l}).time = time;
                end
            otherwise
               data.(noms{l}) = resample_ps_time(data.(noms{l}),time);              
       end
    end
else
    % nothing here
end


