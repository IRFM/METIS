% resmapling in time
function data_out = resample_time_fbe(time_in,data_in,time_out)

% protection against miss oriented vectors
time_in  = time_in(:);
time_out = time_out(:);

indok = find(isfinite(data_in));
if length(indok) > 1
    data_out = interp1(time_in(indok),data_in(indok),time_out,'linear',NaN);
    ind_before = find(time_out < min(time_in(indok)));
    if ~isempty(ind_before)
        data_out(ind_before) = data_in(indok(1));
    end
    ind_after = find(time_out > max(time_in(indok)));
    if ~isempty(ind_after)
        data_out(ind_after) = data_in(indok(end));
    end
elseif length(indok) == 1
    data_out = ones(size(time_out),1) * data_in(indok);
else
    data_out = NaN * ones(size(time_out));
end

if any(~isfinite(data_out(:)))
       warning('NaN or Inf in external equilibrium data after resampling')
       dbstack
end
