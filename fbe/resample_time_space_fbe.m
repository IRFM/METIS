% resampling in time and space
function data_out = resample_time_space_fbe(time_in,x_in,data_in,time_out,x_out)

persistent F_mem
if nargin == 0
    F_mem = [];
    return
end

% protection against miss oriented vectors
time_in  = time_in(:);
time_out = time_out(:);


if (length(time_in) == length(time_out)) && all(time_in == time_out)
    indok = find(all(isfinite(data_in),2));
    if length(indok) > 1
        data_out = NaN * ones(length(time_out),length(x_out));
        if size(x_in,1) == 1
            x_mat = ones(length(time_in(indok)),1) * x_in;
        else
            x_mat = x_in(indok,:);
        end
        data_out(indok,:) = tsplinet(x_mat,data_in(indok,:),ones(length(time_out(indok)),1) * x_out);
        if min(indok) > 1
            data_out(1:min(indok),:) = ones(length(1:min(indok)),1) * data_out(min(indok),:);
        end
        if max(indok) > length(time_out)
            data_out(max(indok):length(time_out),:) = ones(length(max(indok):length(time_out)),1) * data_out(max(indok),:);
        end
    elseif length(indok) == 1
        data_out = ones(length(time_out),1) * pchip(x_in(indok,:),data_in(indok,:),x_out);
    else
        data_out = NaN * ones(length(time_out),length(x_out));
    end    
else
    indok = find(all(isfinite(data_in),2));
    if length(indok) > 1
        time_mat = time_in(indok) * ones(1,size(data_in,2));
        if size(x_in,1) == 1
            x_mat = ones(length(time_in(indok)),1) * x_in;
        else
            x_mat = x_in(indok,:);
        end
        data_mat = data_in(indok,:);
        if (length(indok) == length(time_in)) && ~isempty(F_mem)
            F = F_mem;
            F.Values = data_mat(:);
        else
            F = scatteredInterpolant(time_mat(:),x_mat(:),data_mat(:),'natural','linear');
            if length(indok) == length(time_in)
                F_mem = F;
            end
        end
        data_out = F(time_out * ones(size(x_out)),ones(size(time_out)) * x_out);
        ind_before = find(time_out < min(time_in(indok)));
        if ~isempty(ind_before)
            data_out(ind_before,:) = F(time_in(indok(1)) .* (ones(length(ind_before),1) * ones(size(x_out))),ones(length(ind_before),1) * x_out);
        end
        ind_after = find(time_out > max(time_in(indok)));
        if ~isempty(ind_after)
            data_out(ind_after,:) = F(time_in(indok(end)) .* (ones(length(ind_after),1) * ones(size(x_out))),ones(length(ind_after),1) * x_out);
        end
    elseif length(indok) == 1
        data_out = ones(length(time_out),1) * pchip(x_in(indok,:),data_in(indok,:),x_out);
    else
        data_out = NaN * ones(length(time_out),length(x_out));
    end
end

if any(~isfinite(data_out(:)))
    warning('NaN or Inf in external equilibrium data after resampling')
    dbstack
end
