function [x,data_out] = metis4muscle3_rhotornorm2xli(temps,rho_tor_norm,data_in,coordinates_change,default_value)

% test inputs
if nargin < 5
    default_value = 0;
end

% select time index
ind_time = find(coordinates_change.temps >= temps,1);
if isempty(ind_time)
    ind_time = length(coordinates_change.temps);
end

% keep the number of point from data
x         = coordinates_change.x;
if any(data_in < -9e40)
    % IMAS NaN
    data_in(data_in < -9e40) = NaN;
end
if isempty(data_in)
    data_in = default_value * ones(size(x));
elseif any(~isfinite(data_in)) 
    if sum(double(~isfinite(data_in)))< 3
        data_in(~isfinite(data_in)) = default_value;
    else
        rho_tor_norm = rho_tor_norm(isfinite(data_in));
        data_in      = data_in(isfinite(data_in));
    end
end    
data_out = interp1_imas(rho_tor_norm,data_in,coordinates_change.rho_tor_norm(ind_time,:),'pchip','extrap');



