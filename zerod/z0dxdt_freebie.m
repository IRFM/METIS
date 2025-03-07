function dxdt = z0dxdt_freebie(x,t)

dxdt = (x(end,:) - x(1,:)) ./ (t(end) - t(1));
dxdt = ones(size(x)) * dxdt;
