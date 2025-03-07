function [x,storage] = g2p(x,storage)
if storage(1) == 'G'
 [x{:}] = ndgrid(x{:});
 storage = 'P';
end
