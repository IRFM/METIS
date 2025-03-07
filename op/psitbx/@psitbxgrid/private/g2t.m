function [x,storage] = g2t(x,storage,nt);

[x,storage] = g2p(x,storage);
[x,storage] = p2t(x,storage,nt);
