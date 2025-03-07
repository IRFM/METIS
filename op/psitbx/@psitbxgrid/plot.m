function plot(g,symb,varargin)

% PSITBXGRID/PLOT	PsiTbx-Grid plot

if nargin < 2
 if g.storage(1) == 'G'
  symb = '-r';
 else
  symb = 'or';
 end
end

gg = psitbxg2g(g,'O');
gg.x = gg.x(1:3);
[ndim,knan] = ndims(gg);
if ndim < 2, return, end

x = gg.x;
x(knan) = [];
lab = gg.label;
lab(knan) = [];
if gg.storage(1) == 'G', [x{:}] = ndgrid(x{:}); end

if ndim == 3
 if g.storage(1) == 'G'
  k = 4:ndims(x{1});
  h = plot3(re2(x{1}),re2(x{2}),re2(x{3}),symb,...
   re2(permute(x{1},[2 1 3 k])),re2(permute(x{2},[2 1 3 k])),re2(permute(x{3},[2 1 3 k])),symb,...
   re2(permute(x{1},[3 2 1 k])),re2(permute(x{2},[3 2 1 k])),re2(permute(x{3},[3 2 1 k])),symb);
  if nargin > 3, set(h,varargin{:}), end
 else
  plot3(re2(x{1}),re2(x{2}),re2(x{3}),symb,varargin{:})
 end
 xlabel(gg.label{1}), ylabel(gg.label{2}), zlabel(gg.label{3})
else
 if g.storage(1) == 'G'
  plot(x{1},x{2},symb,squeeze(x{1})',squeeze(x{2})',symb,varargin{:})
 else
  plot(x{1},x{2},symb,varargin{:})
 end
 xlabel(lab{1}), ylabel(lab{2})
end

axis equal, grid on
title(['PsiTbx-Grid object "' inputname(1) '" (' g.type ',' g.storage ')'])

function x = re2(x)
n = size(x);
x = reshape(x,[n(1),prod(n(2:end))]);