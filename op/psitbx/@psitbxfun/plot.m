function plot(f,symb,varargin)

% PSITBXFUN/PLOT

if f.grid.type(1) == 'D' & nargin > 1 & any(symb == '-')
 k = size(f.grid.x{1}); k  = [k(1) prod(k(2:end))];
 plot(reshape(f.grid.x{1},k),reshape(f.x(:,:,:,1),k),symb,varargin{:})
 xlabel(f.grid.label{1})
elseif f.grid.storage(1) == 'G'
 z = f.x;
 x = f.grid.x;
 lab = f.grid.label;
 x = x(1:3);
 for k = 3:-1:1
  if length(x{k}) == 1, x(k) = []; lab(k) = ''; end
 end
 z = squeeze(z);
 if ~isempty(f.t), x = {x{:},f.t}; end
 lab = {lab{:},'t'};
 z = z(:,:,:,1);
 x(4:end) = [];
 [x{:}] = ndgrid(x{:});
 switch length(x)
 case 3
  for k = 1:ceil(size(x{1},3)/5):size(x{1},3)
   surf(x{1}(:,:,k),x{2}(:,:,k),x{3}(:,:,k),f.x(:,:,k))
   hold on
  end
  hold off
  shading interp
  colorbar
  xlabel(lab{1}), ylabel(lab{2}), zlabel(lab{3})
  switch [lab{1},lab{2},lab{3}]
   case {'XYZ'}, axis equal
  end
  switch [lab{1},lab{2}]
   case {'XY','XZ','YZ','RZ'}
    k = get(gca,'dataAspectRatio');
    k(3) = k(3)/sqrt(k(1)*k(2));
    k(1:2) = 1;
    set(gca,'dataAspectRatio',k)
  end
 case 2
  set(mesh(x{1},x{2},z),'facecolor','interp',varargin{:})
  shading interp
  colorbar
  xlabel(lab{1}),ylabel(lab{2})
  switch [lab{1},lab{2}]
   case {'XY','XZ','YZ','RZ'}, axis equal auto
  end
 case 1
  if nargin < 2, symb = '-or'; end
  plot(x{1},z,symb,varargin{:})
  xlabel(lab{1})
 end
else
 if nargin < 2, symb = '.'; end
 x = psitbxg2g(f.grid,'O');
 x = x.x;
 z = f.x;
 lim = [min(z(:)), max(z(:))];
 z = z - lim(1);
 z = z / max(z(:));
 c = colormap;
 z = ceil((size(c,1)-1)*z + 1);
 for k = 1:size(c,1)
  l = rem(find(z == k)-1,prod(size(x{1})))+1;
  plot3(x{1}(l),x{2}(l),x{3}(l),symb,varargin{:},...
   'MarkerFaceColor',c(k,:),'MarkerEdgeColor',c(k,:))
  hold on
 end
 hold off
 axis equal
 xlabel('X'),ylabel('Y'),zlabel('Z')
 set(colorbar,'ylim',lim);
 set(findobj('tag','TMW_COLORBAR'),'ydata',lim)
end
grid on
title(['PsiTbx-Function object "' inputname(1) '" (' f.grid.type ',' f.grid.storage ')'])

