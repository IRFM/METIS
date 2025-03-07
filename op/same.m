% nouvelle version de la fonction same
function cr = same(x,y)

cr = (1>2); % creation d'un logical quelque soit la version de matlab
if ~all(size(x) == size(y))
   return;
elseif isstruct(x) & isstruct(y)
      cr = (2>1);
      noms = fieldnames(x);
      for k = 1:length(noms)
            if ~isfield(y,noms{k})
                  return
            end
            cr = same(getfield(x,noms{k}),getfield(y,noms{k}));
            if cr == 0
                  return
            end
      end
      cr = (2>1);
elseif iscell(x) & iscell(y)
      x = x(:);
      y = y(:);
      for k = 1:length(x)
         cr = same(x{k},y{k});
         if cr == 0
               return
         end
      end
      cr = (2>1);
elseif ischar(x) & ischar(y)
      for k = 1:size(x,1)
            if ~strcmp(x(k,:),y(k,:))
                  return
            end
      end
      cr = (2>1);
elseif isnumeric(x) & isnumeric(y)
      if issparse(x)
         x = full(x);
      end
      if issparse(y)
         y = full(y);
      end
      cr = (all(x(:) == y(:)));
elseif isobject(x) & isobject(y)
   if ~isa(x,class(y)) 
      return
   end
   x = struct(x);
   y = struct(y);
   cr = same(x,y)
end







