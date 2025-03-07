function [p,inter,mask,plot_list] = zplateau(x,t,prec,dplat,low,up)
%function [p,inter,mask,plot_list] = plateau(x,t[,prec[,dplat[,low[,up]]]])
%
% Recherche des plateaux de la fonction x(t)
%
% ENTREES
%
%    x(t,v)             matrice de v * valeurs aux temps t
%    t                  vecteur temps regulierement espace
%    prec [-0.05]       sur un plateau x ne doit pas varier de plus de :
%                          +/- prec si prec > 0 (variation absolue admise)
%                          +/- prec * (max(x) + min(x)) / 2 si prec < 0 (variation relative admise)
%    dplat [0.5 s]      la longueur minimale d' un plateau
%    low   []           limite basse au dela de laquelle x est considere
%    up    []           limite haute en deca de laquelle x est considere
%
% SORTIES
%
%    p(n,v)             valeur moyenne de x sur les n plateaux trouves
%    inter(n,2)         limite des n plateaux trouves
%    mask(n,t)          masque associe i.e. find(mask(n,t)) == [inter(n,1):inter(n,2)]
%    plot_list          liste utilisee pour plotter le resultat
%
% USAGE
%
%    [p,t] = tsbase(19542,'spuiss');
%    [xpl,intpl,mask,liste] = plateau(p,t);
%    plot(t,p,'+',liste{:},'r');
%   
% Th ANIEL (48.95)

if nargin < 6

   low   = [];
   
   if nargin < 5

      up = [];
      
      if nargin < 4
      
         dplat = [];
         
         if nargin < 3

            prec = [];         
            
         end
         
      end
            
   end   
   
end

if isempty(dplat)

   dplat = 0.5;
   
end

if isempty(prec)

   prec = -0.05;
   
end
      
if isempty(low)  

   low  = -Inf;
   
end

if isempty(up)

   up   = Inf;
   
end
      
[nt,nv] = size(x);

if nv > 1

   if length(prec) == 1
   
      prec = prec * ones(1,nv);
      
   end
         
   if length(low) == 1
   
      low  = low * ones(1,nv);
      
   end
         
   if length(up) == 1
   
      up   = up * ones(1,nv);
      
   end
         

end

Nplat   = ceil(dplat / (t(2) - t(1)));

mask = ones(size(t));

for iv=1:nv

   ms   = splat(x(:,iv),prec(iv),Nplat,low(iv),up(iv));
   mask = mask & ms;
   
end

im    = find(mask);
ims   = iconsec(im,Nplat);

nims  = length(ims);
p     = zeros(nims,nv);
inter = zeros(nims,2);
mask  = zeros(nims,nt);

for k=1:nims

   for iv=1:nv
   
      p(k,iv)      = mean(x(ims{k},iv));

   end
   
   inter(k,:)      = [min(t(ims{k})) max(t(ims{k}))];
   mask(k,ims{k})  = ones(size(ims{k}));

end

inter_list = [];

for iv=1:nv

   inter_list = [inter_list; inter];
   
end

plot_list = {inter_list' [p(:) p(:)]'};

function ms = splat(x,prec,Nplat,low,up)

mask = (x > low) & (x < up);

if prec > 0

   step = 2 * prec;

end

im     = find(mask);
ims    = {};

while ~isempty(im)

   
   xmn    = min(x(im));
   xmx    = max(x(im));

   if prec > 0
   
      step = 2 * prec;
   
   else
   
      step = abs(prec * (xmn + xmx));
   
   end
   
   if (xmx - xmn) <= step
      
      if max(im) - min(im) > Nplat
         
         ims     = {ims{:} im};
            
      end
            
      break;

   end
   
   nint   = ceil((xmx - xmn) / step);
   
   [n,y]  = hist(x(im),nint);

   i0     = min(find(n == max(n)));
   x0     = y(i0);

   x1     = x0 - step / 2;
   x2     = x0 + step / 2;
    
   m12    = (x(im) >= x1) & (x(im) <= x2);
   i12    = find(m12);

   if isempty(i12)
   
      break;
      
   end
      
   ims_   = iconsec(im(i12),Nplat);
   ims    = {ims{:} ims_{:}};
   im     = im(find(~m12));

end

ms = zeros(size(x,1),1);

for k=1:length(ims)

   ms(ims{k},1) = ones(length(ims{k}),1);

end
function ims = iconsec(im,N)

ims   = {};
trans = [];

if isempty(im)

   return;
   
end
   
[n,m] = size(im);

if n == 1

   if m == 1
   
      ims   = {im};
      
   else
   
      trans = 0;
      
   end
   
else

   if m == 1
   
      im    = im';
      trans = 1;

   else

      disp('iconsec - 1st argument must be a 1-D matrix');
   
   end
   
end
      
dim   = [Inf diff(im) Inf];
front = find(dim > 1);
dfr   = diff(front);
frpl  = find(dfr >= N);

for k=1:length(frpl)

   ifr    = front(frpl(k));
   ims{k} = [im(front(frpl(k)):front(frpl(k) + 1) - 1)];

   if trans == 1
   
      ims{k} = ims{k}';
   
   end
   
end
