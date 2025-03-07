% fonction de placement auto de la figure
function zuiplace(h)

% regarde si les fenetres sont empilees
hh   = get(0,'children');
if length(hh) == 1
   return
end

% bocule pour le placement
nb = 0;
prem =0;
while(nb <10)
   nb = nb +1;
   
   p    = mesure(h);
   ind  = find(hh == h);
   hh(ind) =[];
   ex = 0;
   ey = 0;
   for k = 1:length(hh)
       pp = mesure(hh(k));
       [mx,my] = posmatch(p,pp);
       if mx >0.8
           ex =1;
       end
       if my >0.8
           ey =1;
       end
   end

   if (ex == 0) & (ey == 0)
      return
   end

   % recupere la position en unites normalisees


   % change le format
   p(3) = p(1) + p(3);
   p(4) = p(2) + p(4);

   % calcul de la nouvelle position
   if prem  == 0
       n   = ceil(h) / 25;
       if n > 1
           n  = rand(1);
       end 
       prem = 1;
   else
       n = rand(1);
   end 

   dx  = 0.95 * (1 - p(3));
   dy  = 0.95 * (1 - p(4));

   px = dx * n;
   py = dy * n;


   % application
   umem = get(h,'units');
   set(h,'units','normalized');
   p    = get(h,'position');
   if ex == 1
      p(1) = px;
   end
   if ey == 1
      p(2) = py;
   end
   set(h,'position',p);
   set(h,'units',umem);
end

% mesure la fenetre en unite normalisee
function p = mesure(h)
	
umem = get(h,'units');
set(h,'units','normalized');
p    = get(h,'position');
set(h,'units',umem);


function [mx,my] = posmatch(p,pp)

lx     = min(p(3),pp(3));
ly     = min(p(4),pp(4));
xmin   = max(p(1),pp(1));
xmax   = min(p(1)+p(3),pp(1)+pp(3));
ymin   = max(p(2),pp(2));
ymax   = min(p(2)+p(4),pp(2)+pp(4));

mx     = (xmax -xmin) ./ lx;
my     = (ymax -ymin) ./ ly;






