% geom(1)  =  hauteur point X haut (m)
% geom(2)  =  decalage point X haut par rapport  au centre magnetique (m)
% geom(3), geom(4)  =  angle des 2 tangentes entourant le point X (degre) 





geom(1) = 1.687;
geom(2) = 0.466;
geom(3) = 0.0;
geom(4) = 0.0;
geom(5) = 2.001;
geom(6) = 0.568;
geom(7) = 22.46;
geom(8) = 67.92;

%
%  portion d'ellipses partout
%
mode = 1;
if mode ==1 
  lim2 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;

  if geom(8) > lim2

    geom(8) = lim2-0.1;

  end
end
[x1,y1]=separatrice(geom,mode);
%
% a petit rayon, ra grand rayon
%
a=2;
ra = 6.2;
rr = x1*a+ra;
zz = y1*a;

% angle
%alpha  = abs(unwrap(angle((rr-ra)+ i .* (zz))));
alpha = unwrap(angle((rr-ra)+ i .* (zz)));
if alpha(end) < -5

  alpha=alpha-alpha(end);

end
%
% reechantillonage par spline sur nbp point
%
nbp = 201;
teta  = linspace(0,2*pi,nbp);
ind = find(diff(alpha) == 0);
alpha(ind) = [];
rr(ind) = [];
zz(ind) = [];
R     = spline(alpha',rr',teta);
Z     = spline(alpha',zz',teta);



geom(1) = 1.687;
geom(2) = 0.466;
geom(3) = 0.0;
geom(4) = 0.0;
geom(5) = 2.001;
geom(6) = 0.568;
geom(7) = 22.46;
geom(8) = 67.92;
%
% 2 portions d'ellipses au dessus, 1 ellipse + 1 hyperbol en dessous
%
mode = 2;
if mode == 2
  lim2 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;

  if geom(8) < lim2

    geom(8) = lim2+0.1;

  end
end
%
% 2 portions d'ellipses au dessus, 1 ellipse + 1 hyperbol en dessous
%

[x2,y2]=separatrice(geom,mode);


geom(1) = 1.687;
geom(2) = 0.466;
geom(3) = 0.0;
geom(4) = 0.0;
geom(5) = 2.001;
geom(6) = 0.568;
geom(7) = 22.46;
geom(8) = 67.92;
%
% 2 portions d'ellipses au dessus, 1 ellipse + 1 hyperbol en dessous
%


mode = 3;

lim4 = atan(0.5*geom(1)/(1-geom(2)))*180/pi;
if geom(4) < lim4

  geom(4) = lim4+0.1;

end
lim8 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;
if geom(8) > lim8

  geom(8) = lim8-0.1;

end

[x3,y3]=separatrice(geom,mode);

geom(1) = 1.687;
geom(2) = 0.466;
geom(3) = 0.0;
geom(4) = 0.0;
geom(5) = 2.001;
geom(6) = 0.568;
geom(7) = 22.46;
geom(8) = 67.92;

mode = 4;
[x4,y4]=separatrice(geom,mode);

%
% cas TS
%

geom(1) = 1.;
geom(2) = 0.;
geom(3) = 0.0;
geom(4) = 0.0;
geom(5) = 1;
geom(6) = 0.;
geom(7) = 0;
geom(8) = 0;

%
%  portion d'ellipses partout
%
mode = 1;
if mode ==1 
  lim2 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;

  if geom(8) > lim2

    geom(8) = lim2-0.1;

  end
end
[x1,y1]=separatrice(geom,mode);
%
% a petit rayon, ra grand rayon
%
a=0.72;
ra = 2.4;
rr = x1*a+ra;
zz = y1*a;

% angle
alpha  = unwrap(angle((rr-ra)+ i .* (zz)));
if alpha(end) < -5

  alpha=alpha-alpha(end);

end
%
% reechantillonage par spline sur nbp point
%
nbp = 201;
teta  = linspace(0,2*pi,nbp);
ind = find(diff(alpha) == 0);
alpha(ind) = [];
rr(ind) = [];
zz(ind) = [];
R     = spline(alpha',rr',teta);
Z     = spline(alpha',zz',teta);
