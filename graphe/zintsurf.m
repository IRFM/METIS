% integralle surfacique (d'une grandeur independante de theta)
%  s = integrale de surface
%  e = valeur a integree
%  x = coordonnees normalisee
%  sp = datak.equi.sp
%  rhomax = datak.equi.rhomax   
function s=zintsurf(e,x,sp,rhomax)   

    s = rhomax .* trapz(x,sp .* e,2);
