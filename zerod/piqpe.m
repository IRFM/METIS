% test de l'inflence du piquage de la source sur le piquage de pe
function [x,pe,pef,ape,aitb] = piqpe(qa,q0,piqq,alpha,piqs)

x = linspace(0,1,21);
ve = ones(size(x));
vt = ones(size(qa));
q = z0qp(x,q0,qa);
s = min(2,max(-1,pdederive(x,q,0,2,2,1) ./ q )); 
fs       = exp(-(s - 0.5 - alpha*ve ./ 2) .^ 2 ./ 2); 

Qe = cumtrapz(x,(1-(vt*x) .^ 2) .^ (piqs*ve) .*(vt*x),2);
inte = -Qe ./ q .^ (piqq*ve) ./ fs;
pe = cumtrapz(x(:,end:-1:1),inte(:,end:-1:1),2);
pe = pe(:,end:-1:1);
pe = pe - pe(:,end)*ve;

pemoy = trapz(x,pe .* (vt*x),2) ./ trapz(x,(vt*x),2);
ape   = pe(:,1) ./pemoy - 1;
pef   = ((pemoy .* (1+ape)) * ve) .* (1- (vt*x) .^ 2) .^ (ape*ve);
aitb   = pe(1) ./ pef(1);

figure(57);
subplot(2,1,1)
plot(x,(pe-pef) ./(pe(:,1) * ve))
subplot(2,1,2)
plot(x,(1-(vt*x) .^ 2) .^ (piqs*ve));

hold on
