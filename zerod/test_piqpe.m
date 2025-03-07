% test de l'inflence du piquage de la source sur le piquage de pe
cf
close all
piqq = [0,1,1.5,2,3];
qa   = (3:0.5:7)';
vq   = ones(size(qa));
q0   = 1 .* vq;
piqs = 1 * vq;
alpha   = 0 .* vq;
ape_tab  = NaN .* ones(length(piqs),length(qa));
aitb_tab  = NaN .* ones(length(piqs),length(qa));
for k=1:length(piqq)
       [x,pe,pef,ape,aitb] = piqpe(qa,q0,piqq(k)*vq,alpha,piqs);
       ape_tab(k,:) = ape';
       aitb_tab(k,:) = aitb';
end

figure(56);
plot(qa,ape_tab);
