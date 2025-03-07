    
% plot de control
h=findobj(0,'type','figure','tag','fp');
if isempty(h)
	h=figure('color',[ 1 1 1],'defaultaxesfontsize',12,'tag','fp');
	hc=uicontrol(h,'style','radio','tag','stop');
else
	figure(h);
end
hc = findobj(h,'tag','stop');
% if nb >30
% 	set(hc,'value',1);
% end

subplot(2,2,1)
ii =pe;
plot(gene.x,fp1(:,ii)-F(:,ii),'r',gene.x,fp2(:,ii)-F(:,ii),'g',gene.x,fp3(:,ii)-F(:,ii),'ob', ...
     gene.x,fpini(:,ii)-F(:,ii),'m',gene.x,fp(:,ii)-F(:,ii),'c');
title('r -> adia, g -> lin, b -> expl., m -> ini, c -> fin')
xlabel(['dt = ',num2str(dt)]);  
ind = 1:iround(gene.x,1);
ymax = max([max(fp1(ind,ii)-F(ind,ii)),max(fp2(ind,ii)-F(ind,ii)), ...
           max(fpini(ind,ii)-F(ind,ii)),max(fp(ind,ii)-F(ind,ii))]);
ymin = min([min(fp1(ind,ii)-F(ind,ii)),min(fp2(ind,ii)-F(ind,ii)), ...
           min(fpini(ind,ii)-F(ind,ii)),min(fp(ind,ii)-F(ind,ii))]);
set(gca,'ylim',[ymin,ymax])

subplot(2,2,2)
ii =pion;
plot(gene.x,fp1(:,ii)-F(:,ii),'r',gene.x,fp2(:,ii)-F(:,ii),'g',gene.x,fp3(:,ii)-F(:,ii),'ob', ...
     gene.x,fpini(:,ii)-F(:,ii),'m',gene.x,fp(:,ii)-F(:,ii),'c');
title('r -> adia, g -> lin, b -> expl., m -> ini, c -> fin')
ind = 1:iround(gene.x,1);
ymax = max([max(fp1(ind,ii)-F(ind,ii)),max(fp2(ind,ii)-F(ind,ii)), ...
           max(fpini(ind,ii)-F(ind,ii)),max(fp(ind,ii)-F(ind,ii))]);
ymin = min([min(fp1(ind,ii)-F(ind,ii)),min(fp2(ind,ii)-F(ind,ii)), ...
           min(fpini(ind,ii)-F(ind,ii)),min(fp(ind,ii)-F(ind,ii))]);
set(gca,'ylim',[ymin,ymax])

subplot(2,2,3)
ii =psi;
plot(gene.x,fp1(:,ii)-F(:,ii),'r',gene.x,fp2(:,ii)-F(:,ii),'g',gene.x,fp3(:,ii)-F(:,ii),'ob', ...
     gene.x,fpini(:,ii)-F(:,ii),'m',gene.x,fp(:,ii)-F(:,ii),'c');
title('r -> adia, g -> lin, b -> expl., m -> ini, c -> fin')
ind = 1:iround(gene.x,1);
ymax = max([max(fp1(ind,ii)-F(ind,ii)),max(fp2(ind,ii)-F(ind,ii)), ...
           max(fpini(ind,ii)-F(ind,ii)),max(fp(ind,ii)-F(ind,ii))]);
ymin = min([min(fp1(ind,ii)-F(ind,ii)),min(fp2(ind,ii)-F(ind,ii)), ...
           min(fpini(ind,ii)-F(ind,ii)),min(fp(ind,ii)-F(ind,ii))]);
           set(gca,'ylim',[ymin,ymax])
           
subplot(2,2,4)
semilogy(eee,'or');

while(get(hc,'value') == 1)
   pause(1)
end
