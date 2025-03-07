    if 1>2
    h=findobj(0,'type','figure','tag','pion');
    if isempty(h)
    	h=figure('tag','pion','color',[1 1 1],'defaultaxesfontsize',12);
        hi=uicontrol(h,'position',[0,0,0.15,0.05],'units','normalized','string', ...
                  '->','style','radio','tag','->');
    else
    	figure(h);
        hi=findobj(h,'type','uicontrol','tag','->');
    end
    subplot(3,2,1)
    plot(gene.x,datakp1.prof.pion,'r',gene.x,last.pion,'r-.',gene.x,datak.prof.pion,'r:', ...
         gene.x,datakp1.prof.pe,'c',gene.x,last.pe,'c-.',gene.x,datak.prof.pe,'c:');
    title(num2str(dt));
    subplot(3,2,2)
    plot(gene.x,fad)
    subplot(3,2,4)
    hist(abs((datakp1.prof.pion-last.pion)./datak.prof.pion),9);
    title(num2str(delta.pion));
    ylabel('ion')
    subplot(3,2,3)
    hist(abs((datakp1.prof.pe-last.pe)./datak.prof.pe),9);
    title(num2str(delta.pe));
    ylabel('el')
    subplot(3,2,5)
    plot(gene.x,datakp1.coef.ee./datak.prof.ne,'c',gene.x,coef_mem.ee./datak.prof.ne,':c', ...
         gene.x,datakp1.neo.coef.ee./datak.prof.ne,'b',gene.x,neo_mem.ee./datak.prof.ne,':b');
    set(gca,'yscale','log');
    title('Chie')
    subplot(3,2,6)
    plot(gene.x,datakp1.coef.ii./datak.prof.ni,'r',gene.x,coef_mem.ii./datak.prof.ni,':r', ...
         gene.x,datakp1.neo.coef.ii./datak.prof.ni,'m',gene.x,neo_mem.ii./datak.prof.ni,':m');
    set(gca,'yscale','log');
    title('Chii')
    %set(hi,'value',1);
    while(get(hi,'value') ==1)
       pause(1)
    end
    end
    
    
    
    
    
% plot de control
if 2>1
h=findobj(0,'type','figure','tag','fp');
if isempty(h)
	h=figure('color',[ 1 1 1],'defaultaxesfontsize',12,'tag','fp');
	hc=uicontrol(h,'style','radio','tag','stop');
else
	figure(h);
end
hc = findobj(h,'tag','stop');
%set(hc,'value',1)
ii =pion;
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

while(get(hc,'value') == 1)
   pause(1)
end
end




% plot de control
h=findobj(0,'type','figure','tag','fp');
if isempty(h)
	h=figure('color',[ 1 1 1],'defaultaxesfontsize',12,'tag','fp');
	hc=uicontrol(h,'style','radio','tag','stop');
else
	figure(h);
end
hc = findobj(h,'tag','stop');
%set(hc,'value',1)
if nb == 0
	subplot(2,1,1);
	hold off
	subplot(2,1,2);
	hold off
	plot(gene.x,datak.coef.ii,'k');
	hold on
end
cl =get(gca,'colororder');
nbk = rem(nb,size(cl,1)-1)+1; 
subplot(2,1,1);
plot(gene.x,datakp1.coef.ii-datak.coef.ii,'color',cl(nbk,:));
hold on
title(['dt = ',num2str(dt)]);  
subplot(2,1,2);
plot(gene.x,datakp1.coef.ii,'color',cl(nbk,:));
title(['nb = ',num2str(nb)]);  
hold on


while(get(hc,'value') == 1)
   pause(1)
end
