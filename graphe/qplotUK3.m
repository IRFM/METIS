% script de plot du profil d q
chargedata
if exist('temps','var')
   tempsold =temps;
   temps = input(['liste des temps a visualiser', mat2str(tempsold)]);
   if isempty(temps )
         temps =tempsold;
    end
 else
    temps =[];
    while isempty(temps)
          temps = input('liste des temps a visualiser ? ');
    end
 end
ind =iround(t,temps);
if ~isempty(jeux1)
  ind1 = iround(t1,temps);
end
texteleg ={};
for k=1:length(temps)
    texteleg{k}=sprintf('t = %4g (s)',temps(k));
end
    

h = findobj(0,'type','figure','tag','qplot2');
if isempty(h)
       h=figure('tag','qplot2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','normal','defaultaxesfontname','times', ...
     'defaultlinelinewidth',2,'color',[1 1 1])
cc=get(gca,'colororder');
set(gca,'colororder',cc);
subplot(1,2,1)
plot(x,data.prof.te(ind,:)/1e3,'-');
legend(texteleg);
hold on
set(gca,'colororder',cc);
if ~isempty(t1)
plot(x(1:3:end),jeux1.data.prof.te(ind1,1:3:end)/1e3,'o');
end
hold off
xlabel('x')
ylabel('Te (keV)')
title('a)')

subplot(1,2,2)

plot(x,data.source.hyb.j(ind,:));

if ~isempty(jeux1)
hold on
plot(x,jeux1.data.source.hyb.j(ind1,:),'o');
end

xlabel('x')
ylabel('A*m^-^2')
title('b)')
