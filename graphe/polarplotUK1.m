% script de plot de la polarimetrie sur ts
t = data.gene.temps;
nl = post.polar.nl;
nlmes = post.polar.nlmes;
af = post.polar.af.*180./pi;
afmes = post.polar.afmes.*180./pi;
if ~isempty(jeux1)
  af1 = jeux1.post.polar.af.*180./pi;
  afmes1 = jeux1.post.polar.afmes.*180./pi;

end
h1 =findobj(0,'type','figure','tag','polarplotnc1');
if isempty(h1)
     h1=figure('color',[1 1 1],'defaultaxesfontsize',14,'defaultaxesfontname', ...
               'times','defaultlinelinewidth',1,'tag','polarplotnc1');
else
     figure(h1);
end
clf

cc=get(gca,'colororder');
plot(t,-af);
inder=iround(t,50);

hold on
errorbar(30,af(inder,4),0.2,0.2)
ind=max(find(~isnan(af(:,1))));

%text('units','data','position',[t(end)-6 af(ind,1)+0.5],'string','chord 1')
%text('units','data','position',[t(end)-6 af(ind,2)+0.7],'string','chord 2')
%text('units','data','position',[t(end)-6 af(ind,3)+0.5],'string','chord 3')
%text('units','data','position',[t(end)-6 af(ind,4)],'string','chord 4')
%text('units','data','position',[t(end)-6 af(ind,5)],'string','chord 5')
if ~isempty(jeux1)
  plot(t1,af1,'--')	
end
set(gca,'colororder',cc);
hp=plot(t,afmes,'o');
set(hp,'markersize',2)

grid
xlabel('time (s)');
ylabel('Faraday angle (°)');
%set(gca,'xlim',[0 20]);
axis([44 51 -2 2])
% correction de l'ecart d'alignement vertical
ind = find(isfinite( prod(af' .* afmes')'));
daf = mean( af(ind,:) -afmes(ind,:));
if ~isempty(jeux1)
  ind1 = find(isfinite( prod(jeux1.post.polar.af' .* jeux1.post.polar.afmes')'));
  daf1= mean( af1(ind1,:) -afmes1(ind1,:));
end
ind = find(isfinite( prod(nl' .* nlmes')'));
dnl = std(nl(ind,:)-nlmes(ind,:))./mean(nlmes(ind,:));
set(gcf,'ResizeFcn','')
title('#59671')

