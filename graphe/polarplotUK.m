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
plot(t,af);
inder=iround(t,4);

hold on
errorbar(4,af(inder,1),0.2,0.2)
ind=max(find(~isnan(af(:,1))));

text('units','data','position',[t(end) af(ind,1)],'string','chord 1')
text('units','data','position',[t(end) af(ind,2)-0.5],'string','chord 2')
text('units','data','position',[t(end) af(ind,3)],'string','chord 3')
text('units','data','position',[t(end) af(ind,4)],'string','chord 4')
text('units','data','position',[t(end) af(ind,5)],'string','chord 5')
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

h2 =findobj(0,'type','figure','tag','polarplotnc2');
if isempty(h2)
     h2=figure('color',[1 1 1],'defaultaxesfontsize',12,'defaultaxesfontname', ...
               'times','defaultlinelinewidth',2,'tag','polarplotnc2');
else
     figure(h2);
end
clf
%figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname','times','defaultlinelinewidth',3);
cc=get(gca,'colororder');
plot(t,af-ones(length(t),1)*daf);
hold on
if ~isempty(jeux1)
  plot(t1,af1-ones(length(t1),1)*daf1,'--');
end
set(gca,'colororder',cc);
plot(t,afmes,'o');
st =sprintf('choc %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d (correction offset)',param.from.machine, ...              
            fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
%title(st);
grid
xlabel('time (s)');
ylabel('degrees');
%legend(num2str(daf(1),2),num2str(daf(2),2),num2str(daf(3),2),num2str(daf(4),2),num2str(daf(5),2));
%gtext(mat2str(dnl));
%set(gca,'xlim',[0 20]);
fin = find(~isnan(af(:,1)));
fin = fin(end);
text('units','data','position',[t(end)-1.5 af(fin,1)-daf(1)+0.3],'string','chord 1')
text('units','data','position',[t(end)-1.5 af(fin,2)-daf(2)],'string','chord 2')
text('units','data','position',[t(end)-1.5 af(fin,3)-daf(3)],'string','chord 3')
text('units','data','position',[t(end)-1.5 af(fin,4)-daf(4)],'string','chord 4')
text('units','data','position',[t(end)-1.5 af(fin,5)-daf(5)],'string','chord 5')

set(gcf,'ResizeFcn','')

h3 =findobj(0,'type','figure','tag','polarplotnc3');
if isempty(h3)
     h3=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
               'times','defaultlinelinewidth',3,'tag','polarplotnc3');
else
     figure(h3);
end
clf

% etalonage de la polarimetrie
afcal = NaN * af;
pp=[];
for k =1:size(af,2)
ind = find(isfinite(afmes(:,k))&isfinite(af(:,k)));
%  pp(k,:)      = polyfit(afmes(ind,k),af(ind,k),1);
%  afmesc(:,k)  = polyval(pp(k,:),afmes(:,k));
%  %legstr{k}    = sprintf(' %4.3g \\alpha ^2  %+4.3g \\alpha %+4.3g',pp(k,:));
   pp(k,2) = mean( af(ind,k) -afmes(ind,k));
   afoff   = af(:,k) - pp(k,2);
   pp(k,1) = mean(afmes(ind,k) ./ afoff(ind));
   afcal(:,k)   = pp(k,1) .* afoff;
  legstr{k}    = sprintf('%5.3g (\\alpha %+5.3g)',pp(k,:));
  ind2 = find(isfinite(nlmes(:,k))&isfinite(nl(:,k)));
  legstr2{k}   = sprintf('%4.2g',sqrt(mean(((nl(ind2,k)-nlmes(ind2,k))./nlmes(ind2,k)).^2)));
end

subplot(2,1,1)
cc=get(gca,'colororder');
plot(t,afcal);
hold on
set(gca,'colororder',cc);
plot(t,afmes,'o');
st =sprintf('choc %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d (calibration)',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);
grid
xlabel('temps (s)');
ylabel('angle Faraday (degres)');

taille = get(gca,'fontsize');
set(gca,'fontsize',12)
legend(gca,legstr,-1);
set(gca,'fontsize',taille)
pos1 =get(gca,'position');

subplot(2,1,2)
set(gca,'colororder',cc);
hold on
plot(t,nl);
hold on
set(gca,'colororder',cc);
plot(t,nlmes,'o');
grid
xlabel('temps (s)');
ylabel('Nl (m^-^3)');
title('Ecart type des Nl')
taille = get(gca,'fontsize');
set(gca,'fontsize',12)
legend(gca,legstr2,-1);
set(gca,'fontsize',taille);
pos2 =get(gca,'position');
pos2(3:4) = pos1(3:4);
set(gca,'position',pos2)
box on
set(gcf,'ResizeFcn','')
