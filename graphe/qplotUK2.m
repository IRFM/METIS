% script de plot du profil d q
chargedata

h = findobj(0,'type','figure','tag','qplot1');
if isempty(h)
       h=figure('tag','qplot1');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',2,'color',[1 1 1])
cc=get(gca,'colororder');
set(gca,'colororder',cc);
subplot(4,1,1)
indx = 1;
plot(t,data.equi.q(:,indx),'-');
hold on
axis([param.gene.tdeb param.gene.t 0.5 3])
set(gca,'colororder',cc);
if ~isempty(t1)
      plot(t1,jeux1.data.equi.q(:,indx),'r--');
end

hold off
grid
ylabel('q')
text('units','normalized','position',[1.01 0.5],'string','x = 0') 
%
%
%
subplot(4,1,2)
indx  = iround(x,0.1);
bile.x=linspace(0,1,21);
indx2 = iround(bile.x,0.1);
plot(t,data.equi.q(:,indx),'-');
hold on
axis([param.gene.tdeb param.gene.t 0.5 2])
set(gca,'colororder',cc);
if ~isempty(t1)
      plot(t1,jeux1.data.equi.q(:,indx),'r--');
end

hold off
grid
ylabel('q')
text('units','normalized','position',[1.01 0.5],'string','x = 0.1') 
%
%
%
subplot(4,1,3)
indx  = iround(x,0.3);
bile.x=linspace(0,1,21);
indx2 = iround(bile.x,0.3);
plot(t,data.equi.q(:,indx),'-');
hold on
axis([param.gene.tdeb param.gene.t 1 2])
set(gca,'colororder',cc);
if ~isempty(t1)
      plot(t1,jeux1.data.equi.q(:,indx),'r--');
end

hold off
grid
ylabel('q')
text('units','normalized','position',[1.01 0.5],'string','x = 0.3') 

%
%
%
subplot(4,1,4)
indx  = iround(x,0.5);
bile.x=linspace(0,1,21);
indx2 = iround(bile.x,0.5);
plot(t,data.equi.q(:,indx),'-');
hold on
axis([param.gene.tdeb param.gene.t 2 3])
set(gca,'colororder',cc);
if ~isempty(t1)
      plot(t1,jeux1.data.equi.q(:,indx),'r--');
end
hold off
grid
ylabel('q')
xlabel('times')
text('units','normalized','position',[1.01 0.5],'string','x = 0.5') 

