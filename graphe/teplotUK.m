% script de plot du profil d q
chargedata

h = findobj(0,'type','figure','tag','Te');
if isempty(h)
       h=figure('tag','Te');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',2,'color',[1 1 1])
cc=get(gca,'colororder');
set(gca,'colororder',cc);
subplot(4,1,1)
indx  = iround(x,0.1);

plot(t,data.prof.te(:,indx)/1e3,'-');
hold on
axis([param.gene.tdeb param.gene.t 4 8])
set(gca,'colororder',cc);
if ~isempty(t1)
      plot(t1,jeux1.data.prof.te(:,indx)/1e3,'r--');
end
plot([0 1],[1 1],'g:');
plot([0 1],[1.5 1.5],'c:');
plot([0 1],[2 2],'b:');
hold off
grid
ylabel('keV')
text('units','normalized','position',[1.01 0.5],'string','x = 0.1') 
%
%
%
subplot(4,1,2)
indx  = iround(x,0.3);
bile.x=linspace(0,1,21);
indx2 = iround(bile.x,0.1);
plot(t,data.prof.te(:,indx)/1e3,'-');
hold on
axis([param.gene.tdeb param.gene.t 2 4])
set(gca,'colororder',cc);
if ~isempty(t1)
      plot(t1,jeux1.data.prof.te(:,indx)/1e3,'r--');
end
plot([0 1],[1 1],'g:');
plot([0 1],[1.5 1.5],'c:');
plot([0 1],[2 2],'b:');
hold off
grid
ylabel('keV')
text('units','normalized','position',[1.01 0.5],'string','x = 0.3') 
%
%
%
subplot(4,1,3)
indx  = iround(x,0.5);
bile.x=linspace(0,1,21);
indx2 = iround(bile.x,0.3);
plot(t,data.prof.te(:,indx)/1e3,'-');
hold on
axis([param.gene.tdeb param.gene.t 0.2 2])
set(gca,'colororder',cc);
if ~isempty(t1)
      plot(t1,jeux1.data.prof.te(:,indx)/1e3,'r--');
end
plot([0 1],[1 1],'g:');
plot([0 1],[1.5 1.5],'c:');
plot([0 1],[2 2],'b:');
hold off
grid
ylabel('keV')
text('units','normalized','position',[1.01 0.5],'string','x = 0.5') 

%
%
%
subplot(4,1,4)
indx  = iround(x,0.7);
bile.x=linspace(0,1,21);
indx2 = iround(bile.x,0.5);
plot(t,data.prof.te(:,indx)/1e3,'-');
hold on
axis([param.gene.tdeb param.gene.t 0 1])
set(gca,'colororder',cc);
if ~isempty(t1)
      plot(t1,jeux1.data.prof.te(:,indx)/1e3,'r--');
end
plot([0 1],[1 1],'g:');
plot([0 1],[1.5 1.5],'c:');
plot([0 1],[2 2],'b:');
hold off
grid
ylabel('keV')
xlabel('times')
text('units','normalized','position',[1.01 0.5],'string','x = 0.7') 

