h = findobj(0,'type','figure','tag','compare');
if isempty(h)
       h=figure('tag','compare');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
paroi = param.from.paroi;


Rext=data.geo.R(200,:);
Zext=data.geo.Z(200,:);
RR=squeeze(double(data.equi.R(200,1:10:101,:)));
ZZ=squeeze(double(data.equi.Z(200,1:10:101,:)));

plot(paroi.R,paroi.Z,Rext,Zext,RR',ZZ')
axis('square')
xlabel('R(m)')
ylabel('Z(m)')



temps = data.gene.temps;
x = param.gene.x;
temps2 = jeux2.data.gene.temps;


%
% predictif seul
%
h1=subplot(4,1,1);
plot(temps2,jeux2.data.gene.nbar,'b')
axis([0 100 1e19 2.5e19])
grid
set(h1,'xticklabel',[]) 
h2=subplot(4,1,2);

plot(temps2,jeux2.data.gene.ihyb,'b')
axis([0 100 2e5 7e5])
grid
set(h2,'xticklabel',[]) 

h3=subplot(4,1,3);

plot(temps2,jeux2.data.gene.li,'b')
axis([0 100 1.2 2])
grid
set(h3,'xticklabel',[]) 

h4=subplot(4,1,4);

plot(temps2,jeux2.data.gene.tau,'b')
axis([0 100 0.04 0.1])

grid
%
% predicitf + experience
%
h1=subplot(4,1,1);
plot(temps,data.gene.nbar,'r',temps2,jeux2.data.gene.nbar,'b')
axis([0 100 1e19 2.5e19])
grid
set(h1,'xticklabel',[]) 
h2=subplot(4,1,2);

plot(temps,data.gene.ihyb,'r',temps2,jeux2.data.gene.ihyb,'b',temps,jeux4.data.gene.ihyb,'m')
axis([0 100 2e5 7e5])
grid
set(h2,'xticklabel',[]) 

h3=subplot(4,1,3);

plot(temps,data.gene.li,'r',temps2,jeux2.data.gene.li,'b',temps,jeux4.data.gene.li,'m')
axis([0 100 1.2 2])
grid
set(h3,'xticklabel',[]) 

h4=subplot(4,1,4);

plot(temps,data.gene.tau,'r',temps2,jeux2.data.gene.tau,'b',temps,jeux4.data.gene.tau,'m')
axis([0 100 0.04 0.1])

grid
