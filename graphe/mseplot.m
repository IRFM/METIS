% plot de la mse de jet
if isempty(post)==1

   disp('You need to run the post processing tool')
   return

end
chargedata
x=param.gene.x;
t=data.gene.temps;

if length(t)> 100
	pas =fix(length(t)/100);
	if pas ==0
		pas =1;
	end
else
	pas =1;
end

ind=param.gene.kmin:min(param.gene.kmax,length(t));
tt=t(ind);
inde = 1:pas:length(t);
t=t(inde);

if exist('jeux1','var')
  if ~isempty(jeux1)
    ind1=jeux1.param.gene.kmin:min(jeux1.param.gene.kmax,length(t1));
    tt1=t1(ind1);
	 if tt1(1) > 35 & tt(1) < 40 & strcmp(param.from.machine,'JET')
	   tt1 = tt1 - 42;
	 end
  else
    clear jeux1
  end
end
h = findobj(0,'type','figure','tag','mse');
if isempty(h)
       h=figure('tag','mse');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

nn = 1:size(post.mse.angle,2);
strchoc = int2str(param.from.shot.num);
if ~exist('post')
  return
end
if isfield(post.mse,'rmse')
   if strcmp(param.from.machine,'JET')
     %mse=load(strcat('/usr/drfc/cgc/matlab5/tcron/JET/data/',strchoc,'/farmse',strchoc));
     mse = load(strcat(getenv('HOME'),'/zineb/data/JET/',strchoc,'/farmse',strchoc));
     post.mse.rmse=mse.R1(1,:);
   end
   if strcmp(param.from.machine,'DIIID')

     mse = post.mse.mesure;
   end 
end
t = data.gene.temps;
r = post.mse.rmse;
%figure('name','MSE JET','color',[1 1 1],'defaultaxesfontsize',12);
plotprof(gca,t,r,post.mse.angle.*180/pi,'color','r','linestyle','*');
plotprof(gca,t,r,post.mse.mesure.*180/pi,'color','k','linestyle','o');
axis([-inf inf -inf inf])
if exist('jeux1','var')
  plotprof(gca,tt1,jeux1.post.mse.rmse,jeux1.post.mse.angle,'color','b','linestyle','--')
  plotprof(gca,tt1,jeux1.post.mse.rmse,jeux1.post.mse.mesure,'color','b','linestyle','p');
  legend('simulation','exp.','ref. sim','ref. mes')
else
  legend('simulation','exp.')
end
st =sprintf('shot %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d',param.from.machine, ...
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);
grid
xlabel('R (m)')
ylabel('Angle (deg.)')
set(gcf,'ResizeFcn','')

h = findobj(0,'type','figure','tag','msetemp');
if isempty(h)
       h=figure('tag','msetemp');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
[rmse,irmse] = sort(post.mse.rmse);
angle = post.mse.angle(:,irmse).*180/pi;
mesure = post.mse.mesure(:,irmse).*180/pi;
for i=1:9
eval(['subplot(3,3,',int2str(i),')'])

if strcmp(param.from.machine,'JET')
   channeltab = [4 7:2:22];
   j = channeltab(i);    % choose 9 channels among the 25 available extreme channels (edge and core are not trusted by the diagnostician)
else
   j = 5*i-4;
   if j == 26
     j = 27;  
   end
   if j == 41
     j = 39;
   end
end

eval(['plot(t,angle(:,',int2str(j),'),t,mesure(:,',int2str(j),'),''.'')'])
title(['R=',num2str(rmse(j),3)])
end


