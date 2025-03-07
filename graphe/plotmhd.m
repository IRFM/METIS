% script de tracer de la stbilite (mishka)
varliste = plotmhd1t;
% callback = fontion a appeler ma_fonction(..,temps_zslider,...,option_zslider...)
callback = 'plotmhd1t(temps_zslider,option_zslider,post.mhd,param);';

% vecteur temps
temps =[];
if isfield(post.mhd,'n1');
   temps = union(temps,post.mhd.n1.t);
end
if isfield(post.mhd,'n2');
   temps = union(temps,post.mhd.n2.t);
end
if isfield(post.mhd,'n2');
   temps = union(temps,post.mhd.n2.t);
end

h = findobj(0,'type','figure','tag','plotmhd');
if isempty(h)
       h=figure('tag','plotmhd');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])


zslider(h,temps,varliste,2,callback);

h = findobj(0,'type','figure','tag','mercier');
if isempty(h)
       h=figure('tag','mercier');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])


subplot(4,1,1)
t   = data.gene.temps;
x   = param.gene.x;
mhd = data.equi.mhd;
zplotprof(gca,t,x,mhd.ballooning,'color','r');
zplotprof(gca,t,x,mhd.mercier,'color','b');
zplotprof(gca,t,x,mhd.ideal,'color','g');
ylabel('su');
ylabel('indice de stabilite (stable si < 0)');
legend('ballooning','Mercier','Ideal');

subplot(4,1,2)
mhd = post.mhd;
if isfield(mhd,'n1')
   zplotprof(gca,mhd.n1.t,x,mhd.n1.FM,'color','r');
elseif isfield(mhd,'n2')
   zplotprof(gca,mhd.n1.t,x,mhd.n2.FM,'color','b');
elseif isfield(mhd,'n2')
   zplotprof(gca,mhd.n1.t,x,mhd.n3.FM,'color','g');
end
ylabel('ballooning');

subplot(4,1,3)
mhd = post.mhd;
if isfield(mhd,'n1')
   zplotprof(gca,mhd.n1.t,x,mhd.n1.DR,'color','r');
elseif isfield(mhd,'n2')
   zplotprof(gca,mhd.n1.t,x,mhd.n2.DR,'color','b');
elseif isfield(mhd,'n2')
   zplotprof(gca,mhd.n1.t,x,mhd.n3.DR,'color','g');
end
ylabel('Mercier');

subplot(4,1,4)
mhd = post.mhd;
if isfield(mhd,'n1')
   zplotprof(gca,mhd.n1.t,x,mhd.n1.DI,'color','r');
elseif isfield(mhd,'n2')
   zplotprof(gca,mhd.n1.t,x,mhd.n2.DI,'color','b');
elseif isfield(mhd,'n2')
   zplotprof(gca,mhd.n1.t,x,mhd.n3.DI,'color','g');
end
ylabel('Ideal');
xlabel('x (su)');
