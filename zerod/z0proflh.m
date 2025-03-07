% visulaisation 3D des profils 0D
zs   = post.zerod;
op0d = post.z0dinput.option;
cons = post.z0dinput.cons;
geo = post.z0dinput.geo;
t    = zs.temps;
% profil de courant et de puissance
swlh = abs(zs.plh ./ max(1,zs.pohm)) > 0.2;
pfweh = max(0,min(1,op0d.fwcd ~=0)) .* zs.picrh_th;
picrh = zs.picrh_th - pfweh;
profli = post.profil0d;
xli    = profli.xli;
t      = profli.temps;
jlh    = profli.jlh ./ max(1,max(profli.jlh,[],2) * ones(size(xli)));
plh    = profli.plh ./ max(1,max(profli.plh,[],2) * ones(size(xli)));
flh    = profli.fxdurplh ./ max(eps,max(profli.fxdurplh,[],2) * ones(size(xli)));

h = findobj(0,'type','figure','tag','z0proflh');
if isempty(h)
       h=figure('tag','z0proflh');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
colormap('hot')


%flag = 0;
%k = 3;
%subplot(k,k-1,1)
zplotprof(gca,t,xli,flh,'color','r','marker','o');
zplotprof(gca,t,xli,jlh,'color','b');
zplotprof(gca,t,xli,plh,'color','g');
if exist('data','var')
	jcr = data.source.hyb.j ./ max(1,max(data.source.hyb.j,[],2) * ones(size(param.gene.x)));
	zplotprof(gca,data.gene.temps,param.gene.x,jcr,'color','m');
end
title('Xdur (r), Jlh (b) & Plh (g) (maximum normalized at 1)')
xlabel('x (normalized radius))');





