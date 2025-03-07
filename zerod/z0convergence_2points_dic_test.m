% test de convergence du modele a 2 points : z0convergence_2points_dic appele avec les 6 methodes
% facteur
zs = post.zerod;
zs.nebord =  15 .* zs.nebord;
zs.pin    =  3 .* zs.pin;
% appel
[tebord0,nelim0,telim0] = z0convergence_2points_dic(post.z0dinput.option,post.z0dinput.cons,post.z0dinput.geo,zs,post.profil0d,0);
[tebord1,nelim1,telim1] = z0convergence_2points_dic(post.z0dinput.option,post.z0dinput.cons,post.z0dinput.geo,zs,post.profil0d,1);
[tebord2,nelim2,telim2] = z0convergence_2points_dic(post.z0dinput.option,post.z0dinput.cons,post.z0dinput.geo,zs,post.profil0d,2);
[tebord3,nelim3,telim3] = z0convergence_2points_dic(post.z0dinput.option,post.z0dinput.cons,post.z0dinput.geo,zs,post.profil0d,3);
[tebord4,nelim4,telim4] = z0convergence_2points_dic(post.z0dinput.option,post.z0dinput.cons,post.z0dinput.geo,zs,post.profil0d,4);
[tebord5,nelim5,telim5] = z0convergence_2points_dic(post.z0dinput.option,post.z0dinput.cons,post.z0dinput.geo,zs,post.profil0d,5);

t = post.zerod.temps;

h = findobj(0,'type','figure','tag','z0plot_test_conv_2pts');
if isempty(h)
       h=figure('tag','z0plot_test_conv_2pts');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(3,1,1)
plot(t,tebord0,'.',t,tebord1,'o',t,tebord2,t,tebord3,'s',t,tebord4,t,tebord5);
legend('dicotomie std.','max','loop','min','random','simulated annealing');
ylabel('Te_a (eV)');
z0loglin(gca);
subplot(3,1,2)
plot(t,telim0,'.',t,telim1,'o',t,telim2,t,telim3,'s',t,telim4,t,telim5);
legend('dicotomie std.','max','loop','min','random','simulated annealing');
ylabel('Te_plate (eV)');
z0loglin(gca);
subplot(3,1,3)
plot(t,nelim0./1e19,'.',t,nelim1./1e19,'o',t,nelim2./1e19,t,nelim3./1e19,'s',t,nelim4./1e19,t,nelim5./1e19);
legend('dicotomie std.','max','loop','min','random','simulated annealing');
ylabel('Ne_plate (1e19 m^-3)');
xlabel('time (s)');
z0loglin(gca);