% script de test de LUKE dans metis	
time = post.profil0d.temps(post.lukeinmetis.computation_time_slice);
ind0d = find(post.zerod.temps >= time,1);		     
% le plot
hz =findobj(0,'type','figure','tag','test_lukeinmetis');
if isempty(hz)
  	  hz=figure('tag','test_lukeinmetis','name','LUKE in METIS');
else
  	  figure(hz);
end
clf
set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(2,2,1);
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.jlh + post.profil0d.jeccd,'color','b');
hold on
plot(post.profil0d.xli, post.lukeinmetis.j_TOT - post.lukeinmetis.j_OHM,'r');
%plot(post.profil0d.xli, -(post.lukeinmetis.j_TOT - post.lukeinmetis.j_OHM),'m');
ylabel('J_{LH} + J_{ECCD} + synergy');
xlabel('x (normalized radius)')
title('Current source');
legend(sprintf('METIS (%g A)',post.zerod.ilh(ind0d) + post.zerod.ieccd(ind0d)), ...
sprintf('LUKE (%g A)',post.lukeinmetis.I_TOT - post.lukeinmetis.I_OHM));

subplot(2,2,2);
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.jeff - post.profil0d.jni,'color','b');
hold on
plot(post.profil0d.xli,post.lukeinmetis.j_OHM,'r');
%plot(post.profil0d.xli,-post.lukeinmetis.j_OHM,'m');
ylabel('J_{Ohm}');
xlabel('x (normalized radius)')
title('Ohmic current');
legend(sprintf('METIS (%g A)',post.zerod.iohm(ind0d)), ...
sprintf('LUKE (%g A)',post.lukeinmetis.I_OHM));

subplot(2,2,3);
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.plh + post.profil0d.pecrh,'color','b');
hold on
plot(post.profil0d.xli,post.lukeinmetis.p_TOT - post.lukeinmetis.p_OHM,'r');
ylabel('P_{LH} + P_{ECCD} ');
xlabel('x (normalized radius)')
title('Heat source');
legend(sprintf('METIS (%g W)',post.zerod.plh(ind0d) + post.zerod.pecrh(ind0d)), ...
sprintf('LUKE (%g W)',post.lukeinmetis.P_TOT - post.lukeinmetis.P_OHM));

subplot(2,2,4);
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.pohm,'color','b');
hold on
plot(post.profil0d.xli,post.lukeinmetis.p_OHM,'r');
ylabel('P_{Ohm}');
xlabel('x (normalized radius)')
title('Ohmic power');
legend(sprintf('METIS (%g W)',post.zerod.pohm(ind0d)), ...
sprintf('LUKE (%g W)',post.lukeinmetis.P_OHM));

htemps=findobj(hz,'tag','temps');
set(htemps,'string',sprintf('%g',post.profil0d.temps(post.lukeinmetis.computation_time_slice)));
zplotprof('temps');
