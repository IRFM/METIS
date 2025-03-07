% scrit pour comparer le Wi_eff cronos et metis
h = findobj(0,'type','figure','tag','z0chieff');
if isempty(h)
       h=figure('tag','z0chieff');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
colormap('hot')


chieff_metis = (post.profil0d.xie .* post.profil0d.nep + post.profil0d.xii .* post.profil0d.nip) ./ ...
              (post.profil0d.nep + post.profil0d.nip);
Qeff_metis       = cumtrapz(post.profil0d.xli,(post.profil0d.source_el + post.profil0d.source_ion) .* post.profil0d.vpr,2) ./ post.profil0d.grho2 ./ max(eps,post.profil0d.vpr_tor);

chieff_cronos = (data.prof.flux.keeff + data.prof.flux.kieff) ./ ...
               data.prof.ne ./ (1 + data.prof.ae);
Qeff_cronos  = data.prof.flux.qe +  data.prof.flux.qi;

subplot(2,2,1)
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,chieff_metis,'color','b');
zplotprof(gca,data.gene.temps,param.gene.x,chieff_cronos,'color','r');
ylabel('Chi_{eff}')
set(gca,'yscale','log')
axis([0,1,0,100]);
title(sprintf('Transport ananlysis : %s @ %d', ...
             post.z0dinput.machine,post.z0dinput.shot));
subplot(2,2,2)
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,Qeff_metis,'color','b');
zplotprof(gca,data.gene.temps,param.gene.x,Qeff_cronos,'color','r');
ylabel('Q_{eff}')
legend('METIS','CRONOS');

%
chieff_metis_f= sgolayfilt(chieff_metis,1,11);
chieff_cronos_f= sgolayfilt(chieff_cronos,1,11);
Qeff_metis_f= sgolayfilt(Qeff_metis,1,11);
Qeff_cronos_f= sgolayfilt(Qeff_cronos,1,11);

subplot(2,2,3)
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,chieff_metis_f,'color','b');
zplotprof(gca,data.gene.temps,param.gene.x,chieff_cronos_f,'color','r');
set(gca,'yscale','log')
axis([0,1,0,100]);
ylabel('Chi_{eff} filtered')
xlabel('\rho')
subplot(2,2,4)
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,Qeff_metis_f,'color','b');
zplotprof(gca,data.gene.temps,param.gene.x,Qeff_cronos_f,'color','r');
ylabel('Q_{eff} filtered')
xlabel('\rho')




