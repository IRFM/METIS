% comparaison des donnees de la chaleur en fonction du temps
chargedata

switch param.from.machine

case 'TS'
     
     ndd_tot = data.equi.rhomax .* trapz(x,data.equi.vpr .* data.source.totale.neutron.dd,2);
     %ndd_tot = medfilt1(ndd_tot,5); disp('attention filtrage provisoire ')
     if ~isempty(t1)
        jeux1.ndd_tot = jeux1.data.equi.rhomax .* trapz(x,jeux1.data.equi.vpr .* jeux1.data.source.totale.neutron.dd,2);
        %jeux1.ndd_tot = medfilt1(jeux1.ndd_tot,5); disp('attention filtrage provisoire ')
     end
     h       = findobj(0,'type','figure','tag','neutron_dd_comparaison');
     if isempty(h)
         h = figure('tag','neutron_dd_comparaison');
     else
         figure(h);
     end
     clf
     set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	   'defaultlinelinewidth',3,'color',[1 1 1])
     
     subplot(3,1,1)
     if ~isempty(t1)
        plot(t,ndd_tot,'k',t1,jeux1.ndd_tot,'bo',tndd_exp,ndd_exp,'r')
        legend('Cronos','Reference','Experience')
     else
        plot(t,ndd_tot,'k',tndd_exp,ndd_exp,'r')
        legend('Cronos','Experience')
     end
     set(gca,'xlim',[min(t),max(t)])
     xlabel('temps (s)')
     ylabel('Neutron (Hz)')
     title('Flux total de neutron dd')

     subplot(3,1,2)
     if ~isempty(t1) & ~isempty(tiprof)
        plot(t,data.prof.ti(:,1),'k',t1,jeux1.data.prof.ti(:,1),'bo',tiprof.t,tiprof.Ti(:,1)*1e3,'r')
        legend('Cronos','Reference','Experience')
     elseif ~isempty(t1) 
        plot(t,data.prof.ti(:,1),'k',t1,jeux1.data.prof.ti(:,1),'bo')
        legend('Cronos','Reference')
     elseif ~isempty(tiprof)
        plot(t,data.prof.ti(:,1),'k',tiprof.t,tiprof.Ti(:,1)*1e3,'r')
        legend('Cronos','Experience')
     else
        plot(t,data.prof.ti(:,1),'k')
        legend('Cronos')
     end
     xlabel('temps (s)')
     ylabel('Ti (eV)')
     set(gca,'xlim',[min(t),max(t)])
     title('Ti centrale')

     subplot(3,1,3)
     indd = find((param.compo.z == 1)&(param.compo.a == 2));
     if ~isempty(indd)
	if ~isempty(t1) & ~isempty(tiprof)
           plot(t,data.impur.impur(:,1,indd),'k',t1,jeux1.data.impur.impur(:,1,indd),'bo',tiprof.t,tiprof.nD,'r', ...
                t,data.prof.ne(:,1),'m:');
           legend('Cronos','Reference','Experience','Ne0')
	elseif ~isempty(t1) 
           plot(t,data.impur.impur(:,1,indd),'k',t1,jeux1.data.impur.impur(:,1,indd),'bo', ...
                t,data.prof.ne(:,1),'m:');
           legend('Cronos','Reference','Ne0')
	elseif ~isempty(tiprof)
           plot(t,data.impur.impur(:,1,indd),'k',tiprof.t,tiprof.nD,'r', ...
                t,data.prof.ne(:,1),'m:');
           legend('Cronos','Experience','Ne0')
	else
           plot(t,data.impur.impur(:,1,indd),'k', ...
                t,data.prof.ne(:,1),'m:');
           legend('Cronos','Ne0')
	end
     end
     xlabel('temps (s)')
     ylabel('nD (m^-3)')
     title('Densite de deuterium')
     set(gca,'xlim',[min(t),max(t)])
     
otherwise

    disp('Machine inconnue')
end
