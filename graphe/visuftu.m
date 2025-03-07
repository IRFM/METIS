h = findobj(0,'type','figure','tag','flux neutrons');
if isempty(h)
       h=figure('tag','flux neutrons');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
if strcmp(param.from.machine,'TS')

  [fnexp,tfnexp]=tsbase(fix(param.from.shot.num),'gfluntn');

else

  fnexp = [];
  tfnexp = [];

end
 
% coordonnees cronos r/a 
rsa              =  data.equi.a ./ (data.equi.a(:,end) * ones(1,size(data.equi.a,2)));
indd             = min(find((param.compo.z == 1) & (param.compo.a == 2)));

xfn              = param.gene.x;
xxfn             = ones(length(data.gene.temps),1) * xfn;
nd               = squeeze(data.impur.impur(:,:,indd));
nd               = tsplinet(rsa,nd,xxfn);
rm               = tsplinet(rsa,data.equi.raxe,xxfn);
ti               = tsplinet(rsa,data.prof.ti,xxfn);
flux.temps       = data.gene.temps;
sv               = sigvddn(data.prof.ti/1e3); 
flux.neutron     = 1/2 * data.equi.rhomax .* trapz(xfn,data.equi.vpr .* sv .* nd .* nd,2);


flux.neutronthe  = data.equi.rhomax .* trapz(xfn,data.equi.vpr .* data.source.fus.neutron.dd,2);
flux.neutrontot  = data.equi.rhomax .* trapz(xfn,data.equi.vpr .* data.source.totale.neutron.dd,2);


subplot(2,1,1)
if ~isempty(fnexp)
  plot(flux.temps,flux.neutron,'ro',flux.temps,flux.neutronthe,'b',...
     flux.temps,flux.neutrontot,'m+',tfnexp,10.^fnexp(:,1),'k:')
  legend('formule','thermique','total','exp')
else
  plot(flux.temps,flux.neutron,'ro',flux.temps,flux.neutronthe,'b',...
     flux.temps,flux.neutrontot,'m+')
  legend('formule','thermique','total')

end
xlabel('temps');
ylabel('flux neutron')


subplot(2,1,2)

plotprof('gca',flux.temps,xfn,data.source.fus.neutron.dd,'color','b')
plotprof('gca',flux.temps,xfn,data.source.totale.neutron.dd,'color','g','linestyle','--')

title('profil flux neutrons')
