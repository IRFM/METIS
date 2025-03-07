h = findobj(0,'type','figure','tag','neutron flux');
if isempty(h)
       h=figure('tag','neutron flux');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
if strcmp(param.from.machine,'TS')

  [fnexp,tfnexp]=tsbase(fix(param.from.shot.num),'gfluntn');
  fnexp         = 10.^smooth(fnexp,10);

elseif  strcmp(param.from.machine,'JET')
  try
    racine = evalin('base','prepare.racine');
  catch
	 if strcmp(getenv('USER'),'cgc')
		racine ='/usr/drfc/cgc/cgc_data/jet/data/';
	 else
		racine = strcat(getenv('HOME'),'/zineb/data/JET/');
	end
%	racine ='/usr/drfc/cgc/matlab5/tcron/JET/data/';
  end

  chemin = [racine,int2str(param.from.shot.num),'/'];
  temp = load([chemin,'temp',int2str(param.from.shot.num)]);
  if ~isempty(temp.rnt)
    fnexp = temp.rnt;
    tfnexp = temp.trnt;
  end

elseif  strcmp(param.from.machine,'DIIID')

  racine = strcat(getenv('HOME'),'/zineb/data/diiid/');
  chemin = [racine,int2str(param.from.shot.num),'/'];
  load([chemin,'diiiddiag']);
  if ~isempty(diiiddiag.neut) 

    fnexp = diiiddiag.neut;
    tfnexp = diiiddiag.tneut;

  end
end
 
% coordonnees cronos r/a 
rsa              =  data.equi.a ./ (data.equi.a(:,end) * ones(1,size(data.equi.a,2)));
indd             =  min(find((param.compo.z == 1) & (param.compo.a == 2)));

xfn              = param.gene.x;
xxfn             = ones(length(data.gene.temps),1) * xfn;
nd               = squeeze(data.impur.impur(:,:,indd));
nd               = tsplinet(rsa,nd,xxfn);
rm               = tsplinet(rsa,data.equi.raxe,xxfn);
ti               = tsplinet(rsa,data.prof.ti,xxfn);
flux.temps       = data.gene.temps;
sv               = sigvddn(data.prof.ti/1e3); 
flux.neutron     = 1/2 * data.equi.rhomax .* trapz(xfn,data.equi.vpr .* sv .* nd .* nd,2);
flux.nd          = nd;
flux.ti          = ti;
data.source.fus.neutron.dd(data.source.fus.neutron.dd>1e20) = 0;
data.source.totale.neutron.dd(data.source.totale.neutron.dd>1e20) = 0;

flux.neutronthe  = data.equi.rhomax .* trapz(xfn,data.equi.vpr .* data.source.fus.neutron.dd,2);
flux.neutrontot  = data.equi.rhomax .* trapz(xfn,data.equi.vpr .* data.source.totale.neutron.dd,2);


subplot(2,1,1)
if ~isempty(fnexp)
  plot(flux.temps,flux.neutron,'ro',flux.temps,smooth(flux.neutronthe,3),'b',...
     flux.temps,smooth(flux.neutrontot,3),'m+',tfnexp,fnexp(:,1),'k:')
  legend('from Ti','thermal','total','exp')
  axis([data.gene.temps(1) data.gene.temps(param.gene.k) 0 max(max(fnexp(:,1)),max(flux.neutrontot))])
else
  plot(flux.temps,flux.neutron,'ro',flux.temps,flux.neutronthe,'b',...
     flux.temps,flux.neutrontot,'m+')
  legend('from Ti','thermal','total')
  axis([data.gene.temps(1) data.gene.temps(param.gene.k) 0 max(fnexp(:,1))])

end
xlabel('time (s)');
ylabel('neutron flux')


subplot(2,1,2)

plotprof('gca',flux.temps,xfn,data.source.fus.neutron.dd,'color','b')
plotprof('gca',flux.temps,xfn,data.source.totale.neutron.dd,'color','g','linestyle','--')

title('neutron flux profile')
