% verification de l'equation 14 (diffusion du courant)
% FI, 24/04/2002

temps=input('temps ? ');

ii=iround(data.gene.temps,temps)

terme1 = data.equi.grho2r2(ii,:)./param.phys.mu0.*data.neo.eta(ii,:)./data.equi.rhomax(ii).^2./ ...
         data.equi.r2i(ii,:).*data.prof.psid2(ii,:);
x   = param.gene.x;
%a   = (data.equi.vpr(ii,3)-data.equi.vpr(ii,2))/(x(3)-x(2));
%b   = data.equi.vpr(ii,3)-a*x(3);
%data.equi.vpr(ii,1)=a*x(1)+b;
lnc = log(data.equi.vpr(ii,:).*data.equi.grho2r2(ii,:)./data.equi.F(ii,:));
dlnc = gradient(lnc,x);

%terme2 = data.equi.grho2r2(ii,:)./param.phys.mu0.*data.neo.eta(ii,:)./data.equi.rhomax(ii).^2./ ...
%         data.equi.r2i(ii,:).*(dlnc+param.gene.x./data.equi.rhomax(ii).* ...
%			data.equi.drhomaxdt(ii)).*data.prof.psid1(ii,:);
terme2 = (data.equi.grho2r2(ii,:)./param.phys.mu0.*data.neo.eta(ii,:)./data.equi.rhomax(ii).^2./ ...
         data.equi.r2i(ii,:).*dlnc + param.gene.x./data.equi.rhomax(ii).* ...
			data.equi.drhomaxdt(ii)).*data.prof.psid1(ii,:);
terme3 = data.geo.b0(ii).*data.neo.eta(ii,:)./data.equi.F(ii,:)./data.equi.r2i(ii,:).* ...
         data.source.totale.j(ii-1,:);

h = findobj(0,'type','figure','tag','eq_14');
if isempty(h)
       h=figure('tag','eq_14');
else
       figure(h);
end   

subplot(3,1,1)
%plot(param.gene.x,terme1,param.gene.x,terme2,param.gene.x,terme3)
%subplot(2,1,2)
plot(param.gene.x,data.prof.dpsidt(ii,:),param.gene.x,terme1+terme2+terme3,'ro',param.gene.x,data.prof.dpsidt(ii,:),'m-.')
legend('terme gauche eq 14','terme droit eq 14')
subplot(3,1,2)
% ecart homogene a une densite de courant
jequierr = (terme1+terme2+terme3-data.prof.dpsidt(ii,:)) ./(data.geo.b0(ii).*data.neo.eta(ii,:)./data.equi.F(ii,:)./data.equi.r2i(ii,:));
plot(param.gene.x,data.prof.jeff(ii,:),param.gene.x,data.source.totale.j(ii-1,:)+data.prof.epar(ii,:)./data.neo.eta(ii,:),param.gene.x,jequierr)
legend('jmoy','jni+E/eta','j erreur eq 14')
subplot(3,1,3)
jerr    = abs(data.prof.jeff(ii,:) - data.source.totale.j(ii-1,:)-data.prof.epar(ii,:)./data.neo.eta(ii,:));
semilogy(param.gene.x,abs(data.prof.dpsidt(ii,:)-terme1-terme2-terme3)./ ...
    abs(diff(data.prof.psi(ii,[1,end]))),'ro',x,abs(jerr./data.prof.jmoy(ii,:)),'b')
