% useful tools for Frederic (April 2006)
 
choc = 123042;
[shoto,status]=mdsopen(['atlas.gat.com::EFIT03'],choc);
efit.q=mdsvalue('\EFIT03::qpsi');
efit.q0=mdsvalue('\EFIT03::q0');
efit.qmin=mdsvalue('\EFIT03::qmin');
efit.rho=mdsvalue('\EFIT03::rhovn');
efit.t = mdsvalue('dim_of(\EFIT03::q0)');

% To set the reference value of qmin to the EFIT one
data.cons.asser.qmin=interp1(efit.t/1000,efit.qmin,data.gene.temps);
data.cons.asser.qmin(end)=2.4052;   % small bug for 123042, fix this later

temps=2.8; %(s)
figure
itc=iround(data.gene.temps,temps);
ite=iround(efit.t,temps*1000);
plot(param.gene.x,data.prof.q(itc,:),'b-.',param.gene.x,jeux1.data.prof.q(itc,:),'g-',efit.rho(:,ite),efit.q(:,ite),'r--')
xlabel('\rho')
ylabel('q')
title(['Shot ',int2str(choc),', t = ',num2str(temps),' s'])
legend('CRONOS, modified NBI','CRONOS, original NBI','EFIT')

qmin=data.gene.temps.*0;
%jeux1qmin = qmin;
for i=1:length(data.gene.temps)
   qmin(i) = min(data.prof.q(i,:));
%   jeux1qmin(i) = min(jeux1.data.prof.q(i,:));
end

figure
plot(data.gene.temps,qmin,'b-.',data.gene.temps,data.cons.asser.qmin,'g-')
%plot(data.gene.temps,qmin,'b-.',data.gene.temps,jeux1qmin,'g-',efit.t/1000,efit.qmin,'r--')
%title(['Shot ',int2str(choc)])
xlabel('Time (s)')
ylabel('qmin')
legend('result','target')
%legend('CRONOS, modified NBI','CRONOS, original NBI','EFIT')

figure
plot(data.gene.temps,data.gene.paddidn,'b-.',data.gene.temps,sum(jeux1.data.cons.idn')','g-')
xlabel('Time (s)')
ylabel('NBI power (W)')
legend('result','feed forward')
