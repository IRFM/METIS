% script de verification de la coherence entre les profils de puissances et les donnees 0d
pradc  = trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.prad,2);
pbremc = trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.pbrem,2);
pcycloc = trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.pcyclo,2);
pfusc  = trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.pfus,2);
picrhc = trapz(post.profil0d.xli,post.profil0d.vpr .* (post.profil0d.picrh +post.profil0d.pfweh) ,2);
plhc = trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.plh,2);
pnbic = trapz(post.profil0d.xli,post.profil0d.vpr .* (real(post.profil0d.pnbi) + imag(post.profil0d.pnbi)),2);
pecrhc = trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.pecrh,2);
pohmc = trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.pohm,2);

% current sources
inbicdc = trapz(post.profil0d.xli,post.profil0d.spr .* (real(post.profil0d.jnbicd) + imag(post.profil0d.jnbicd)),2);
ilhc    = trapz(post.profil0d.xli,post.profil0d.spr .* post.profil0d.jlh,2);
ieccdc  = trapz(post.profil0d.xli,post.profil0d.spr .* post.profil0d.jeccd,2);
ifwcdc  = trapz(post.profil0d.xli,post.profil0d.spr .* post.profil0d.jfwcd,2);
irunc   = trapz(post.profil0d.xli,post.profil0d.spr .* post.profil0d.jrun,2);
ibootc  = trapz(post.profil0d.xli,post.profil0d.spr .* post.profil0d.jboot,2);


h = findobj(0,'type','figure','tag','z0plotcohe');
if isempty(h)
       h=figure('tag','z0plotcohe');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(4,1,1);
plot(post.profil0d.temps,pradc./1e6,'or',post.zerod.temps,post.zerod.prad./1e6,'r', ...
     post.profil0d.temps,pbremc./1e6,'ob',post.zerod.temps,post.zerod.pbrem./1e6,'b', ...
     post.profil0d.temps,pcycloc./1e6,'og',post.zerod.temps,post.zerod.pcyclo./1e6,'g');
z0loglin(gca);
title(sprintf('METIS : %s@%d/coherence puissance 0D/1D ', post.z0dinput.machine,post.z0dinput.shot));
legend('Line 1D','Line 0D','Brem 1D','Brem 0D','Cyclo 1D','Cyclo 0D')
ylabel('MW')
subplot(4,1,2);
plot(post.profil0d.temps,pfusc./1e6,'or',post.zerod.temps,post.zerod.pfus_th./1e6,'r', ...
     post.profil0d.temps,pnbic./1e6,'ob',post.zerod.temps,real(post.zerod.pnbi_th) ./1e6 + imag(post.zerod.pnbi_th)./1e6,'b', ...
     post.profil0d.temps,picrhc./1e6,'og',post.zerod.temps,post.zerod.picrh_th./1e6,'g');
legend('Alpha 1D','Alpha 0D','NBI 1D','NBI 0D','ICRH 1D','ICRH 0D');
z0loglin(gca);
ylabel('MW')
subplot(4,1,3);
plot(post.profil0d.temps,pecrhc./1e6,'or',post.zerod.temps,post.zerod.pecrh./1e6,'r', ...
     post.profil0d.temps,plhc./1e6,'ob',post.zerod.temps,post.zerod.plh./1e6,'b', ...
     post.profil0d.temps,pohmc./1e6,'og',post.zerod.temps,post.zerod.pohm./1e6,'g');
legend('ECRH 1D','ECRH 0D','LH 1D','LH 0D','Ohm 1D','Ohm 0D');
z0loglin(gca);
ylabel('MW')
subplot(4,1,4);
plot(post.profil0d.temps,ieccdc./1e6,'or',post.zerod.temps,post.zerod.ieccd./1e6,'r', ...
     post.profil0d.temps,ilhc./1e6,'ob',post.zerod.temps,post.zerod.ilh./1e6,'b', ...
     post.profil0d.temps,inbicdc./1e6,'ok',post.zerod.temps,(real(post.zerod.inbicd) + imag(post.zerod.inbicd)) ./1e6,'k', ...
     post.profil0d.temps,ibootc./1e6,'oc',post.zerod.temps,post.zerod.iboot./1e6,'c', ...
     post.profil0d.temps,irunc./1e6,'om',post.zerod.temps,post.zerod.irun./1e6,'m');
legend('ECRH 1D','ECRH 0D','LH 1D','LH 0D','NBI 1D','NBI 0D','Boot 1D','Boot 0D','Runaway 1D','Runaway 0D');
z0loglin(gca);
ylabel('MA')
xlabel('time (s)');
joint_axes(h,4);
edition2
