d   =abs(data.gene.temps-14);
ind = min(find(min(d)==d));
t =data.gene.temps(1:ind);
ip =data.cons.ip(1:ind)./1e6;
nbar = data.gene.nbar(1:ind)./1e19;
pfci = data.gene.paddfci(1:ind)./1e6;
iboot = data.gene.iboot(1:ind)./1e6;
ilh = data.gene.ihyb(1:ind)./1e6;
te    = data.prof.te(1:ind,:)./1e3;
ti    = data.prof.ti(1:ind,:)./1e3;

vcal  = data.gene.vres(1:ind);
li = data.gene.li(1:ind);
betap = data.gene.betap(1:ind);

%vmes  = data.exp.vloop(1:ind);
%belical  = smooth(data.gene.betap(1:ind)+data.gene.li(1:ind)/2,2);
belimes  = data.exp.li(1:ind)/2+data.exp.betadia(1:ind);
jlh = data.source.hyb.j/1e6;
jb = data.source.jboot/1e6;
jtot = data.prof.jmoy/1e6;
x = param.gene.x;
%af = post.polar.af.*180./pi;
%afmes = post.polar.afmes.*180./pi;
%ind2 = find(isfinite( prod(af' .* afmes')'));
%daf = mean( af(ind2,:) -afmes(ind2,:));
%afcal = af(1:ind,:)-ones(ind,1)*daf;
%afmes = post.polar.afmes(1:ind,:).*180./pi;
%gou=cat(2,t,ip,nbar,pfci,iboot,te,vcal,vmes,belical,belimes);
%goufara=cat(2,t,afcal,afmes); 
%save gou gou -ASCII
%save goufara goufara -ASCII
%gou=cat(2,t,ip,nbar,pfci,iboot,te,vmes,belimes);
%save gou gou -ASCII
