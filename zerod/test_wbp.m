% test du calcul de wbp
t= z0dinput.cons.temps;
wbp1  = (4*pi*1e-7) .* z0dinput.exp.ip .^ 2 .* z0dinput.geo.R ./ 4 .* z0dinput.exp.li;
wbp2  = (4*pi*1e-7) .* z0dinput.exp.ip .^ 2 .* z0dinput.geo.R ./ 2 .* (z0dinput.exp.li ./2 + log(8.*z0dinput.geo.R./z0dinput.geo.a) - 2);

rp    = z0dinput.exp.pohm ./ z0dinput.exp.ip .^ 2;
v1    = rp .* z0dinput.exp.ip + zdxdt(wbp1,t)./z0dinput.exp.ip;
v2    = rp .* z0dinput.exp.ip + zdxdt(wbp2,t)./z0dinput.exp.ip; 

figure(20);plot(t,rp .* z0dinput.exp.ip,'r',t,v1,'m',t,v2,'c',t,z0dinput.exp.vloop,'k');
