function f=elliE(k)
% ELLIE		Complete elliptic integral of the second kind E(k)
% 0<k<1, k=m^2  :  eps(m) < 2e-8
% Source: Abramowitz 17.3.34 (page 591)
% J-P Roubin
a1=0.44325141463;
a2=0.06260601220;
a3=0.04757383546;
a4=0.01736506451;
b1=0.24998368310;
b2=0.09200180037;
b3=0.04069697526;
b4=0.00526449639;
y=1-k.*k;
a=(((a4*y+a3).*y+a2).*y+a1).*y+1.;
b=(((b4*y+b3).*y+b2).*y+b1).*y;
f=a+b.*log(1../y);
