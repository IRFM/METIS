% ref : P.B. Snyder PoP 2009 vol 16  p  056118-
% ref : Onjun PoP 2002 vol 9 p -5022-
% model reference :F.D. Halpern et al, PoP 15 (2008) p 062505
function [wped_onjun,pped_onjun] = z0pped_onjun(zs,profli,geo,cons,option)

if isfield(profli,'kx')
	k95 = profli.kx(:,end - 1);
else
	k95 = 0.95 .* geo.K; 
end
if isfield(profli,'dx')
	d95 = profli.dx(:,end - 1);
else
	d95 = 0.95 .* geo.d; 
end
%
vt = ones(size(zs.ip,1),1);
x  = [0.94 0.95 0.96];
ve = ones(size(x));
%
a    = geo.a;
R    = geo.R;
Bt   = geo.b0 ;
ip   = (zs.ip) / 1e6;
if isfield(profli,'nep')
	nped = profli.nep(:,end - 1);
else
	nped = zs.nem;
end 

% il est difficile de connaitre q et s dans le piedestal de METIS
ee   =   1.602176462e-19;          % electron charge (C)   (+/- 0.000000063e-19)
mu0  =   4*pi*1e-7;                % permeablity of vacuum (H/m) (definition)
% utilisation d'un formule analytique :
gsan   = ((1 + k95 .^2 .* (1 + 2 .* d95 .^ 2 - 1.2 .* d95 .^ 3)) .* (1.17 - 0.65 .* geo.a ./ geo.R)) ./  ...
         (1 - (geo.a ./ geo.R) .^ 2) .^ 2 ./ 2;
qcyl   = 5  .* geo.a .^ 2 .* geo.b0 ./ zs.ip ./ mu0 ./ geo.R;
%qrad   = ((1 + (0.95 ./ 1.4 ) .^ 2) .^ 2 + 0.27 .* abs(log(1 - 0.95))) ./ 2.5;
qan    = qcyl  .* gsan; 
san    = 2.68 .* ones(size(qan)); 
if isfield(profli,'qjli')
	q  = pchip(profli.xli,profli.qjli,x);
	gs = profli.qjli(:,end-1) ./ qcyl;
	x    = vt*x;
	s    = x(:,2) ./ q(:,2) .* (q(:,3) - q(:,1)) ./ (x(:,3) - x(:,1)); 
	q    = q(:,2);
%  figure(21);clf;
%  subplot(2,2,1);plot(zs.temps,s,'b',zs.temps,san,'r');
%  subplot(2,2,2);plot(zs.temps,q,'b',zs.temps,qan,'r');
%  subplot(2,2,3);plot(zs.temps,gs,'b',zs.temps,gsan,'r');
%  drawnow
else
	q    = qan;
	gs   = gsan;
	s    = san;
end
ac   = 0.4 .* s .* (1 + k95 .^ 2 .* (1 + 5 .* d95 .^2));
c6   = 0.025;
%
tped = c6 .^ 2 ./ (4 .* mu0 .* ee) .* (geo.b0 ./ q .^ 2) .^ 2 .* (geo.R ./ geo.a) .^ 2 .* (ac .^ 2 ./ nped)  .* (pi .* q .* (1 + k95) ./ 5 ./ gs) .^ 2;
tped = max(zs.tebord,min(option.te_max,tped .* zs.modeh)); 
if isfield(profli,'nip') && isfield(profli,'vpr')
	tpped      =  tped * ones(size(profli.xli));
	tpped(:,end) = zs.tebord;
	wped_onjun = (3./2) .* 1.602176462e-19 .* trapz(profli.xli,(profli.nep + profli.nip) .*  tpped .* profli.vpr,2);
else
	wped_onjun = (3./2) .* zs.vp .* (zs.nem + zs.nim) .* 1.602176462e-19 .* tped;
end 

pped_onjun=2.*(tped.*nped).*ee;

