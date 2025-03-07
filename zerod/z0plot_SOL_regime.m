% plot the SOL collisionallity regime 
% references:
%  D. Carralero et al, Nucl. Fus. 54 (2014) 123005
%  J. R Myra and D. A. D'Ippolito, PoP 13 (2006) 092509
%  J. R. Myra et al, PoP 13 (2006) 112502
%  D. A. Russel, PoP 14 (2007) 102307

%
% these E. Tsitrone
lclim = pi .* post.z0dinput.geo.R .* post.zerod.qa;
lcpol = pi .* post.z0dinput.geo.R;
lcx = sqrt(post.zerod.peri .^ 2  + (pi .* post.z0dinput.geo.R .* post.z0dinput.option.lcx .* post.zerod.qa) .^ 2);  
switch post.z0dinput.option.configuration
case 0
	lc = lcpol;
case 1
	lc = lclim;
case 2
	lc  = post.zerod.xpoint .* lcx + (~post.zerod.xpoint) .* lcpol;
case 3
	lc  = post.zerod.xpoint .* lcx + (~post.zerod.xpoint) .* lclim;
otherwise
	lc  = lcx;
end
lc   = interp1(post.zerod.temps,lc,post.profil0d.temps,'linear','extrap');
dsol = interp1(post.zerod.temps,post.zerod.dsol,post.profil0d.temps,'linear','extrap');
a    = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'linear','extrap');
rext = post.profil0d.Raxe(:,end) + a; 


lambda = lc ./ cs(:,end) ./ tei(:,end) .* wci(:,end) ./ wce(:,end);
rhos_i = cs ./ wci;
a_hat  = dsol .* rext .^ (1/5) ./ lc .^ (2/5) ./ rhos_i(:,end) .^ (4/5);
theta  = a_hat .^(5/2);


ind_RB = double((lambda > theta) & (lambda > 1));
ind_RX = double(((lambda < theta) & (lambda > 1)) | ((lambda > (0.5 .* theta)) & (lambda > 1)));
ind_CS = double((theta > 2) & (lambda <= 1));
ind_CI = double((ind_RB == 0) & (ind_RX == 0) & (ind_CS == 0));

ind_RB(ind_RB==0) = NaN;
ind_RX(ind_RX==0) = NaN;
ind_CS(ind_CS==0) = NaN;
ind_CI(ind_CI==0) = NaN;



hz =findobj(0,'type','figure','tag','SOL_regime');
if isempty(hz)
  	  hz=figure('tag','SOL_regime','name','SOL regime');
else
  	  figure(hz);
end
clf
set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
subplot(3,1,1)
plot(post.profil0d.temps,lambda,'r',post.profil0d.temps,ones(size(post.profil0d.temps)),'b');
xlabel('time (s)');
ylabel('Lambda');
title(sprintf('METIS : %s@%d/SOL regime', ...
             post.z0dinput.machine,post.z0dinput.shot));
legend('Lambda','Lamda=1 limit')
z0loglin(gca);
subplot(3,1,2);
semilogy(post.profil0d.temps,lambda,'r',post.profil0d.temps,theta,'b');
xlabel('time (s)');
legend('Lambda','Theta')
z0loglin(gca);
subplot(3,1,3);
plot(post.profil0d.temps,ind_RB,'*',post.profil0d.temps,ind_RX,'x',post.profil0d.temps,ind_CS,'o',post.profil0d.temps,ind_CI,'+');
xlabel('time (s)');
set(gca,'ylim',[0,1.5]);
legend('Resistive balooning','Resistive X-point','Sheath interchange','Ideal interchange');
joint_axes(hz,2);
edition2



% alpha SOL and MHD in SOL
% minimal stable lambda_q in SOL
mu0 = 4*pi*1e-7;
ee = 1.602176462e-19;
alpha_sol = 2 * mu0 * (z0dinput.geo.a + z0dinput.geo.R) .^ 3 .* post.zerod.qeff .^ 2 ./ (z0dinput.geo.R .* z0dinput.geo.b0) .^ 2 .* ...
           ee .* (post.zerod.nebord .* post.zerod.tebord + post.zerod.nibord .* post.zerod.tibord) ./ post.zerod.dsol;
shear  = pdederive(post.profil0d.xli,post.profil0d.qjli,0,2,2,1)  .* (ones(size(post.profil0d.qjli,1),1) * post.profil0d.xli) ./ post.profil0d.qjli;
alpha_sol_critical = 0.4 .* shear(:,end) .* (1 + post.profil0d.kx(:,end) .^ 2 .* (1 + 5 .*   post.profil0d.dx(:,end) .^ 2));
alpha_sol_critical = interp1(post.profil0d.temps,alpha_sol_critical,post.zerod.temps,'linear','extrap');
dsol_sol_min = 2 * mu0 * (z0dinput.geo.a + z0dinput.geo.R) .^ 3 .* post.zerod.qeff .^ 2 ./ (z0dinput.geo.R .* z0dinput.geo.b0) .^ 2 .* ...
           ee .* (post.zerod.nebord .* post.zerod.tebord + post.zerod.nibord .* post.zerod.tibord) ./ alpha_sol_critical;
% banana width and larmor radius
switch post.z0dinput.option.gaz
    case 3
        ag = (2 + 3 .* post.z0dinput.cons.iso) ./ (1 + post.z0dinput.cons.iso);
        zg = 1;
        lg = (6.46e-3 + 7.92e-3 .* post.z0dinput.cons.iso) ./ (1 + post.z0dinput.cons.iso);
    case 4
        ag = 4;
        zg = 2;
        lg = 4.55e-3;
    case 2
        ag = 2;
        zg = 1;
        lg = 6.46e-3;
    case 5
        ag = (2 + 3 .* real(post.z0dinput.cons.iso)) ./ (1 + real(post.z0dinput.cons.iso));
        zg = 1;
        lg = 4.576e-3 .* (sqrt(2) ./ 1  + sqrt(3) ./ 2 .* real(post.z0dinput.cons.iso)) ./ (1 + real(post.z0dinput.cons.iso));
        warning('nHe3onD & nTonD not yet implemented !');

    case 11
        ag = (1 + 11 .* post.z0dinput.cons.iso) ./ (1 + post.z0dinput.cons.iso);
        zg = 1;
        lg = 4.576e-3 .* (sqrt(1) ./ 1  + sqrt(11) ./ 5 .* post.z0dinput.cons.iso) ./ (1 + post.z0dinput.cons.iso);       
    otherwise
        ag = 1;
        zg = 1;
        lg = 4.576e-3;
end
rloc = z0dinput.geo.a + z0dinput.geo.R;
ral    = lg .* sqrt(post.zerod.tibord ./ 1e3) ./ (z0dinput.geo.R .* z0dinput.geo.b0) .* rloc; % en m
qp     =  post.zerod.qeff;
dp1    = rloc .* (2 .* qp .* ral ./ rloc) .^ (2/3);
% banana
dp2     = sqrt(post.z0dinput.geo.a  ./ rloc) .* ral .* qp;
dp      = dp2 .* (dp2 < post.z0dinput.geo.a) + dp1 .* (dp2 >= post.z0dinput.geo.a);
%
dsol_sol_min = max(dsol_sol_min,max(dp,ral));
%
hz =findobj(0,'type','figure','tag','alpha_SOL');
if isempty(hz)
  	  hz=figure('tag','alpha_regime','name','SOL MHD');
else
  	  figure(hz);
end
clf
set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
subplot(2,1,1)
plot(post.zerod.temps,alpha_sol,'r',post.zerod.temps,alpha_sol_critical,'b');
xlabel('time (s)');
ylabel('\alpha');
legend('\alpha_{SOL}','\alpha_{critic}');
title(sprintf('METIS : %s@%d/SOL MHD',post.z0dinput.machine,post.z0dinput.shot));
subplot(2,1,2)
plot(post.zerod.temps,post.zerod.dsol,'r',post.zerod.temps,dsol_sol_min,'b');
xlabel('time (s)');
ylabel('\lambda_q (m)');
legend('d_{SOL}','d_{SOl,min}');
z0loglin(gca);
joint_axes(hz,2);
edition2

