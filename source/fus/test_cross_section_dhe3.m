% function to compare Nevins and Tentori + validation of rate
function test_cross_section_dhe3

if isappdata(0,'USED_ZSIGMAVFUSION_INSTEAD_ZSIGMAVFUSION2')
    rmappdata(0,'USED_ZSIGMAVFUSION_INSTEAD_ZSIGMAVFUSION2')
end

x = logspace(3,7,1e4);

% cross section
sigma = z0csdhe3(x/1e3);

% reaction rate
[dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion(x);
sv1 = dhe3.sv;
[dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion2(x);
sv2 = dhe3.sv;

% recompute rate
rate_recomputed = 1e-31 .* compute_rate_dhe3(x,sigma);


figure;

loglog(x/1e3,sigma)
set(gca,'xlim',[0,4800]);
xlabel('E (keV)');
ylabel('sigma (MeV b)');


figure

semilogx(x/1e3,sv1,x/1e3,sv2, x/1e3,rate_recomputed);
set(gca,'xlim',[10 1000]);
xlabel('Ti (keV)');
ylabel('Rate (m^3 / s)');
legend('Bosh & Gale','Mavrin','Recomputed'); 




function rate = compute_rate_dhe3(x,sigma)



% physical constants
phys = cphys;

% reduce masse in eV
% for pure D and He3
mu = 3.016029322 * 2.0141017778 / (3.016029322 + 2.0141017778) * phys.ua * phys.c ^ 2  / phys.e;

ss = size(x);
x = x(:);
sigma = ones(size(x)) * sigma(:)';
xx = x * ones(1,size(sigma,2));
E  = xx';
sigma(~isfinite(sigma)) = 0;

rate = phys.c .* sqrt(8/pi/mu) ./ x .^ (3/2) .* trapz(E(1,:),E .* sigma .* exp(-E ./ xx),2);
rate = reshape(rate,ss);

