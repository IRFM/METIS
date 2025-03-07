% function to compare Nevins and Tentori + validation of rate
function test_cross_section_pB11

x = logspace(3,7,1e4);

% cross section
[sigma_old,S1_old,S2_old,S3_old,S_old] = pb11_cross_section_nevins(x);
[sigma_new,S1_new,S2_new,S3_new,S_new] = pb11_cross_section_tentori(x);

% data extract from figure of Sikora, M.H., Weller, H.R. A New Evaluation of the 
%B(p,3*alpha) Reaction Rates. 
%J Fusion Energ 35, 538543 (2016). https://doi.org/10.1007/s10894-016-0069-y
%xsw      = [0.15    0.3     0.4     0.5     0.640       1       1.4     1.7     2.55    2.85    3.6] .* 1e6;
%sigma_sw = [10      150     400     1000    1400        300     350     250     600     350     600] .* 1e-31;

% reaction rate
[rate_old,rate_NRL_old,rate_NRH_old,rate_R_old] = pb11_fusion_rate_nevins(x);
[rate_new,rate_NRL_new,rate_NRH_new,rate_R_new] = pb11_fusion_rate_tentori(x);

% recompute rate
rate_recomputed_old = compute_rate(x,sigma_old);
rate_recomputed_new = compute_rate(x,sigma_new);


figure;

subplot(2,2,1)
plot(x/1e3,S1_old/1e6,x/1e3,S1_new/1e6)
set(gca,'xlim',[0,400]);
xlabel('E (keV)');
ylabel('S1 (MeV b)');

subplot(2,2,2)
plot(x/1e3,S2_old/1e6,x/1e3,S2_new/1e6)
set(gca,'xlim',[400 700]);
xlabel('E (keV)');
ylabel('S2 (MeV b)');

subplot(2,2,3)
plot(x/1e3,S3_old/1e6,x/1e3,S3_new/1e6)
set(gca,'xlim',[600 5000]);
xlabel('E (keV)');
ylabel('S3 (MeV b)');

subplot(2,2,4)
loglog(x/1e3,S_old/1e6,x/1e3,S_new/1e6)
set(gca,'xlim',[100 10000]);
xlabel('E (keV)');
ylabel('S (MeV b)');
legend('Nevins','Tentori');


figure;
%loglog(x/1e3,sigma_old * 1e28,x/1e3,sigma_new *1e28,xsw /1e3,sigma_sw*1e28,'o');
loglog(x/1e3,sigma_old * 1e28,x/1e3,sigma_new *1e28);
set(gca,'xlim',[100 10000]);
xlabel('E (keV)');
ylabel('sigma (b)');
%legend('Nevins','Tentori','Sikora');
legend('Nevins','Tentori');

figure;
subplot(2,2,1)
plot(x/1e3,rate_R_old,x/1e3,rate_R_new)
set(gca,'xlim',[0 1000]);
xlabel('E (keV)');
ylabel('Rate_R (m^3 / s)');

subplot(2,2,2)
plot(x/1e3,rate_NRL_old,x/1e3,rate_NRL_new)
set(gca,'xlim',[0 70]);
xlabel('E (keV)');
ylabel('Rate_{NR,low} (m^3 / s)');

subplot(2,2,3)
plot(x/1e3,rate_NRH_old,x/1e3,rate_NRH_new)
set(gca,'xlim',[50 600]);
xlabel('E (keV)');
ylabel('Rate_{NR,high} (m^3 / s)');

subplot(2,2,4)
%xsw = [0  120   200     300     400     500     600     700     800     900] .* 1e3;
%rsw = [0  1     3       4.45    5.2     5.6     5.9     6.2      6.3    6.4] .* 1e-22;   
semilogx(x/1e3,rate_old,x/1e3,rate_new, ...
    x/1e3,rate_recomputed_old,x/1e3,rate_recomputed_new,[500 500],[0 8e-22]);%,xsw/1e3,rsw,'o');
set(gca,'xlim',[10 1000]);
xlabel('E (keV)');
ylabel('Rate (m^3 / s)');
legend('Nevins','Tentori','Recomputed from sigma Nevins','Recomputed from sigma Tentori'); 




function rate = compute_rate(x,sigma)



% physical constants
phys = cphys;

% reduce masse in eV
% for pure B11 and proton
mu = 11.0093054 * 1.00782503207 / (11.0093054 + 1.00782503207) * phys.ua * phys.c ^ 2  / phys.e;
% for natural Boron and natural hydrogen 
%mu = 859.526e6;

ss = size(x);
x = x(:);
sigma = ones(size(x)) * sigma(:)';
xx = x * ones(1,size(sigma,2));
E  = xx';
sigma(~isfinite(sigma)) = 0;

rate = phys.c .* sqrt(8/pi/mu) ./ x .^ (3/2) .* trapz(E(1,:),E .* sigma .* exp(-E ./ xx),2);
rate = reshape(rate,ss);

