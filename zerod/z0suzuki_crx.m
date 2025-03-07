% Attenuation of high-energy neutral hydrogen beams in high-density plasmas
% S Suzuki†§, T Shirai‡, M Nemoto†, K Tobita†, H Kubo†, T Sugie†, A Sakasai† and Y Kusama†
function sv = z0suzuki_crx(A,E,ne,te,nhe,nlow,nhigh,nw,zlow,zhigh,Sn_fraction)

% changement d'unite
E      = E ./ A ./ 1e3;
epsi   = log(E);
N      = ne ./ 1e19;
te     = te ./ 1e3;
U      = log(te);

%'He','Li','Be','B','C','N','O','Fe'
zimp = [2,3,4,5,6,7,8,26];
sz_w = sz_one(8,E,epsi,N,U); % in fact Fe as it is the heavier element avalailable in the table
sz_he = sz_one(1,E,epsi,N,U);

if zlow < 8
    sz_low = sz_one(min(8,max(1,ceil(zlow) - 1)),E,epsi,N,U);
else
    fz     = max(0,min(1,(zlow - 26) ./ (8 - 26)));
    sz_low = fz .* sz_one(7,E,epsi,N,U) + fz .* sz_one(8,E,epsi,N,U);
end
if zhigh < 8
    sz_high = sz_one(min(8,max(1,ceil(zhigh) - 1)),E,epsi,N,U);
else
    fz     = max(0,min(1,(zhigh - 26) ./ (8 - 26)));
    sz_high = fz .* sz_one(7,E,epsi,N,U) + fz .* sz_one(8,E,epsi,N,U);
end
% compute equivalent factor with sn
factor_w_sn =  (1 - Sn_fraction) .* nw .* 74 .* 73 + Sn_fraction .* nw .* 50 .* 49;

% Janev R K, Boley C D and Post D E 1989 Nucl. Fusion 12 2125
szsum  = (nhe .* 2 .* 1 .* sz_he + factor_w_sn .* sz_w + nlow .* zlow .* (zlow - 1) .* sz_low + nhigh .* zhigh .* (zhigh - 1) .* sz_high) ./ ne;
% cm ^ 2 - > m^2
sv = 1e-4 .* sigma_h(A,E,epsi,N,U) .* (1 + szsum);


function sv = sigma_h(A,E,epsi,N,U)

% security
A = max(1,min(3,A));

% low energy part
if E < 100
    E = max(10,E);
    if A < 2
       sv = (2 - A) .* sigma_h_tab(3,1,E,epsi,N,U) + (A - 1) .* sigma_h_tab(3,2,E,epsi,N,U);
    else
       sv = (3 - A) .* sigma_h_tab(3,2,E,epsi,N,U) + (A - 2) .* sigma_h_tab(3,3,E,epsi,N,U);   
    end
else
% high energy part
    E = min(1e4,E);
    if A < 2
       sv = (2 - A) .* sigma_h_tab(2,1,E,epsi,N,U) + (A - 1) .* sigma_h_tab(2,2,E,epsi,N,U);
    else
       sv = (3 - A) .* sigma_h_tab(2,2,E,epsi,N,U) + (A - 2) .* sigma_h_tab(2,3,E,epsi,N,U);    
    end
end

function sv = sigma_h_tab(tab,col,E,epsi,N,U)

switch tab
case 2
    d = table2_data;
case 3
    d = table3_data;
otherwise
    error('undefined table');
end
% in cm^2
sv = d.A1(col) .* 1e-16 ./ E .* (1 + d.A2(col) .* epsi + d.A3(col) .* epsi .^ 2) .* ...
     (1 + (1- exp(-d.A4(col) * N)) .^ d.A5(col) .* (d.A6(col) + d.A7(col) .* epsi + d.A8(col) .* epsi .^ 2)) .* ...
     (1 + d.A9(col) .* U + d.A10(col) .* U .^ 2);
     

function sz = sz_one(col,E,epsi,N,U)

% low energy part
if E < 100
    d = table3_data;
else
% high energy part
    d = table2_data;
end

n  = log(N);
sz = d.B111(col) .* epsi .^ 0  .* n .^ 0  .* U .^ 0 + ...
     d.B112(col) .* epsi .^ 0  .* n .^ 0  .* U .^ 1 + ...
     d.B121(col) .* epsi .^ 0  .* n .^ 1  .* U .^ 0 + ...
     d.B122(col) .* epsi .^ 0  .* n .^ 1  .* U .^ 1 + ...
     d.B211(col) .* epsi .^ 1  .* n .^ 0  .* U .^ 0 + ...
     d.B212(col) .* epsi .^ 1  .* n .^ 0  .* U .^ 1 + ...
     d.B221(col) .* epsi .^ 1  .* n .^ 1  .* U .^ 0 + ...
     d.B222(col) .* epsi .^ 1  .* n .^ 1  .* U .^ 1 + ...
     d.B311(col) .* epsi .^ 2  .* n .^ 0  .* U .^ 0 + ...
     d.B312(col) .* epsi .^ 2  .* n .^ 0  .* U .^ 1 + ...
     d.B321(col) .* epsi .^ 2  .* n .^ 1  .* U .^ 0 + ...
     d.B322(col) .* epsi .^ 2  .* n .^ 1  .* U .^ 1;
     
     

function table1 = table1_data
% table 1
table1.Aijk = {'H','D','T','D@1T'};
table1.A111 = [4.28,4.27,4.27,4.32]; 
table1.A112 = [-3.91e-2,-3.22e-2,-3.09e-2,-3.02e-2];
table1.A121 = [9.43e-2,9.41e-2,9.39e-2,7.98e-2];
table1.A122 = [-4.07e-3,-3.47e-3,-3.17e-3,-3.62e-3];
table1.A131 = [6.40e-3,6.40e-3,6.39e-3,7.35e-3];
table1.A132 = [1.56e-4,1.99e-4,2.28e-4,1.50e-4];
table1.A211 = [2.51e-1,2.52e-1,2.53e-1,2.42e-1];
table1.A212 = [-7.63e-3,-8.46e-3,-8.61e-3,-8.85e-3];
table1.A221 = [-6.57e-3,-6.55e-3,-6.52e-3,-3.51e-3];
table1.A222 = [2.37e-4,1.63e-4,1.25e-4,1.98e-4];
table1.A231 = [9.84e-4,9.84e-4,9.86e-4,7.50e-4];
table1.A232 = [-4.18e-5,-4.68e-5,-5.05e-5,-3.86e-5];


table1.Bijk = {'He','Li','Be','B','C','N','O','Fe'};

table1.B111 = [-1.05,-1.27,-1.28,-1.32,-1.54,-1.50,-1.46,-4.27e-1]; 
table1.B112 = [1.41e-1,-1.41e-3,-4.39e-2,-6.00e-2,-8.68e-2,-8.83e-2,-8.50e-2, 4.39e-2];
table1.B121 = [-3.75e-1,-2.84e-1,-2.55e-1,-2.25e-1,-1.83e-1,-1.47e-1,-1.05e-1, 1.03e-1];
table1.B122 = [-1.55e-2,-1.84e-2,-1.76e-2,-1.85e-2,-1.53e-2,-1.19e-2,-8.88e-3, 1.24e-2];
table1.B211 = [5.31e-1,5.52e-1,5.29e-1,5.19e-1,5.67e-1,5.37e-1,5.10e-1,8.27e-2];
table1.B212 = [-3.09e-2,1.02e-2,2.09e-2,2.42e-2,3.13e-2,3.00e-2,2.78e-2,2.24e-2];
table1.B221 = [1.05e-1,7.81e-2,6.82e-2,5.94e-2,4.60e-2,3.55e-2,2.37e-2,3.78e-2];
table1.B222 = [5.03e-3,5.63e-3,5.20e-3,5.32e-3,4.21e-3,3.12e-3,2.21e-3,4.72e-3];
table1.B311 = [-4.17e-2,-4.08e-2,-3.77e-2,-3.59e-2,-3.86e-2,-3.54e-2,-3.29e-2, 2.84e-3];
table1.B312 = [2.58e-3,-2.95e-4,-9.73e-4,-1.14e-3,-1.60e-3,-1.41e-3,-1.20e-3, 2.91e-3];
table1.B321 = [-7.02e-3,-5.09e-3,-4.27e-3,-3.67e-3,-2.68e-3,-1.97e-3,-1.18e-3, 3.19e-3];
table1.B322 = [-3.47e-4,-3.75e-4,-3.27e-4,-3.29e-4,-2.41e-4,-1.61e-4,-9.79e-5, 4.32e-4];


% table 1 (continued)
table1.Bijk_continued = {'C','N','O','Fe'};


table1.B111_continued = [-1.59,-1.60,-1.59,-1.30];
table1.B112_continued = [-7.52e-2,-6.98e-2,-6.18e-2,1.41e-2];
table1.B121_continued = [-1.88e-1,-1.58e-1,-1.24e-1,1.83e-2];
table1.B122_continued = [-1.57e-2,-1.29e-2,-1.07e-2,1.99e-3];
table1.B211_continued = [5.85e-1,5.68e-1,5.51e-1,3.64e-1];
table1.B212_continued = [2.77e-2,2.44e-2,2.07e-2,-1.10e-2];
table1.B221_continued = [4.80e-2,3.97e-2,3.01e-2,-9.37e-3];
table1.B222_continued = [4.37e-3,3.52e-3,2.87e-3,-1.03e-3];
table1.B311_continued = [-4.00e-2,-3.77e-2,-3.60e-2,-1.89e-2];
table1.B312_continued = [-1.35e-3,-1.01e-3,-7.01e-4,1.88e-3];
table1.B321_continued = [-2.84e-3,-2.31e-3,-1.69e-3,9.32e-4];
table1.B322_continued = [-2.57e-4,-1.99e-4,-1.57e-4,1.19e-4];

function table2 = table2_data
% table 2

table2.Ai = {'H','D','T','D@1T'};

table2.A1  = [1.27e1,1.41e1,1.27e1,-2.47e1];
table2.A2  = [1.25,1.11,1.26,-1.30];
table2.A3  = [4.52e-1,4.08e-1,4.49e-1,-1.54e-1];
table2.A4  = [1.05e-2,1.05e-2,1.05e-2,8.02e-3];
table2.A5  = [5.47e-1,5.47e-1,5.47e-1,4.81e-1];
table2.A6  = [-1.02e-1,-4.03e-2,-5.77e-3,-1.49e-1];
table2.A7  = [3.60e-1,3.45e-1,3.36e-1,3.92e-1];
table2.A8  = [-2.98e-2,-2.88e-2,-2.82e-2,-2.99e-2];
table2.A9  = [-9.59e-2,-9.71e-2,-9.74e-2,-9.76e-2];
table2.A10 = [4.21e-3,4.74e-3,4.87e-3,4.79e-3];


table2.Bijk = {'He','Li','Be','B','C','N','O','Fe'};

table2.B111 = [2.31e-1,-4.41e-1,-6.13e-1,-7.32e-1,-1.00,-1.00,-9.89e-1,-1.93e-1];
table2.B112 = [3.43e-1,1.29e-1,5.52e-2,1.83e-2,-2.55e-2,-4.15e-2,-4.98e-2,5.03e-3];
table2.B121 = [-1.85e-1,-1.70e-1,-1.67e-1,-1.55e-1,-1.25e-1,-9.76e-2,-6.36e-2,1.06e-1];
table2.B122 = [-1.62e-2,-1.62e-2,-1.59e-2,-1.72e-2,-1.42e-2,-1.10e-2,-7.75e-3,1.25e-2];
table2.B211 = [1.05e-1,2.77e-1,3.04e-1,3.21e-1,3.88e-1,3.70e-1,3.52e-1,2.36e-3];
table2.B212 = [-7.03e-2,-1.56e-2,1.54e-3,9.46e-3,2.06e-2,2.28e-2,2.34e-2,-7.41e-3];
table2.B221 = [5.31e-2,4.66e-2,4.36e-2,3.97e-2,2.97e-2,2.15e-2,1.17e-2,-3.88e-2];
table2.B222 = [3.42e-3,3.79e-3,3.78e-3,4.20e-3,3.26e-3,2.33e-3,1.42e-3,-4.60e-3];
table2.B311 = [-8.38e-3,-1.93e-2,-2.01e-2,-2.04e-2,-2.46e-2,-2.22e-2,-2.04e-2,9.39e-3];
table2.B312 = [4.15e-3,7.53e-4,-2.16e-4,-6.19e-4,-1.31e-3,-1.33e-3,-1.29e-3,1.58e-3];
table2.B321 = [-3.35e-3,-2.86e-3,-2.51e-3,-2.24e-3,-1.48e-3,-9.25e-4,-2.75e-4,3.32e-3];
table2.B322 = [-2.21e-4,-2.39e-4,-2.27e-4,-2.54e-4,-1.80e-4,-1.15e-4,-5.45e-5,3.88e-4];

% table 2 (continued)
table2.Bijk_continued = {'C','N','O','Fe'};


table2.B111_continued = [-1.01,-1.02,-1.02,-8.20e-1];
table2.B112_continued = [-8.65e-3,-1.39e-2,-1.48e-2,-6.36e-3];
table2.B121_continued = [-1.24e-1,-9.79e-2,-6.74e-2,5.42e-2];
table2.B122_continued = [-1.45e-2,-1.17e-2,-9.17e-3,3.95e-3];
table2.B211_continued = [3.91e-1,3.75e-1,3.59e-1,2.02e-1];
table2.B212_continued = [1.61e-2,1.56e-2,1.43e-2,8.06e-4];
table2.B221_continued = [2.98e-2,2.24e-2,1.39e-2,-2.00e-2];
table2.B222_continued = [3.32e-3,2.54e-3,1.84e-3,-1.78e-3];
table2.B311_continued = [-2.48e-2,-2.26e-2,-2.09e-2,-6.10e-3];
table2.B312_continued = [-1.04e-3,-8.89e-4,-7.32e-4,6.51e-4];
table2.B321_continued = [-1.52e-3,-1.04e-3,-5.02e-4,1.75e-3];
table2.B322_continued = [-1.89e-4,-1.39e-4,-9.49e-5,1.46e-4];



function table3 = table3_data
% table 3

table3.Ai = {'H','D','T','D@1T'};

table3.A1  = [-5.29e1,-6.79e1,-7.42e1,-6.98e1];
table3.A2  = [-1.36,-1.22,-1.18,-1.21];
table3.A3  = [7.19e-2,8.14e-2,8.43e-2,8.33e-2];
table3.A4  = [1.37e-2,1.39e-2,1.39e-2,1.35e-2];
table3.A5  = [4.54e-1,4.54e-1,4.53e-1,4.45e-1];
table3.A6  = [4.03e-1,4.65e-1,4.91e-1,4.89e-1];
table3.A7  = [-2.20e-1,-2.73e-1,-2.94e-1,-2.90e-1];
table3.A8  = [6.66e-2,7.51e-2,7.88e-2,7.86e-2];
table3.A9  = [-6.77e-2,-6.30e-2,-6.12e-2,-6.30e-2];
table3.A10 = [-1.48e-3,-5.08e-4,-1.85e-4,-5.12e-4];


table3.Bijk = {'He','Li','Be','B','C','N','O','Fe'};

table3.B111 = [-7.92e-1,1.12e-1,1.12e-1,1.22e-1,1.58e-1,1.34e-1,1.07e-1,-4.65e-5];
table3.B112 = [4.20e-2,4.95e-2,4.95e-2,5.27e-2,5.54e-2,5.24e-2,4.31e-2,-7.29e-4];
table3.B121 = [5.30e-2,1.16e-2,1.16e-2,-4.30e-4,-4.31e-3,-4.69e-3,-2.83e-3,-3.10e-3];
table3.B122 = [-1.39e-2,-2.86e-3,-2.86e-3,-3.18e-3,-3.35e-3,-3.06e-3,-1.82e-3,-5.17e-4];
table3.B211 = [3.01e-1,-1.49e-1,-1.49e-1,-1.51e-1,-1.55e-1,-1.30e-1,-1.06e-1,-1.34e-2];
table3.B212 = [-2.64e-2,-3.31e-2,-3.31e-2,-3.64e-2,-3.74e-2,-3.52e-2,-2.89e-2,-5.06e-5];
table3.B221 = [-2.99e-2,-4.26e-3,-4.26e-3,3.43e-3,5.37e-3,5.31e-3,3.87e-3,2.87e-3];
table3.B222 = [6.07e-3,9.80e-4,9.80e-4,1.51e-3,1.74e-3,1.64e-3,9.10e-4,2.74e-4];
table3.B311 = [2.72e-4,4.47e-2,4.47e-2,4.20e-2,3.88e-2,3.31e-2,2.77e-2,5.04e-3];
table3.B312 = [6.11e-3,6.52e-3,6.52e-3,6.92e-3,6.83e-3,6.35e-3,5.27e-3,2.13e-4];
table3.B321 = [3.47e-3,-3.56e-4,-3.56e-4,-1.41e-3,-1.60e-3,-1.51e-3,-1.23e-3,-6.65e-4];
table3.B322 = [-9.19e-4,-2.03e-4,-2.03e-4,-2.90e-4,-3.22e-4,-3.04e-4,-1.89e-4,-5.83e-5];


% table 3 (continued)
table3.Bijk_continued = {'C','N','O','Fe'};


table3.B111_continued = [1.61e-1,1.39e-1,1.11e-1,-1.10e-2];
table3.B112_continued = [5.98e-2,6.06e-2,5.41e-2,2.02e-2];
table3.B121_continued = [-3.36e-3,-3.06e-3,-3.46e-4,9.46e-4];
table3.B122_continued = [-4.26e-3,-4.55e-3,-3.68e-3,-4.09e-3];
table3.B211_continued = [-1.57e-1,-1.33e-1,-1.08e-1,-6.66e-3];
table3.B212_continued = [-3.96e-2,-3.94e-2,-3.47e-2,-1.17e-2];
table3.B221_continued = [4.60e-3,3.99e-3,1.93e-3,-2.36e-4];
table3.B222_continued = [2.19e-3,2.36e-3,1.81e-3,2.02e-3];
table3.B311_continued = [3.91e-2,3.35e-2,2.80e-2,4.08e-3];
table3.B312_continued = [7.11e-3,6.90e-3,6.04e-3,1.85e-3];
table3.B321_continued = [-1.44e-3,-1.24e-3,-8.41e-4,-6.48e-5];
table3.B322_continued = [-3.85e-4,-4.05e-4,-3.17e-4,-3.13e-4];








