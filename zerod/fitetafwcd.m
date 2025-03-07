% fit des mesure d'efficacite de generation de courant FWCD
clear te0 zeff eta
k = 0;

k = k+1;te0(k) =6.5;zeff(k) = 4;eta(k) = 0.030;
k = k+1;te0(k) =6.5;zeff(k) = 2;eta(k) = 0.051;
k = k+1;te0(k) =6.5;zeff(k) = 3;eta(k) = 0.042;
k = k+1;te0(k) =8.25;zeff(k) = 4;eta(k) = 0.037;
k = k+1;te0(k) =8.25;zeff(k) = 2;eta(k) = 0.063;
k = k+1;te0(k) =8.25;zeff(k) = 3;eta(k) = 0.050;
k = k+1;te0(k) =10.0;zeff(k) = 4;eta(k) = 0.044;
k = k+1;te0(k) =10.0;zeff(k) = 2;eta(k) = 0.075;
k = k+1;te0(k) =10.0;zeff(k) = 3;eta(k) = 0.061;
k = k+1;te0(k) =5.5;zeff(k) = 3;eta(k) = 0.045;
k = k+1;te0(k) =6;zeff(k) = 3;eta(k) = 0.045;
k = k+1;te0(k) =5.5;zeff(k) = 3;eta(k) = 0.04;
k = k+1;te0(k) =5.2;zeff(k) = 3;eta(k) = 0.033;
k = k+1;te0(k) =5;zeff(k) = 3;eta(k) = 0.033;
k = k+1;te0(k) =5;zeff(k) = 3;eta(k) = 0.035;
k = k+1;te0(k) =4.5;zeff(k) = 3;eta(k) = 0.029;
k = k+1;te0(k) =4;zeff(k) = 3;eta(k) = 0.022;
k = k+1;te0(k) =3;zeff(k) = 3;eta(k) = 0.02;
k = k+1;te0(k) =2;zeff(k) = 3;eta(k) = 0.013;
k = k+1;te0(k) =1.5;zeff(k) = 3;eta(k) = 0.007;

% ajout pour le fit avec saturation
k = k+1;te0(k) =0;zeff(k) = 1;eta(k) =0;
k = k+1;te0(k) =0;zeff(k) = 2;eta(k) =0;
k = k+1;te0(k) =0;zeff(k) = 3;eta(k) =0;
k = k+1;te0(k) =0;zeff(k) = 4;eta(k) =0;


eta0 = (5+zeff) ./ 6 .* eta;
[pp,ss] = polyfit(te0,eta0,1)
tex = linspace(0,10,21);
[etax,errx] = polyval(pp,tex,ss);
figure(51);clf
plot(te0,eta,'o',te0,eta0,'*',tex,polyval(pp,tex),'r');
hold on
errorbar(tex,etax,errx,'r');
