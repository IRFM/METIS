%  NOM DE LA FONCTION  courte description  
%------------------------------------------------------------------------------- 
% fichier : nom_du_fichier -> nom_fonction_principale, nom_fonction_locale 1,... 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function .... 
%  
% entrees :  
%  
%  
% sorties :  
%  
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  4.1  du  08/04/2008  
%  
%*@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
a1 = load('iter_ramp_up_lh_5MW_37GHz_2npar_Dir75.mat');
a2 = load('iter_ramp_up_lh_5MW_ECRH.mat');
a3 = load('iter_ramp_up_lh_5MW_ECCD.mat');
figure
plot(a1.post.zerod.temps,a1.post.zerod.vloop,'r')
hold on
plot(a2.post.zerod.temps,a2.post.zerod.vloop,'b')
plot(a3.post.zerod.temps,a3.post.zerod.vloop,'c')
legend('LHCD 5MW','ECRH (@0) 5MW','ECCD (@0.6) 5MW');
ylabel('Vloop (V)');
xlabel('time (s)');
set(gca,'xlim',[0,120]);
edition2
vloop = a1.post.zerod.vloop;
vloop(find(~isfinite(vloop))) = sqrt(-1);
flux = cumtrapz(t,vloop,1);
flux(find(imag(vloop))) = NaN;
flux1 = real(flux);
vloop = a2.post.zerod.vloop;
vloop(find(~isfinite(vloop))) = sqrt(-1);
flux = cumtrapz(t,vloop,1);
flux(find(imag(vloop))) = NaN;
flux2 = real(flux);
vloop = a3.post.zerod.vloop;
vloop(find(~isfinite(vloop))) = sqrt(-1);
flux = cumtrapz(t,vloop,1);
flux(find(imag(vloop))) = NaN;
flux3 = real(flux);
figure
plot(a1.post.zerod.temps,flux1,'r')
hold on
plot(a2.post.zerod.temps,flux2,'b')
plot(a3.post.zerod.temps,flux3,'c')
legend('LHCD 5MW','ECRH (@0) 5MW','ECCD (@0.6) 5MW');
ylabel('Flux (V s)');
xlabel('time (s)');
set(gca,'xlim',[0,120]);
edition2
figure
zs = a1.post.zerod
plot(t,zs.etalh0./1e19,'m',t,zs.etalh1./1e19,'c',t,zs.etalh./1e19,'r');
ylabel('eta_L_H (1e19 A W^ -1 m^-^2)');
xlabel('time (s)');
legend('eta_L_H @ Vloop = 0','hot conductivity','eta_L_H');
set(gca,'xlim',[0,120]);
edition2
