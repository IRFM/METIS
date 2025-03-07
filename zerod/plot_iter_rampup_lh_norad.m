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
a1 = load('iter_ramp_up_lh_5MW_37GHz_LH_norad.mat');
t= a1.post.zerod.temps;
a2 = load('iter_ramp_up_lh_OH_norad.mat');
figure
plot(a1.post.zerod.temps,a1.post.zerod.vloop,'r')
hold on
plot(a2.post.zerod.temps,a2.post.zerod.vloop,'b')
legend('LH 5MW','OH startup');
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
figure
plot(a1.post.zerod.temps,flux1,'r')
hold on
plot(a2.post.zerod.temps,flux2,'b')
legend('LH 5MW','OH startup');
ylabel('Flux (Vs)');
xlabel('time (s)');
set(gca,'xlim',[0,120]);
edition2
figure
zs = a1.post.zerod;
plot(t,zs.etalh0./1e19,'m',t,zs.etalh1./1e19,'c',t,zs.etalh./1e19,'r');
ylabel('eta_L_H (1e19 A W^ -1 m^-^2)');
xlabel('time (s)');
legend('eta_L_H @ Vloop = 0','hot conductivity','eta_L_H');
set(gca,'xlim',[0,120]);
edition2


figure
zplotprof(gca,a1.post.profil0d.temps,a1.post.profil0d.xli,a1.post.profil0d.qjli,'color','r')
zplotprof(gca,a2.post.profil0d.temps,a2.post.profil0d.xli,a2.post.profil0d.qjli,'color','b')
legend('LH 5MW','OH startup');
ylabel('q ');
xlabel('x (su)');
