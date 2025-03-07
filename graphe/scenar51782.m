load /usr/drfc/artaud/zineb/data/zineb51782_20010503.mat
load /usr/drfc/cgc/matlab5/tcron/JET/51782/temp51782
%t  = data.gene.temps;
%pidn = data.gene.paddidn;
%plh  = data.gene.paddhyb;

%beta = data.gene.beta;
%ip   = data.gene.ip;

%wdia = data.gene.wdia;
%wth  = data.gene.wth;
%nbar = data.gene.nbar;

h    = figure;
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])

subplot(4,1,1)

plot(tpnbi,pnbi/1e6,tplh,plh/1e6)
ylabel('MW')

legend('P_N_B_I','P_L_H')
axis([1 12 0 20])
subplot(4,1,2)

plot(tip,-ip/1e6,tnem,nem/1e19)

legend('I_p (MA)','<n_e> (10^1^9 m^-^2')
axis([1 12 0 2])
subplot(4,1,3)

plot(twdia,wdia/1e6,twe,we/1e6)

legend('W_d_i_a','W_e_l_e_c')
ylabel('MJ');
axis([1 12 0 inf])
subplot(4,1,4)

plot(tbetap,betap)

xlabel('times (s)')

axis([1 12 0 inf])
