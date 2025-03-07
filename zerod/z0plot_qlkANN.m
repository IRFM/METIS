function z0plot_qlkANN(internal)
mode = [];
noms__ = fieldnames(internal);
for k__ =1:length(noms__)
  eval(sprintf('%s = internal.%s;',noms__{k__},noms__{k__}));
end

switch mode
case 'full'

      hz =findobj(0,'type','figure','tag','flux_net_kin');
      if isempty(hz)
		hz=figure('tag','flux_net_kin','name','Flux comparison (ITG part only)');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',3,'color',[1 1 1])

      subplot(2,2,1)
      zplotprof(gca,profil.temps,profil.xli,profil.qe,'color','b');
      zplotprof(gca,profil.temps,profil.xli,flux_qe,'color','r');
      zplotprof(gca,profil.temps,profil.xli,(flux_qe -  flux_qe_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(flux_qe +  flux_qe_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('flux_qe (W/m^{-3})');
      legend('METIS','QLKNN-4Dkin');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,2)
      zplotprof(gca,profil.temps,profil.xli,profil.qi,'color','b');
      zplotprof(gca,profil.temps,profil.xli,flux_qi,'color','r');
      zplotprof(gca,profil.temps,profil.xli,(flux_qi -  flux_qi_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(flux_qi -  flux_qi_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('flux_qi (W/m^{-3})');
      legend('METIS','QLKNN-4Dkin');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,3)
      zplotprof(gca,profil.temps,profil.xli,profil.ge,'color','b');
      zplotprof(gca,profil.temps,profil.xli,flux_ne,'color','r');
      zplotprof(gca,profil.temps,profil.xli,(flux_ne - flux_ne_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(flux_ne + flux_ne_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('flux_ne (e^-/m^{-3}/s)');
      legend('METIS','QLKNN-4Dkin');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      hz =findobj(0,'type','figure','tag','Chi_net_kin');
      if isempty(hz)
		hz=figure('tag','Chi_net_kin','name','Transport coefficients comparison (ITG part only)');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',3,'color',[1 1 1])

      subplot(2,2,1)
      zplotprof(gca,profil.temps,profil.xli,chie_ref,'color','b');
      zplotprof(gca,profil.temps,profil.xli,chie,'color','r');
      zplotprof(gca,profil.temps,profil.xli,chie_loop,'color','k');
      zplotprof(gca,profil.temps,profil.xli,(chie-chie_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(chie+chie_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('Chi_e (m^2/s)');
      legend('METIS','QLKNN-4Dkin with METIS profiles','QLKNN-4Dkin prediction');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,2)
      zplotprof(gca,profil.temps,profil.xli,chii_ref,'color','b');
      zplotprof(gca,profil.temps,profil.xli,chii,'color','r');
      zplotprof(gca,profil.temps,profil.xli,chii_loop,'color','k');
      zplotprof(gca,profil.temps,profil.xli,chii_neo,'color','g');
      zplotprof(gca,profil.temps,profil.xli,(chii-chii_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(chii+chii_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('Chi_i (m^2/s)');
      legend('METIS','QLKNN-4Dkin with METIS profiles','QLKNN-4Dkin prediction','Neo');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,3)
      zplotprof(gca,profil.temps,profil.xli,D_ref,'color','b');
      zplotprof(gca,profil.temps,profil.xli,D,'color','r');
      zplotprof(gca,profil.temps,profil.xli,D_loop,'color','k');
      zplotprof(gca,profil.temps,profil.xli,(D-D_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(D+D_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('D (m^2/s)');
      legend('METIS','QLKNN-4Dkin with METIS profiles','QLKNN-4Dkin prediction');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,4)
      VsD_err = abs(V_err ./ D) + abs(V.* D_err ./ D .^ 2);
      zplotprof(gca,profil.temps,profil.xli,V_ref ./ D_ref,'color','b');
      zplotprof(gca,profil.temps,profil.xli,V ./ D,'color','r');
      zplotprof(gca,profil.temps,profil.xli,V_loop ./ D_loop,'color','k');
      zplotprof(gca,profil.temps,profil.xli,(V ./ D - VsD_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(V ./ D + VsD_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('V / D (m)');
      legend('METIS','QLKNN-4Dkin with METIS profiles','QLKNN-4Dkin prediction');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon
end

hz =findobj(0,'type','figure','tag','Prof_net_kin');
if isempty(hz)
  	  hz=figure('tag','Prof_net_kin','name','Profiles prediction qlkAnn_kin_e (ITG part only, fixed sources but equipartition, no time dependance)');
else
  	  figure(hz);
end
clf
set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(2,2,1)
zplotprof(gca,profil.temps,profil.xli,profil.tep ./ 1e3,'color','b');
zplotprof(gca,profil.temps,profil.xli,te_loop ./ 1e3,'color','r');
zplotprof(gca,profil.temps,profil.xli,mask,'color','g');
set(gca,'YScale','linear');
ylabel('T_e (keV)');
%legend('METIS','QLKNN-4Dkin','interval of computation with QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,2)
zplotprof(gca,profil.temps,profil.xli,profil.tip ./ 1e3,'color','b');
zplotprof(gca,profil.temps,profil.xli,ti_loop ./ 1e3,'color','r');
zplotprof(gca,profil.temps,profil.xli,mask,'color','g');
set(gca,'YScale','linear');
ylabel('T_i (keV)');
%legend('METIS','QLKNN-4Dkin','interval of computation with QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,3)
zplotprof(gca,profil.temps,profil.xli,profil.nep ./ 1e19,'color','b');
zplotprof(gca,profil.temps,profil.xli,ne_loop ./ 1e19,'color','r');
zplotprof(gca,profil.temps,profil.xli,mask,'color','g');
set(gca,'YScale','linear');
ylabel('n_e (10^{19} m^{-3})');
%legend('METIS','QLKNN-4Dkin','interval of computation with QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,4)
zplotprof(gca,profil.temps,profil.xli,profil.qei ./ 1e6,'color','b');
zplotprof(gca,profil.temps,profil.xli, Qei ./ 1e6,'color','r');
zplotprof(gca,profil.temps,profil.xli,mask,'color','g');
set(gca,'YScale','linear');
ylabel('Q_{e,i} (MW)');
legend('METIS','QLKNN-4Dkin','interval of computation with QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

switch mode
case 'full'
    hz =findobj(0,'type','figure','tag','Wth_net_kin');
    if isempty(hz)
	      hz=figure('tag','Wth_net_kin','name','Energy content prediction qlkAnn_kin_e (ITG part only, fixed sources but equipartition, no time dependance)');
    else
	      figure(hz);
    end
    clf
    set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',3,'color',[1 1 1])

    subplot(2,1,1)
    plot(profil.temps,wth ./ 1e6,'b',profil.temps,wth_nn ./ 1e6,'r');
    xlabel('time (s)');
    ylabel('W_{th} (MJ)');
    legend('METIS','QLKNN-4Dkin');
    zoom on
    subplot(2,1,2);
    semilogy(profil.temps,wth_nn ./ max(1,wth),'b',profil.temps,factor_time,'r',profil.temps,ones(size(profil.temps)),'g');
    xlabel('time (s)');
    legend('W_{th,QLKNN-4Dkin} / W_{th,METIS}','time factor');
    z0loglin(gca);
    zoom on
    joint_axes(hz,2);
end
