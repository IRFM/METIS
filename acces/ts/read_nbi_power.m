function [tcronos_out,pnbi_cronos,beam_energy,frac,amass] = read_nbi_power(shot,tcronos,idisplay)

%% -------------------------------------------------------------------------------------
%% CALCULATE THE NBI INJECTED POWER FOR TORE SUPRA, READING THE TS DATABASE
%% -------------------------------------------------------------------------------------
%% HOW TO CALL IT FROM A CRONOS WORKSPACE:
%% [tcronos_out,data.cons.idn,ene,frac,amass]=read_nbi_power(param.from.shot.num,data.gene.temps);
%% OPTIONAL ARGUMENTS: 
%%    - TCRONOS  = VECTOR FOR TIMEBASE (DEFAULT = TCRONOS = TS DATABASE TIME)
%%    - IDISPLAY = 0 OR 1 IN ORDER TO DISPLAY THE NBI POWER VERSUS TIME (DEFAULT = 0)
%% -------------------------------------------------------------------------------------
%% OUTPUT:
%%    - TCRONOS_OUT:     OUTPUT TIME VECTOR (S)
%%    - PNBI_CRONOS: NBI POWER VERSUS TIME (W)
%%    - BEAM_ENERGY: NEUTRAL BEAM ENERGY (KEV)
%%    - FRAC(3):     PARTICLE FRACTIONS (-)
%%    - AMASS:       MASS NUMBER OF THE INTECTED ATOM (-)
%% -------------------------------------------------------------------------------------
%% SINCE RELEVENT SIGNALS ARE NOT CALIBRATED IN THE DATABASE, WE
%% USE OTHER SIGNALS, LESS APPROPRIATE BUT RE-CALIBRATED, THEREFORE
%% CLOSER TO THE REALITY
%% -------------------------------------------------------------------------------------

shot = floor(shot);

if nargin < 3
  idisplay = 0;
end

%% -----------------------------------------------
%% READ NEUTRAL BEAM CURRENT FROM TS DATABASE (A)
%% -----------------------------------------------
[sidn,tbase] = tsbase(shot,'sidn');

if nargin == 1
  tcronos = tbase;
end
tcronos_out = tcronos;

%% SECURITY TO AVOID CRASHED WHEN DATABASE IS NOT FILLED
if length(sidn)>1

  sidn = abs(sidn);
  
%% -------------------------------------------------------------
%% LIMITS OF THE COMMON TIME INTERVAL BETWEEN TCRONOS AND TBASE
%% -------------------------------------------------------------
tmin = max(tcronos(1),tbase(1));
tmax = min(tcronos(end),tbase(end));

%% --------------------------------------------
%% CORRESPONDING INDICES FOR TCRONOS AND TBASE
%% --------------------------------------------
idxtminc = closest(tcronos,tmin);
idxtmaxc = closest(tcronos,tmax);
idxtminb = closest(tbase,tmin);
idxtmaxb = closest(tbase,tmax);

%% -----------------------------------------------------------------
%% CHECK THE NEUTRON RATE TO EVALUATE IF THIS IS A HYDROGEN OR
%% DEUTERIUM BEAM (NEUTRON RATE INCREASES DURING D-INJECTION DUE TO
%% D-D REACTIONS, WHILE IT IS CONSTANT DURING H-INJECTION)
%% -----------------------------------------------------------------
[fn,tfn,rfn,cfn]= tsbase(shot,'gfluntn');

if length(tfn)>0
  tfn = tfn(:,2);
  flux_neutron = 10.^fn(:,2);

  %% REFERENCE NEUTRON RATE = THE ONE AT T=3 SECOND
  neutron_ref = flux_neutron(closest(tfn,3));

  %% NEUTRON RATE DURING NBI:
  t_nbi_max = tbase(find(sidn==max(sidn)));
  t_nbi_max = t_nbi_max(1);
  idx_nbi   = closest(tfn,t_nbi_max);
  neutron_ratio = flux_neutron(idx_nbi)/neutron_ref;
  
else
  
  neutron_ratio = 1;
  disp('----------------------------------------------------')
  disp('Neutron rate not available => HYDROGEN beam assumed')
  disp('----------------------------------------------------')
  
end

%% ------------------------------------------------------------
%% APPROXIMATE NEUTRALIZER VOLTAGE (= BEAM ENERGY) (KEV OR KV)
%% ------------------------------------------------------------
if shot < 45200
  ene_h = 55;
  ene_d = 65;
else
  ene_h = 60;
  ene_d = 80;
end
if neutron_ratio > 10
  disp('---------------------')
  disp('DEUTERIUM BEAM')
  disp(['Neut ratio = ',num2str(neutron_ratio)])
  disp('---------------------')
  beam_energy = ene_d; % DEUTERIUM
  amass = 2;
else
  if length(tfn)>0
    disp('---------------------')
    disp('HYDROGEN BEAM')
    disp(['Neut ratio = ',num2str(neutron_ratio)])
    disp('---------------------')
  end
  beam_energy = ene_h; % HYDROGEN
  amass = 1;
end

%% ------------------------------------------------------
%% NEUTRALIZATION COEFFICIENT FOR HYDROGEN AND DEUTERIUM
%% ------------------------------------------------------
[coef_h,coef_d,enevec] = coef_neutralisation;
idx_ene = closest(enevec,beam_energy);
coefh = coef_h(idx_ene);
coefd = coef_d(idx_ene);
if neutron_ratio > 10
  coef = coefd; % DEUTERIUM
else
  coef = coefh; % HYDROGEN
end

%% ---------------------
%% DEDUCE NBI POWER (W)
%% ---------------------
pnbi = sidn.*beam_energy*1e3*0.85*coef; % (85% = RE-IONIZATION, TRANSMISSION ETC)
pnbi(find(pnbi<0)) = 0;

%% --------------------------------------------------------------
%% INTERPOLATE NBI POWER FROM TS DATABASE TO CRONOS TIME VECTOR
%% --------------------------------------------------------------
pnbi_cronos = zeros(size(tcronos));
if idxtmaxb > idxtminb
  pnbi_cronos(idxtminc:idxtmaxc) = interp1(tbase(idxtminb:idxtmaxb),pnbi(idxtminb:idxtmaxb),tcronos(idxtminc:idxtmaxc),'nearest','extrap');
end
  
%% -----------------------------
%% PARTICLE AND POWER FRACTIONS
%% -----------------------------
fracn(1) = 0.45;
fracn(2) = 0.12;
fracn(3) = 0.43;
fracp(1) = 0.81;
fracp(2) = 0.06;
fracp(3) = 0.13;

frac = fracn; % WE WANT PARTICLE FRACTIONS IN CRONOS

%% --------------------------------------------------------------
%% DISPLAY THE NBI POWER ON THE TORE SUPRA AND CRONOS TIME BASE,
%% FOR VERIFICATIONS
%% --------------------------------------------------------------
if idisplay == 1

  linewidth = 2;
  fontsize  = 15;

  clf
  h=axes;
  set(h,'FontSize',fontsize)
  hold on ; grid on

  idx  = find(pnbi_cronos>0);
  idx1 = idx(1)-10;
  idx2 = idx(end)+10;
  if idx1 < 1
    idx1 = 1;
  end
  if idx2 > length(tcronos)
    idx2 = idx(end);
  end

  plot(tbase,pnbi*1e-3,'LineWidth',linewidth)
  plot(tcronos,pnbi_cronos*1e-3,'r','LineWidth',linewidth)

  warning('off','MATLAB:dispatcher:InexactMatch')
  
  xlabel('Time (s)','FontSize',fontsize)
  ylabel('NBI injected power (kW)','FontSize',fontsize)
  h=legend('TS timebase','CRONOS timebase','Location','NorthEast');
  set(h,'FontSize',fontsize)
  Title(['Shot #',num2str(shot)],'FontSize',fontsize)
  
  set(gca,'Xlim',[tcronos(idx1) tcronos(idx2)])
  hold off

end % TEST IDISPLAY

else

  disp('NO AVAILABLE NBI DATA FOR THIS SHOT')
  pnbi_cronos = zeros(size(tcronos_out));
  beam_energy = 0;
  frac(1:3)   = 0;
  amass       = 0;
  
end % TEST LENGTH(SIDN)>1


%eval(sprintf('print -djpeg100 nbi_power_%s.jpeg',num2str(shot)))

%[svhtneut,tbase] = tsbase(shot,'svhtneut');
%[svht,tbase]     = tsbase(shot,'svht');
%[sineut,tbase]   = tsbase(shot,'sineut');
%[sirecup,tbase]  = tsbase(shot,'sirecup');
%[sig2,tbase]     = tsbase(shot,'sig2');
%cur_ions = abs(sineut-sineut(1))+abs(sirecup-sirecup(1));
%beam_energy = abs(mean(svhtneut));

