%%%%%%%%  ARTIFICIAL NEURAL NETWORK BASED ON QUALIKIZ REGRESSION
%%%%%%%   ITG WITH KINETIC ELECTRONS. 4D INPUTS
%%%%%%%   APPROX 6 ORDERS OF MAGNITUDE FASTER THAN QUALIKIZ
%%%%%%%   ALSO CALCULATES GRADIENTS WITH RESPECT TO INPUTS
%%%%%%%%  JONATHAN CITRIN, JUAN REDONDO, SARAH BRETON, CLARISSE BOURDELLE, FREDERIC IMBEAUX

%%%%%%%%% WARNING, NETWORK VALIDITY RANGE: Ati [2 12]; ti/te [0.3 3]; q [1 5]; s [0.1 3]
%%%%%%%%% IN PRACTICE THE REGULARIZED NETWORKS CAN OFTEN REASONABLY EXTRAPOLATE, BUT CAVEAT EMPTOR!

%%%%%%%%% TRANSPORT MODEL DEVELOPED FOR USE IN RAPTOR (F. FELICI et al) and METIS (J.F. Artaud et al)

function varargout = qlkANN4Dkin(te,ti,tepa,tipa,q,shear,R0,a0,B0,Amain,pp,transflag)

%INPUTS
%te - electron temperature in eV 
%ti - ion temperature in eV
%tepa - electron temperature gradient: dte/da
%tipa - ion temperature gradient: dti/da
%q - safety factor
%shear - magnetic shear
%R0 - device major radius [m]
%a0 - device minor radius (Rmin-Rmax)/2 [m]
%B0 - vacuum B0 [T]
%Amain - Atomic number of main ion species
%pp - parameter structure for the module, see chi_qlkANN
%transport coefficent flag. If = 1, heat transport. If = 2, particle transport  
%NOTE: te, ti, tegrad, tigradcan all be profiles (same length each)

%OUTPUTS
%chi(e,i) - transport coefficient output. SI units: m^2/s
%dchi(e,i)_dte - chi(e,i) derivative with respect to te (te in eV)
%dchi(e,i)_dti - chi(e,i) derivative with respect to ti (ti in eV)
%dchi(e,i)_dtepa - chi(e,i) derivative with respect to dte (te in eV)
%dchi(e,i)_dtipa - chi(e,i) derivative with respect to dti (ti in eV)
%dchi(e,i)_dq - chi(e,i) derivative with respect to q
%dchi(e,i)_dshear - chi(e,i) derivative with respect to magnetic shear

%D - electron diffusivity. SI units: m^2/s
%V - electron pinch. SI units: m/s
%dD_dte - D derivative with respect to te (te in eV)
%dD_dti - D derivative with respect to ti (ti in eV)
%dD_dtepa - D derivative with respect to dte (te in eV)
%dD_dtipa - D derivative with respect to dti (ti in eV)
%dD_dq - D derivative with respect to q
%dD_dshear - D derivative with respect to magnetic shear
%dV_dte - V derivative with respect to te (te in eV)
%dV_dti - V derivative with respect to ti (ti in eV)
%dV_dtepa - V derivative with respect to dte (te in eV)
%dV_dtipa - V derivative with respect to dti (ti in eV)
%dV_dq - V derivative with respect to q
%dV_dshear - V derivative with respect to magnetic shear

  
%%%%%%BEGINNING OF CODE BODY

%NN trained for dimensionless R/LTi = -R*grad(Ti)/Ti , Ti/Te , q, magnetic shear.
%Not some as function input, thus, must initialize input vector

qel = 1.6e-19;
%Rescale temperatures to SI units
teSI  = te*qel; 
tiSI  = ti*qel;
teSIp = tepa*qel; 
tiSIp = tipa*qel;

%Define gyroBohm normalization. NN outputs are in GB normalized units
%chiGB = sqrt(mi)/(a*e^2*B^2)*Te^1.5
chifac = sqrt(Amain*1.67e-27)./(qel^2*B0.^2.*a0);

rlti = -R0 .* (tiSIp) ./ tiSI; 
tite = tiSI./teSI ;

%% saturate values if exceed training range
marg = 0.95; % margin: set to
if pp.constrains
    shearmin = pp.shearmin/marg; shearmax = pp.shearmax*marg;
    ishear = shear<shearmin | shear>shearmax;
    shear = min(max(shear,shearmin),shearmax);
else
    ishear = false(size(shear));
end

if pp.constrainq
    qmin=pp.qmin/marg; qmax = pp.qmax*marg;
    iq = q>qmax | q<qmin;
    q = min(max(q,qmin),qmax);
else
    iq = false(size(q));
end

if pp.constrainrlti
    rltimin=pp.rltimin/marg; rltimax = pp.rltimax*marg;
    irlti = rlti>rltimax | rlti<rltimin;
    rlti = min(max(rlti,rltimin),rltimax);
else
    irlti = false(size(rlti));
end

if pp.constraintite
    titemin=pp.titemin/marg; titemax = pp.titemax*marg;
    itite = tite>titemax | tite<titemin;
    tite = min(max(tite,titemin),titemax);
else
    itite = false(size(tite));
end

invec=[q,shear,tite,rlti]';
len=length(rlti);
%initialize input matrix in order expected by the trained network

%Normalize input matrix to match NN normalization
%NOTE: This is hardwired to match specific NN trained for this problem
nrhopoints = numel(q);

scalefac =[2 ; 1.45 ; 1.35 ; 5] ; 
dispfac  =[3 ; 1.55 ; 1.65 ; 7] ;

%Te powers needed for chiGB normalization
teSI12 = sqrt(teSI);
teSI32 = teSI12.*teSI;

%normalize input vector to match input normalization in NN
unorm = bsxfun(@times,bsxfun(@minus,invec,dispfac),1./scalefac);

%%shift R/LTi thgresholds for threshold matching
pp.b1_eef = pp.b1_eef + pp.shift_chi_e*pp.IW_eef(:,4)./scalefac(4);
pp.b1_ief = pp.b1_ief + pp.shift_chi_i*pp.IW_ief(:,4)./scalefac(4);

pp.b1_dfe = pp.b1_dfe + pp.shift_DV*pp.IW_dfe(:,4)./scalefac(4);
pp.b1_vce = pp.b1_vce + pp.shift_DV*pp.IW_vce(:,4)./scalefac(4);
pp.b1_vte = pp.b1_vte + pp.shift_DV*pp.IW_vte(:,4)./scalefac(4);



%%%Now the magic begins. Calculate separately the transport coefficients from the Neural Networks

if transflag == 1 %heat transport

%calculate chie regression neural network piece by piece, for all radii
  g_e=bsxfun(@plus,pp.IW_eef *unorm,pp.b1_eef);
  sg_e = sigmo(g_e);
  f_e=bsxfun(@plus,pp.L1W_eef *sg_e,pp.b2_eef);
  sf_e = sigmo(f_e);

  chieGB = (pp.L2W_eef*sf_e + pp.b3_eef)';

	       % if the NN output normalization was to between 0 and 1
  if pp.prepros_eef==0
    chieGB= chieGB.*(pp.max_qlk_eef-pp.min_qlk_eef)+ones(len,1)*pp.min_qlk_eef;
  end 

	       % if the NN output normalization was to between -1 a +1
  if pp.prepros_eef==-1
    chieGB= (chieGB+ones(size(chieGB))).*(ones(len,1)*pp.max_qlk_eef-ones(len,1)*pp.min_qlk_eef)/2 + ones(len,1)*pp.min_qlk_eef';
  end

				%set output vectors
  chie = chieGB.*teSI32.*chifac;

%%%%%%%%
%calculate chii regression neural network piece by piece, for all radii
  g_i=bsxfun(@plus,pp.IW_ief *unorm,pp.b1_ief);
  sg_i = sigmo(g_i);
  f_i=bsxfun(@plus,pp.L1W_ief *sg_i,pp.b2_ief);
  sf_i = sigmo(f_i);

  chiiGB = (pp.L2W_ief*sf_i + pp.b3_ief)';

	       % if the NN output normalization was to between 0 and 1
  if pp.prepros_ief==0
    chiiGB= chiiGB.*(pp.max_qlk_ief-pp.min_qlk_ief)+ones(len,1)*pp.min_qlk_ief;
  end 

	       % if the NN output normalization was to between -1 a +1
  if pp.prepros_ief==-1
    chiiGB= (chiiGB+ones(size(chiiGB))).*(ones(len,1)*pp.max_qlk_ief-ones(len,1)*pp.min_qlk_ief)/2 + ones(len,1)*pp.min_qlk_ief';
  end

				%set output vectors
  chii = chiiGB.*teSI32.*chifac;

  varargout{1} = chie; %chie output
  varargout{2} = chii; %chii output


  if nargout>3 
    
			      % derivatives w.r.t. NN input parameters
			      % evalation by chain rule
    
						   %First for chie
    dsgdg = (sg_e+1).^2.*exp(-2*g_e); % d(sig(g))/dg
    dsfdf = (sf_e+1).^2.*exp(-2*f_e); % d(sig(f))/df
    dgdu = bsxfun(@times,pp.IW_eef ,1./scalefac'); % dg/duphys
    
    dchieGB = zeros(nrhopoints,4);
    for iu=1:4 % loop over NN inputs
               % chain rule
      if pp.prepros_eef==0
        dchieGB(:,iu) = (pp.max_qlk_eef-pp.min_qlk_eef).*pp.L2W_eef*...
			bsxfun(@times,dsfdf,pp.L1W_eef *bsxfun(@times,dsgdg,dgdu(:,iu)));
      end
      if pp.prepros_eef==-1
        dchieGB(:,iu) = (pp.max_qlk_eef-pp.min_qlk_eef)./2.*pp.L2W_eef*...
			bsxfun(@times,dsfdf,pp.L1W_eef*bsxfun(@times,dsgdg,dgdu(:,iu)));
      end

    end
    
    dchieGB(itite,3) = 0; % mask for saturation
			  % derivatives w.r.t. physical inputs
    dchie_dte = (1.5*chieGB.*teSI12 - ...
                 tiSI./teSI12.*dchieGB(:,3))*chifac*qel;
    
    dchie_drlti = dchieGB(:,4);
    dchie_drlti(irlti) = 0; % mask for saturation
    
    dchie_dti = (teSI12.*dchieGB(:,3) - ...
                 teSI32.*(rlti)./tiSI.*dchie_drlti)*chifac*qel;

    dchie_dtepa = zeros(nrhopoints,1)*qel ;
    
    dchie_dtipa = -chifac.*R0./tiSI.*teSI32.*dchie_drlti*qel ;

    dchie_dq   = chifac.*teSI32.*dchieGB(:,1);
    dchie_dshear = chifac.*teSI32.*dchieGB(:,2);
    
    %% fixes for saturation
    dchie_dshear(ishear) = 0;
    dchie_dq(iq) = 0;
    
						   %Now for chii
    dsgdg = (sg_i+1).^2.*exp(-2*g_i); % d(sig(g))/dg
    dsfdf = (sf_i+1).^2.*exp(-2*f_i); % d(sig(f))/df
    dgdu = bsxfun(@times,pp.IW_ief ,1./scalefac'); % dg/duphys
    
    dchiiGB = zeros(nrhopoints,4);
    for iu=1:4 % loop over NN inputs
               % chain rule
      if pp.prepros_ief==0
        dchiiGB(:,iu) = (pp.max_qlk_ief-pp.min_qlk_ief).*pp.L2W_ief*...
			bsxfun(@times,dsfdf,pp.L1W_ief *bsxfun(@times,dsgdg,dgdu(:,iu)));
      end
      if pp.prepros_ief==-1
        dchiiGB(:,iu) = (pp.max_qlk_ief-pp.min_qlk_ief)./2.*pp.L2W_ief*...
			bsxfun(@times,dsfdf,pp.L1W_ief*bsxfun(@times,dsgdg,dgdu(:,iu)));
      end

    end
    
    dchiiGB(itite,3) = 0; % mask for saturation
    
    dchii_drlti = dchiiGB(:,4);
    dchii_drlti(irlti) = 0; % mask for saturation
    
				% derivatives w.r.t. physical inputs
    dchii_dte = (1.5*chiiGB.*teSI12 - ...
                 tiSI./teSI12.*dchiiGB(:,3))*chifac*qel;
    
    dchii_dti = (teSI12.*dchiiGB(:,3) - ...
                 teSI32.*(rlti)./tiSI.*dchii_drlti)*chifac*qel;

    dchii_dtepa = zeros(nrhopoints,1)*qel ;
    
    dchii_dtipa = -chifac.*R0./tiSI.*teSI32.*dchii_drlti*qel ;

    dchii_dq   = chifac.*teSI32.*dchiiGB(:,1);
    dchii_dshear = chifac.*teSI32.*dchiiGB(:,2);
    
    %% fixes for saturation
    dchii_dshear(ishear) = 0;
    dchii_dq(iq) = 0;
       
				%set output chie derivatives
    varargout{3} = dchie_dte;
    varargout{4} = dchie_dti;    
    varargout{5} = dchie_dtepa;
    varargout{6} = dchie_dtipa;
    varargout{7} = dchie_dq;
    varargout{8} = dchie_dshear;
				%set output chii derivatives
    varargout{9}  = dchii_dte;
    varargout{10} = dchii_dti;
    varargout{11} = dchii_dtepa;
    varargout{12} = dchii_dtipa;
    varargout{13} = dchii_dq;
    varargout{14} = dchii_dshear;
  end

elseif transflag == 2 %particle transport

% Particle transport coefficients were (accidentally) trained in SI units
% We must revert back to GB units using the dimensional units employed in the training, and then revert to SI
% using the dimensional units in our current simulation
% NN training set values: Te = 8 keV , a = 1 m, B0 = 3 T, Amain = 2, n = 5e19 [m^-3], R/Ln=2
  chiGBNN = (8e3.*qel).^1.5.*sqrt(2.*1.67e-27)./(qel.^2.*3.^2.*1);
  
%calculate De regression neural network piece by piece, for all radii
  g_de=bsxfun(@plus,pp.IW_dfe *unorm,pp.b1_dfe);
  sg_de = sigmo(g_de);
  f_de=bsxfun(@plus,pp.L1W_dfe *sg_de,pp.b2_dfe);
  sf_de = sigmo(f_de);

  deSI = (pp.L2W_dfe*sf_de + pp.b3_dfe)';

	       % if the NN output normalization was to between 0 and 1
  if pp.prepros_dfe==0
    deSI= deSI.*(pp.max_qlk_dfe-pp.min_qlk_dfe)+ones(len,1)*pp.min_qlk_dfe;
  end 

	       % if the NN output normalization was to between -1 a +1
  if pp.prepros_dfe==-1
    deSI= (deSI+ones(size(deSI))).*(ones(len,1)*pp.max_qlk_dfe-ones(len,1)*pp.min_qlk_dfe)/2 + ones(len,1)*pp.min_qlk_dfe';
  end

  %set output vectors
  deGB = deSI./chiGBNN; %convert from SI to GB with dimensional units used in training
  
  de=deGB.* teSI32.*chifac ; % convert to SI with actual simulation units

%%%%%%%%
%calculate Ve regression neural network piece by piece, for all radii
  g_vte=bsxfun(@plus,pp.IW_vte *unorm,pp.b1_vte);
  sg_vte = sigmo(g_vte);
  f_vte=bsxfun(@plus,pp.L1W_vte *sg_vte,pp.b2_vte);
  sf_vte = sigmo(f_vte);

  vteSI = (pp.L2W_vte*sf_vte + pp.b3_vte)';

  % if the NN output normalization was to between 0 and 1
  if pp.prepros_vte==0
    vteSI= vteSI.*(pp.max_qlk_vte-pp.min_qlk_vte)+ones(len,1)*pp.min_qlk_vte;
  end 

  % if the NN output normalization was to between -1 a +1
  if pp.prepros_vte==-1
    vteSI= (vteSI+ones(size(vteSI))).*(ones(len,1)*pp.max_qlk_vte-ones(len,1)*pp.min_qlk_vte)/2 + ones(len,1)*pp.min_qlk_vte';
  end

  %set output vectors

  vteGB = vteSI./chiGBNN.*3;  %VGB=VSI*R/chiGB;
  vte=vteGB.*teSI32.*chifac./R0;

  g_vce=bsxfun(@plus,pp.IW_vce *unorm,pp.b1_vce);
  sg_vce = sigmo(g_vce);
  f_vce=bsxfun(@plus,pp.L1W_vce *sg_vce,pp.b2_vce);
  sf_vce = sigmo(f_vce);

  vceSI = (pp.L2W_vce*sf_vce + pp.b3_vce)';

  % if the NN output normalization was to between 0 and 1
  if pp.prepros_vce==0
    vceSI= vceSI.*(pp.max_qlk_vce-pp.min_qlk_vce)+ones(len,1)*pp.min_qlk_vce;
  end 

  % if the NN output normalization was to between -1 a +1
  if pp.prepros_vce==-1
    vceSI= (vceSI+ones(size(vceSI))).*(ones(len,1)*pp.max_qlk_vce-ones(len,1)*pp.min_qlk_vce)/2 + ones(len,1)*pp.min_qlk_vce';
  end

  %set output vectors
  vceGB = vceSI./chiGBNN.*3;  %VGB=VSI*R/chiGB;
  vce=vceGB.*teSI32.*chifac./R0;


  if ~pp.allowpospinch
      varargout{1} = de;   %de output
      varargout{2} = min(0,vce+vte); %Ve output  
  else
      varargout{1} = de;   %de output
      varargout{2} = vce+vte; %Ve output
  end
  
  if nargout>3 
    
    % derivatives w.r.t. NN input parameters
    % evalation by chain rule
    
    %First for chie
    dsgdg = (sg_de+1).^2.*exp(-2*g_de); % d(sig(g))/dg
    dsfdf = (sf_de+1).^2.*exp(-2*f_de); % d(sig(f))/df
    dgdu = bsxfun(@times,pp.IW_dfe ,1./scalefac'); % dg/duphys
    
    ddeGB = zeros(nrhopoints,4);
    for iu=1:4 % loop over NN inputs
               % chain rule
      if pp.prepros_dfe==0
        ddeGB(:,iu) = ((pp.max_qlk_dfe-pp.min_qlk_dfe).*pp.L2W_dfe*...
			bsxfun(@times,dsfdf,pp.L1W_dfe *bsxfun(@times,dsgdg,dgdu(:,iu))) )./chiGBNN;
      end
      if pp.prepros_dfe==-1
        ddeGB(:,iu) = ((pp.max_qlk_dfe-pp.min_qlk_dfe)./2.*pp.L2W_dfe*...
			bsxfun(@times,dsfdf,pp.L1W_dfe*bsxfun(@times,dsgdg,dgdu(:,iu))))./chiGBNN;
      end

    end
    
    ddeGB(itite,3) = 0; % mask for saturation
			  % derivatives w.r.t. physical inputs
    dde_dte = (1.5*deGB.*teSI12 - ...
                 tiSI./teSI12.*ddeGB(:,3))*chifac*qel;
    
    dde_drlti = ddeGB(:,4);
    dde_drlti(irlti) = 0; % mask for saturation
    
    dde_dti = (teSI12.*ddeGB(:,3) - ...
                 teSI32.*(rlti)./tiSI.*dde_drlti)*chifac*qel;

    dde_dtepa = zeros(nrhopoints,1)*qel ;
    
    dde_dtipa = -chifac.*R0./tiSI.*teSI32.*dde_drlti*qel ;

    dde_dq   = chifac.*teSI32.*ddeGB(:,1);
    dde_dshear = chifac.*teSI32.*ddeGB(:,2);
    
    %% fixes for saturation
    dde_dshear(ishear) = 0;
    dde_dq(iq) = 0;
    
    %Now for vte
    dsgdg = (sg_vte+1).^2.*exp(-2*g_vte); % d(sig(g))/dg
    dsfdf = (sf_vte+1).^2.*exp(-2*f_vte); % d(sig(f))/df
    dgdu = bsxfun(@times,pp.IW_vte ,1./scalefac'); % dg/duphys
    
    dvteGB = zeros(nrhopoints,4);
    for iu=1:4 % loop over NN inputs
               % chain rule
      if pp.prepros_vte==0
        dvteGB(:,iu) = ((pp.max_qlk_vte-pp.min_qlk_vte).*pp.L2W_vte*...
			bsxfun(@times,dsfdf,pp.L1W_vte *bsxfun(@times,dsgdg,dgdu(:,iu))) )./chiGBNN.*3;
      end
      if pp.prepros_vte==-1
        dvteGB(:,iu) = ((pp.max_qlk_vte-pp.min_qlk_vte)./2.*pp.L2W_vte*...
			bsxfun(@times,dsfdf,pp.L1W_vte*bsxfun(@times,dsgdg,dgdu(:,iu))))./chiGBNN.*3;
      end

    end
    
    dvteGB(itite,3) = 0; % mask for saturation
			  % derivatives w.r.t. physical inputs
    dvte_dte = (1.5*vteGB.*teSI12 - ...
                 tiSI./teSI12.*dvteGB(:,3))*chifac./R0*qel;
    
    dvte_drlti = dvteGB(:,4);
    dvte_drlti(irlti) = 0; % mask for saturation
    
    dvte_dti = (teSI12.*dvteGB(:,3) - ...
                 teSI32.*(rlti)./tiSI.*dvte_drlti)*chifac./R0*qel;

    dvte_dtepa = zeros(nrhopoints,1)*qel ;
    
    dvte_dtipa = -chifac./R0.*R0./tiSI.*teSI32.*dvte_drlti*qel ;

    dvte_dq   = chifac./R0.*teSI32.*dvteGB(:,1);
    dvte_dshear = chifac./R0.*teSI32.*dvteGB(:,2);
    
    %% fixes for saturation
    dvte_dshear(ishear) = 0;
    dvte_dq(iq) = 0;

    %Now for vce
    dsgdg = (sg_vce+1).^2.*exp(-2*g_vce); % d(sig(g))/dg
    dsfdf = (sf_vce+1).^2.*exp(-2*f_vce); % d(sig(f))/df
    dgdu = bsxfun(@times,pp.IW_vce ,1./scalefac'); % dg/duphys
    
    dvceGB = zeros(nrhopoints,4);
    for iu=1:4 % loop over NN inputs
               % chain rule
      if pp.prepros_vce==0
        dvceGB(:,iu) = ((pp.max_qlk_vce-pp.min_qlk_vce).*pp.L2W_vce*...
			bsxfun(@times,dsfdf,pp.L1W_vce *bsxfun(@times,dsgdg,dgdu(:,iu))) )./chiGBNN.*3;
      end
      if pp.prepros_vce==-1
        dvceGB(:,iu) = ((pp.max_qlk_vce-pp.min_qlk_vce)./2.*pp.L2W_vce*...
			bsxfun(@times,dsfdf,pp.L1W_vce*bsxfun(@times,dsgdg,dgdu(:,iu))))./chiGBNN.*3;
      end

    end
    
    dvceGB(itite,3) = 0; % mask for saturation
			  % vcerivatives w.r.t. physical inputs
    dvce_dte = (1.5*vceGB.*teSI12 - ...
                 tiSI./teSI12.*dvceGB(:,3))*chifac./R0*qel;
    
    dvce_drlti = dvceGB(:,4);
    dvce_drlti(irlti) = 0; % mask for saturation
    
    dvce_dti = (teSI12.*dvceGB(:,3) - ...
                 teSI32.*(rlti)./tiSI.*dvce_drlti)*chifac./R0*qel;

    dvce_dtepa = zeros(nrhopoints,1)*qel ;
    
    dvce_dtipa = -chifac./R0.*R0./tiSI.*teSI32.*dvce_drlti*qel ;

    dvce_dq   = chifac./R0.*teSI32.*dvceGB(:,1);
    dvce_dshear = chifac./R0.*teSI32.*dvceGB(:,2);
    
    %% fixes for saturation
    dvce_dshear(ishear) = 0;
    dvce_dq(iq) = 0;

    %set output de derivatives
    varargout{3} = dde_dte;
    varargout{4} = dde_dti;    
    varargout{5} = dde_dtepa;
    varargout{6} = dde_dtipa;
    varargout{7} = dde_dq;
    varargout{8} = dde_dshear;

    %set output V derivatives
    varargout{9}  = dvte_dte+dvce_dte;
    varargout{10} = dvte_dti+dvce_dti;
    varargout{11} = dvte_dtepa+dvce_dtepa;
    varargout{12} = dvte_dtipa+dvce_dtipa;
    varargout{13} = dvte_dq+dvce_dq;
    varargout{14} = dvte_dshear+dvce_dshear;
  end
end

return

%nonlinear neuron function
function s = sigmo(x)
s= 2./(1+exp(-2.*x))-1; %define sigmoid function (nonlinear transfer functions in NN)
return

