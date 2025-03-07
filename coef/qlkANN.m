%%%%%%%%  ARTIFICIAL NEURAL NETWORK BASED ON QUALIKIZ REGRESSION
%%%%%%%   ITG WITH ADIABATIC ELECTRONS
%%%%%%%   APPROX 5 ORDERS OF MAGNITUDE FASTER THAN QUALIKIZ
%%%%%%%   ALSO CALCULATES GRADIENTS WITH RESPECT TO INPUTS
%%%%%%%%  SARAH BRETON,  CLARISSE BOURDELLE, JONATHAN CITRIN, FREDERIC IMBEAUX

%%%%%%%%% WARNING, NETWORK VALIDITY RANGE: Ati [2 12]; ti/te [0.3 3]; q [1 5]; s [0.1 3]
%%%%%%%%% IN PRACTICE THE REGULARIZED NETWORKS CAN OFTEN REASONABLY EXTRAPOLATE, BUT CAVEAT EMPTOR!

%%%%%%%%% TRANSPORT MODEL DEVELOPED FOR USE IN RAPTOR (F. FELICI et al) and CRONOS (J.F. Artaud et al)

%NOTES:
%1. Since typically qi/qe~2-3 , and kinetic electron qi is typically ~3 more than adiabatic electron qi,
%   we make the assumption that the flux output here can be considered as qe from kin ele ITG 
%2. This function is hardwired for 2 hidden layers and 3 biases. Can generalize if necessary in future

function varargout = qlkANN(te,ti,tep,tip,q,shear,R0,a0,B0,Amain,pp)

%INPUTS
%te - electron temperature in eV 
%ti - ion temperature in eV
%tegrad - electron temperature gradient. In general geometry, recommended to take LFS midplane grad
%tigrad - ion temperature gradient. In general geometry, recommended to take LFS midplane grad
%q - safety factor
%shear - magnetic shear
%pp - parameter structure for the module, see chi_qlkANN
%NOTE: te, ti, tegrad, tigradcan all be profiles (same length each)

%OUTPUTS
%chie - electron heat conductivity [m^2/s]. Note: electron heat flux = chie * ne *tegrad [W/m^2]
%dchie_dte - chie derivative with respect to te (te in eV)
%dchie_dti - chie derivative with respect to ti (ti in eV)
%dchie_dtegrad - chie derivative with respect to dte (te in eV)
%dchie_dtigrad - chie derivative with respect to dti (ti in eV)
%dchie_dq - chie derivative with respect to q
%dchie_dshear - chie derivative with respect to shear

%%%%%%BEGINNING OF CODE BODY

%NN trained for dimensionless R/LTi = -R*grad(Ti)/Ti , Ti/Te , q, shear.
%Not some as function input, thus, must initialize input vector

qel = 1.6e-19;
%Rescale temperatures to SI units
teSI  = te*qel; 
tiSI  = ti*qel;
teSIp = tep*qel; 
tiSIp = tip*qel;

%%

%Define gyroBohm normalization. NN outputs are in GB normalized units
%chiGB = sqrt(mi)/(a*e^2*B^2)*Te^1.5
chifac = sqrt(Amain*1.67e-27)/(qel^2*B0^2*a0);

%rlti = -R0 .* (tiSIp/a0) ./ tiSI ;
rlti = -R0 .* (tiSIp) ./ tiSI ;
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
%

invec=[rlti,q,tite,shear]' ;  %initialize input matrix in order expected by the trained network

%Normalize input matrix to match NN normalization
%NOTE: This is hardwired to match specific NN trained for this problem
nrhopoints = numel(q);

scalefac = [5;2;1.35;1.45];
dispfac  = [7;3;1.65;1.55];


%OUTPUT NORMALIZATION (hardwired to the training set used)
max_qlk=28.378817;

%normalize input vector to match input normalization in NN
unorm = bsxfun(@times,bsxfun(@minus,invec,dispfac),1./scalefac);

%calculate classification neural network piece by piece, for all radii
g=bsxfun(@plus,pp.IWcla*unorm,pp.b1cla);
sg = sigmo(g);
f=bsxfun(@plus,pp.L1Wcla*sg,pp.b2cla);
sf = sigmo(f);

class = sigmo((pp.L2Wcla*sf + pp.b3cla)');

class(class <= 0.5) = 0;
class(class > 0.5) = 1;

%calculate regression neural network piece by piece, for all radii
g=bsxfun(@plus,pp.IWreg*unorm,pp.b1reg);
sg = sigmo(g);
f=bsxfun(@plus,pp.L1Wreg*sg,pp.b2reg);
sf = sigmo(f);

chieGB = max_qlk*class.*(pp.L2Wreg*sf + pp.b3reg)';

teSI12 = sqrt(teSI);
teSI32 = teSI12.*teSI;

%set output vectors
chie = chieGB.*teSI32.*chifac;

varargout{1} = chie; %adiabatic electron assumption for chie
varargout{2} = 3*chie; % adiabatic electron assumption for chii

if nargout>2 %advanced use if gradients needed
    
    % drivatives w.r.t. NN input parameters
    % evalation by chain rule
    dsgdg = (sg+1).^2.*exp(-2*g); % d(sig(g))/dg
    dsfdf = (sf+1).^2.*exp(-2*f); % d(sig(f))/df
    dgdu = bsxfun(@times,pp.IWreg,1./scalefac'); % dg/duphys
    
    dchieGB = zeros(nrhopoints,4);
    for iu=1:4 % loop over NN inputs
        % chain rule
        dchieGB(:,iu) = max_qlk*pp.L2Wreg*...
        bsxfun(@times,dsfdf,pp.L1Wreg*bsxfun(@times,dsgdg,dgdu(:,iu)));
    end
    
    % derivatives w.r.t. physical inputs
    dchie_dte = class.*(1.5*chieGB.*teSI12 - ...
                  tiSI./teSI12.*dchieGB(:,3))*chifac*qel;
    dchie_dti = class.*(teSI12.*dchieGB(:,3) - ...
                  teSI32.*(rlti)./tiSI.*dchieGB(:,1))*chifac*qel;

    dchie_dtep = class.*zeros(nrhopoints,1)*qel ;
    dchie_dtip = class.*-chifac.*R0./tiSI.*teSI32.*dchieGB(:,1)*qel ;

    dchie_dq   = class.*chifac.*teSI32.*dchieGB(:,2);
    dchie_dshear = class.*chifac.*teSI32.*dchieGB(:,4);
    
    %% fixes for saturation
    dchie_dshear(ishear) = 0;
    dchie_dq(iq) = 0;

    varargout{3} = dchie_dte;
    varargout{4} = dchie_dti;
    varargout{5} = dchie_dtep;
    varargout{6} = dchie_dtip;
    varargout{7} = dchie_dq;
    varargout{8} = dchie_dshear;
    
end

return

%nonlinear neuron function
function s = sigmo(x)
s= 2./(1+exp(-2.*x))-1; %define sigmoid function (nonlinear transfer functions in NN)
return

