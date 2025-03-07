function [rgaz,rdcx,rvis] = lit_nhnd_ts(shot)

 
 
% lecture injection de gaz 
try
	[debit,tg]=tsbase(shot,'gmesdeb'); 
	tg=tg(:,1);
	debD2=debit(:,6)+debit(:,7);debD2=smooth(debD2,5)/1000;
	debH2=debit(:,1)+debit(:,2);debH2=smooth(debH2,5)/1000;
	% nl
	[gnl,tnl] = tsbase(shot,'gnl');
	nl = max(gnl,[],2);
	nl = interp1(tnl,nl,tg,'nearest');
	nl(~isfinite(nl)) = 0;
	% injection de gaz
	intD2 = trapz(tg,2 .* debD2 .* (debD2 >0) .* 3e20 + nl ./ 7e-3)+eps;
	intH2 = trapz(tg,2.* debH2.* (debH2 >0) .* 3e20);
	rgaz = intH2 ./ (intH2 + intD2);
catch
	rgaz =[];
end

try
	% rapport isotopique de l'ANR 2
	[dcx,tdcx]=tsbase(shot,'SRI2');
	[pfci,tpfci]=tsbase(shot,'spuiss');
	pfci        = interp1(tpfci,pfci,tdcx,'nearest');
	pfci(~isfinite(pfci)) = 0;
	dcx(pfci > 0.5) = 0; 
	indnok = find(dcx < 0.005);
	dcx(indnok) =[];
	tdcx(indnok) =[];
	rdcx = mean(dcx);
catch
	rdcx = [];
end

% rapport isotopique de DVIS
if nargout >2
	try
		[H_HD,H_HDcorr,tdata]=TS_H_DH(shot,5);
		x = 0:0.005:0.3;
		n=hist(H_HDcorr,x);
		rvis = x(min(find(n == max(n))));
	catch
		rvis = [];
	end
end

 

function [H_HD,H_HDcorr,tdata]=TS_H_DH(shot,iview)

global H_HD H_HDcorr tdata  dh_ratio_disp_var


%shot=32543;
[cad,cad1,nbtr1,nbexpo,nbexpo1,lgonde,reseau,fente] = ...
      tsmat(shot,'DVISR;GENE_FONCT;tps_expo1',...
                 'DVISR;GENE_FONCT;tps_expo2',...
                 'DVISR;GENE_FONCT;nb_train',...
                 'DVISR;GENE_FONCT;nb_expo1',...
                 'DVISR;GENE_FONCT;nb_expo2',...
		 'DVISR;MM4006;lg_onde',...
       	         'DVISR;MM4006;reseau',...		     
		 'DVISR;MM4006;fente');
		 
if lgonde~=6570
    disp('   Wrong wavelength set on spectrometer')    
    return
end;
		 
%[voie1,voie2,voie3,voie4,voie5,voie6,voie7,voie8] = ...
%      tsmat(shot,'DVISR;POS_FIBRE;voie_1',...
%                 'DVISR;POS_FIBRE;voie_2',...
%		 'DVISR;POS_FIBRE;voie_3',...
%		 'DVISR;POS_FIBRE;voie_4',...
%		 'DVISR;POS_FIBRE;voie_5',...
%		 'DVISR;POS_FIBRE;voie_6',...
%		 'DVISR;POS_FIBRE;voie_7',...		     
%		 'DVISR;POS_VIBRE;voie_8');
		 
         npix = 512;
         cad  = cad/1e6;
         cad1 = cad1/1e6;
%        Vecteur de temps d'expo;
         texp = ones(1,nbexpo)*cad;
         texp1 = ones(1,nbexpo1)*cad1;
         if nbtr1 == 2
            texp = [texp texp1];
         elseif nbtr1 == 3
            texp = [texp texp1 texp];
         elseif nbtr1 == 4
            texp = [texp texp1 texp texp1];
         elseif nbtr1 == 5
            texp = [texp texp1 texp texp1 texp];
         elseif nbtr1 >= 6
            texp = [texp texp1 texp texp1 texp texp1];
         elseif nbtr1 == 7
            texp = [texp texp1 texp texp1 texp texp1 texp];
         elseif nbtr1 == 8
            texp = [texp texp1 texp texp1 texp texp1 texp texp1];
         elseif nbtr1 == 9
            texp = [texp texp1 texp texp1 texp texp1 texp texp1 texp];
         elseif nbtr1 == 10
            texp = [texp texp1 texp texp1 texp texp1 texp texp1 texp texp1];
         elseif nbtr1 >= 11
            disp('get_SpMeterSetup - nbtr1>=11 - Le vecteur de temps d''expo ne peut pas �re construit - Modifiez la routine')
         end
         texp = ones(npix,1)*texp;
	 
	 rdata=tsbase(shot,'scampir1');
	 for ivoie=2:8;
	 rdata=[rdata,tsbase(shot,['scampir',int2str(ivoie)])];
	 end
	 spectra=reshape(rdata,npix,(length(rdata)./npix),8);
	 
if reseau==600, dispersion = -2.3784e-6*lgonde + 0.2200;end % Mesures 09/2001
if reseau == 1200, dispersion = 1.3213e-14*lgonde.^3 - 5.6398e-10*lgonde.^2 ...
 + 0.0200e-4*lgonde + 0.1008;end    % Mesures 09/2001

wlength=lgonde-255*dispersion;
for iwl=2:512;
wlength=[wlength, wlength(1)+(iwl-1)*dispersion];
end
  	 
	 %READ TIME CODE
	 
t=tsbase(shot,'GDATEVISR');
t=(t(:,1)*65536+t(:,2))*1e-6;
ignitron=tsmat(shot,'IGNITRON|1');
t=t-ignitron;

fp=1;
lp=250;
for ipix=1:(lp-fp+1);
spec(ipix,:,:)=spectra(ipix-1+fp,:,:);
wl(ipix)=wlength(ipix-1+fp);
pix(ipix)=ipix;
end;
for ik=1:length(t);
[h_dh(ik),h_dh_pks(ik)]=dh_ratio_disp(pix,reshape(spec(:,ik,iview),1,250),0);
end;
Sh_dh=transpose(smooth(transpose(h_dh),5));
Sh_dh_pks=transpose(smooth(transpose(h_dh_pks),5));

%figure(3);
%plot(t,h_dh,t,h_dh_pks);
%legend('algorithm','peaks');
% title([num2str(shot) '   H/H+D']);

H_HD=h_dh_pks;
H_HDcorr=h_dh;
tdata=t;

function [H_HD,H_HD_Peaks]=dh_ratio_disp(x,y,ipflg)


% ----------------
% x: array of lambda (or simply in pixel)
% y: array of counts . y(x) is the spectrum
%
% Output variables:
%------------------
% H_HD: is the isotopic ratio H/(H+D), where H and D are the intensities of
%       Dalpha and Halpha lines in counts. Of course it is 0 <= H_HD <= 1.
%
% spec_ok: this is a integer flag that specifies if the input 
%          spectrum is good for fitting:
%       spec_ok=1  spectrum ok, go and fit!
%       spec_ok=0  no signal over background limit: do not fit, 
%                  return zero curves plus mean bakground.
%       spec_ok=-1 bad spectrum (negative peaks, bground steps, 
%                  etc.): return zero curves plus background.
%
% fmin_ok: this is a integer flag that shows how the minimisation 
%          procedure has ended:
%       fmin_ok=1  minimun found.
%       fmin_ok=0  minimun not found, maximum number of iterations reached.
%
% fit_ok: this is a integer flag that shows if the best fit 
%         coefficients are phisically consistent:
%       fit_ok=1   fit coefficients Ok
%       fit_ok=0   problems... minimisation converged into a physically 
%                  insignificant solution.
%
%
%
%
% SUBFUNCTIONS:
% 1) A number of subfunctions are called for the execution of the main code:
%       a) check_spec
%       b) weight_res
%       c) Ferr_4g_8p_BL
%       d) check_4gfit
%       e) ii2ab
%    The subfunction  are all appended to this file, each with an 
%    explanation of their aim and specifications.
%==========================================================================


global ii_bground
ii_bground=[1:50]; % Background Range 
                        

global ii_xD
ii_xD=[150:170]; % Dalpha Range 
                        
                        
global dx_DH
dx_DH=19; % 'dx_DH' Distance (in pixel) between H and D

                         
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
global x_d y_d
x_d=x; y_d=y; % this is the experimental spectrum

global n_pixel
n_pixel=length(x_d); % the number of abscissae points (pixel)


global xDm yDm xHm yHm  % this are positions and heights of D and H peaks roughly determined
			% by 'check_spec'. 

global bground sd_bground % this are the constant background value (as a mean in a x_d range without
			  % Ha or Da signal) and his standard deviation, determined by 'check_spec'.
			   
global cres
cres=ones(1,n_pixel); % this is an array of weigth factors to calculate the 
                                   % weighted residuals in the minimisation procedure.
				   
global Hpro
Hpro=zeros(1,n_pixel);

global Dpro
Dpro=zeros(1,n_pixel);

global diffs neg pos;

Abf=0;         % this is the array of the 7 free minimisation parameters: 
			%3x2 (gauss)+1(bg).
			
H_HD=0;                % initial value of isotopic ratio.

Bbf=0;       % this is the array of the best fit coefficients (3x2 gaussians + background)

spec_ok=-100;          % default values of the others output variables
fmin_ok=-100;
fit_ok=-100;
%==========================================================================                                          
                                              

%==========================================================================
% Check Spectrum
%-----------------------------------------------------------------------
[xDm,yDm,xHm,yHm,bground,sd_bground,spec_ok]=check_spec(x_d,y_d);
%
% H/H+D for Peaks

H_HD_Peaks=yHm/(yHm+yDm);

%--------------------------------------------------------------------------

%if yDm<5000; H_HD=0; return; end;

if spec_ok==-1 % escape operations in case of bad experimental spectrum in input
   
    disp('   Bad experimental spectrum found!')
    %return
elseif spec_ok==0 % escape operations in case of too small line signal
    
    disp('   Bad experimental spectrum found!')    
    %return
end
%==========================================================================

%==========================================================================
% Set weights 
%-----------------------------------------------------------------------

weigth_res; 

%-----------------------------------------------------------------------          
options=optimset('Display','off','TolX',.1,'MaxIter',5000,'MaxFunEvals',100); 
%-----------------------------------------------------------------------

% Starting Fit Spectra


%Astart(1)=yHm/yDm;
%Astart(2)=xDm;
Astart=0;

%-----------------------------------------------------------------------
% Minimizing target (sub)function 'ferr_disp' we will get the
% best fit parameters.

%[Abf,fff,fmin_ok]=fminsearch('ferr_disp',Astart,options);
[Abf,fff,fmin_ok]=fminbnd('ferr_disp',0,1,options);



%==========================================================================
% From the 9 free parameters let's now calculate the 13 (12+1) parameters
% describing the 4 gaussian curves plus background, that is the 'Bbf'
% output variable.
%-----------------------------------------------------------------------

Bbf=Abf; 

%--------------------------------------------------------------------------
% Check for good fit, and identify components

%[Bbf(1:12),fit_ok]=check_4gfit(Bbf(1:12));

% y_fit=gauss1(B(1:3),x_d)+gauss1(B(4:6),x_d)+gauss1(B(7:9),x_d)+gauss1(B(10:12),x_d)+B(13);
% hold on, plot(x_d,y_fit,'r')
%--------------------------------------------------------------------------

%if fit_ok<0 % escape operations in case of bad fit
    
%    disp('   Bad resulting fit!')
%    %return
%end

%==========================================================================

%==========================================================================
% Calculate ratio
%--------------------------------------------------------------------------

%D=gauss_area(Bbf(1:3)); % intensity of Dalpha 
%H=gauss_area(Bbf(4:6)); % intensity of Halpha 
%H_HD=H/(H+D) % H/(T+H+D) intensity ratio
%fff
%Abf
if Abf~=0;
 H_HD=1/(1+1/Abf);
 else
 H_HD=0;
 end;

%------------------------------------------------------------------------
%
% Plot Spectra
%
%------------------------------------------------------------------------

 
gtot=Hpro+Dpro+bground;


if ipflg==1
%figure(2)
plot(x_d,Hpro,x_d,Dpro,x_d,gtot,x_d,y_d)
end;

%==========================================================================
% END OF MAIN
%==========================================================================


%--------------------------------------------------------------------------
function [x1m,y1m,x2m,y2m,bg,stdevbg,ok]=check_spec(x,y)
% It verifies if the input spectrum is ok for fitting with a number of tests.
% If all conditions are satisfied it returns a rough estimation of some 
% important quantities of the spectrum.
%
% Input variables:
% ----------------
% x: array of lambda (or simply in pixel)
% y: array of counts. y(x) is the spectrum
%
% Output variables:
%------------------
% x1m: position of first maximum
% y1m: height of first maximum
% x2m: position of second maximum
% y2m: height of second maximum
% bg:  background value
% devstbg: background fluctuation 
% ok:  exit flag
%       ok=1  spectrum ok for fitting.
%       ok=0  no signal over background limit. It returns  y1m=y2m=0.
%       ok=-1 bad spectrum (negative peaks, bground steps, 
%             etc.). It returns y1m=y2m=0.
%--------------------------------------------------

global ii_xD dx_DH ii_bground 

x1m=0;y1m=0;x2m=0;y2m=0;bg=0;devstbg=0;ok=1; % default output values

[y1m,i]=max(y([ii_xD]));

i1m=i+ii_xD(1)-1; x1m=x(i+ii_xD(1)-1);

i2m=i1m+dx_DH; y2m=y(i2m);

x2m=x(i2m);

bg=mean(y(ii_bground));

stdevbg=std(y(ii_bground));


%--------------------------------------------------
%if y1m<stdevbg,     y1m=1; end
%if y2m<stdevbg,     y2m=1; end
%if y1m<-5*stdevbg,  y1m=1; y2m=1; ok=-1,end
%if y2m<-5*stdevbg,  y2m=1; y1m=1; ok=-1,end
%--------------------------------------------------------------------------

function weigth_res();
% It provides the suitable array of factors for weigthing the resisuals
% in the minimisation procedure. It works on global variable so no
% input or output variable are needed.

global x_d y_d xDm yDm xHm yHm bground sd_bground cres

m=min(yDm,yHm);
M=max(yDm,yHm);
r=M/m;
ccc=M*5;%1e9 % default value

if r>3, ccc=bground+sd_bground*10; end

cres=1./(ccc+abs(y_d));
