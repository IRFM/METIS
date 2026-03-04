% FACIT interface with matlab
% --------------------------------------------------------------------------------
% out = FACIT_IO(in)
% Inputs
%

function out = FACIT_IO(in,tools_path);

FDIR = tools_path;

Ai = in.Ai;
Zi = in.Zi;
Aimp = in.Aimp;
Zimp = in.Zimp;
Te = in.Te;
Ti = in.Ti;
Ne = in.Ne;
Ni = in.Ni;
Nimp = in.Nimp;
x = in.x;
nx  = length(x);

if isfield(in,'gradTi')==1
	grad_ln_Ti = in.gradTi./Ti;
	gradTi = in.gradTi;
else
	grad_ln_Ti = in.grad_ln_Ti;
	gradTi = in.Ti.*in.grad_ln_Ti;
end
if isfield(in,'gradNi')==1
	grad_ln_ni = in.gradNi./Ni;
	gradNi = in.gradNi;
else
	grad_ln_ni = in.grad_ln_Ni;
	gradNi = in.Ni.*in.grad_ln_Ni;
end
if isfield(in,'gradNimp')==1
	grad_ln_na = in.gradNimp./Nimp;
	gradNimp = in.gradNimp;
else
	grad_ln_na = in.grad_ln_Nimp;
	gradNimp = in.Nimp.*in.grad_ln_Nimp;
end
if isfield(in,'dpsidx')==1
	dpsidx = in.dpsidx;
else
	dpsidx = 0*x;
end
if isfield(in,'FV')==1
	FV = in.FV;
else
	FV = in.B0.*in.R0*ones(size(in.x));
end
if isfield(in,'pol_asym')==1
	pol_asym = in.pol_asym;
else
	pol_asym = 1;
end
if isfield(in,'rotation')==1
	rotation = in.rotation;
else
	rotation = 1;
end
if isfield(in,'theta')==1
	theta = in.theta;
else
	theta = 0;
end
nth = length(theta);
qmag = in.qmag;
Machi = in.Machi;
invaspct = in.invaspct;
B0 = in.B0;
R0 = in.R0;
ionregime = in.ionregime;
geom = in.geom;
if geom == 1
	jacob = in.jacob;
	BV = in.BV;
	RV = in.RV;
	NV = in.NV;
	PhiV = in.PhiV;
	regulopt = in.param;
	AsymN = in.AsymN;
	AsymPhi = in.AsymPhi;
	if nth ~= 64
		disp('With in.geom == 1, nth should be equal to 64')
		stop
	end
else
	AsymN = in.AsymN;
	AsymPhi = in.AsymPhi;
	jacob = 0*x;
	BV = zeros(nx,nth);
	RV = zeros(nx,nth);
	NV = ones(nx,nth);
	PhiV = zeros(nx,nth);
	regulopt = zeros(4,1);
end
if length(Zimp)==1
	Zimp = Zimp*ones(size(x));
end

if ionregime == 1
 VS = zeros(1,max(18,nth));
 VS(1:15) = [nx, nth, Aimp,Ai,Zi,B0,R0,invaspct,pol_asym,geom,rotation,regulopt(1),regulopt(2),regulopt(3),regulopt(4)];
 VP = zeros(nx,max(18,nth));
 for ii=1:nx
	VP(ii,1:18) = [ x(ii),Zimp(ii),Te(ii),Ti(ii),Ne(ii),Ni(ii),Nimp(ii),qmag(ii),Machi(ii),gradNi(ii),gradTi(ii),gradNimp(ii),FV(ii),dpsidx(ii),AsymPhi(ii,1),AsymPhi(ii,2),AsymN(ii,1),AsymN(ii,2)];
 end
 VV = [VS;VP];
   if geom>0
	VF = zeros(1,max(18,nth));
	VFX = zeros(nx,max(18,nth));
	for k=1:nth
		VF(1,k) = [theta(k)];
	end
	VFX(:,1:nth) = [BV];
	VF = [VF;VFX];
	VFX(:,1:nth) = [RV];
	VF = [VF;VFX];
	VFX(:,1:nth) = [jacob];
	VF = [VF;VFX];
	VFX(:,1:nth) = [PhiV];
	VF = [VF;VFX];
	VFX(:,1:nth) = [NV];
	VF = [VF;VFX];
	VV = [VV;VF];
 end

 eval(['save facit_input.dat VV -ASCII'])
 
%Check FACIT interface is present in current directory
currentpath = pwd;
filename = [pwd '/FACIT_interface'];
if exist(filename, 'file') ~= 2
	eval(['!cp ',FDIR,'/runs/FACIT_interface .'])
end

disp(['==============        RUNS FACIT               =================='])
[sout,tout] = unix(sprintf(['module purge \n module load inteloneapi/2025.1.0 lapack/3.9.1-inteloneapi-2025.1.0 \n ./FACIT_interface']));
%disp(tout);
disp('==============        END OF FACIT RUN        ==================')

 %disp(['Read facit_output file'])
 eval(['load facit_output.dat -ASCII'])
 offset = 20;
 Va = facit_output(1:nx,offset);
 Da = facit_output(1:nx,offset+1);
 dmin = facit_output(1:nx,offset+2);
 dmaj = facit_output(1:nx,offset+3);
 Da_PS = facit_output(1:nx,offset+4);
 Da_BP = facit_output(1:nx,offset+5);
 Da_cl = facit_output(1:nx,offset+6);
 Ka_PS = facit_output(1:nx,offset+7);
 Ka_BP = facit_output(1:nx,offset+8);
 Ka_CL = facit_output(1:nx,offset+9);
 Ha_PS = facit_output(1:nx,offset+10);
 Ha_BP = facit_output(1:nx,offset+11);
 Ha_CL = facit_output(1:nx,offset+12);
 Va_PS = facit_output(1:nx,offset+13);
 Va_BP = facit_output(1:nx,offset+14);
 Va_cl = facit_output(1:nx,offset+15);
 Da_neo = Da - Da_cl;
 Va_neo = Va - Va_cl;
 timeexec = facit_output(1:nx,offset+16);
 if geom == 1
	nn = facit_output(nx+1:2*nx,1:nth);
 end

 if geom == 0
	out = struct('rhon',x, ...
	'dmin',dmin, ...
	'dmaj',dmaj, ...
	'Da',Da, ...
	'Da_neo',Da_neo, ...
	'Da_cl',Da_cl, ...
	'Da_PS',Da_PS, ...
	'Da_BP',Da_BP, ...
	'Va',Va, ...
	'Va_neo',Va_neo, ...
	'Va_cl',Va_cl, ...
	'Va_PS',Va_PS, ...
	'Va_BP',Va_BP, ...
	'Ka_PS',Ka_PS, ...
	'Ka_BP',Ka_BP, ...
	'Ka_CL',Ka_CL, ...
	'Ha_PS',Ha_PS, ...
	'Ha_BP',Ha_BP, ...
	'Ha_CL',Ha_CL, ...
	'info',timeexec);
else
	out = struct('rhon',x, ...
	'dmin',dmin, ...
	'dmaj',dmaj, ...
	'Da',Da, ...
	'Da_neo',Da_neo, ...
	'Da_cl',Da_cl, ...
	'Da_PS',Da_PS, ...
	'Da_BP',Da_BP, ...
	'Va',Va, ...
	'Va_neo',Va_neo, ...
	'Va_cl',Va_cl, ...
	'Va_PS',Va_PS, ...
	'Va_BP',Va_BP, ...
	'Ka_PS',Ka_PS, ...
	'Ka_BP',Ka_BP, ...
	'Ka_CL',Ka_CL, ...
	'Ha_PS',Ha_PS, ...
	'Ha_BP',Ha_BP, ...
	'Ha_CL',Ha_CL, ...
	'nn',nn, ...
	'info',timeexec);
end

else
	out = FACIT(in);
end

