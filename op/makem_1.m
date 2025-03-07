function M=makem_1(Y,T,xpix,ypix,flat)

% function M=makem_1(Y,T,xpix,ypix,flat) %
% set up flat default model with zero border % (used for different Xtomo algorithms, regulo_ , max_ent) %
%	input: Y	line int data
%	T	transfer matrix
%	xpix
%	ypix pixel coordinates
%	flat 1: totally flat model, just borders set zero
%	0: simple estimation used (see HOlland and Navratil)
%
%	output M	default model for the emissivities
%
%-------------MA 2/12/94 ------------------- 

if nargin==4
flat=1;
elseif nargin ~=4 & nargin ~=5
error('makem_1: wrong number of input arguments');
end;



epsilon=sqrt(eps);

[nl,npix]=size(T);
nx=length(xpix);
ny=length(ypix);
iborder=[1:ny,(1:(nx-2))*ny+1,(2:(nx-1))*ny,(1+(nx-1)*ny):nx*ny];
nborder=length(iborder);

% ---- data scaling

Ymax=max(Y);
D=max(Y/Ymax,epsilon*ones(nl,1));


% --- set default model --------------------- 

Tbar=sum(T(:))/(nl+npix);
Ep=D./(sum(T')');
nn=find(sum(T)~=0);
%Eca=dxml('*',Ep,T(:,nn),'T')./sum(T(:,nn));
Eca=Ep'*T(:,nn)./sum(T(:,nn));


%for c=1:npix
%	nn=find(T(:,c));
%	Ecg(c)=(prod(Ep(nn).*T(nn,c))).^(1/sum(T(nn,c)));
%end
%I0=sum(Eca+Ecg)/2;


if flat

% -- totally flat model with zero border 

I0=sum(Eca);
Emean=I0/(npix-nborder);
M=Emean*ones(npix,1);


else

M=Eca;

end

M(iborder)=1e-9*ones(size(iborder));
M=M*Ymax;
