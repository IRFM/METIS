function [X,chi2] = maxentyp(mode,Y,dY,T,zrs,display,lambda,eps1,eps2,gamma,Xinit)
%
%	Maximum entropy algorithm
%
%	Input:
%
%		- mode: maximum entropy algorithm 	[1,1]
%					- (1) Paraboloid
%					- (2) Davidon-Fletcher-Powell
%					- (3) Gull-Daniell-Delsuc
%					- (4) Gull-Daniell-Wu
%		- Y: normalized line-integrated emission [nchord,1]
%		- dY: relative error of the line-integrated emission [nchord,1]
%		- T: matrix of transfer between the local emissivity and the line
%		  integrated emission [nchord,nlayers]
%		- zrs: a fictive chord is added to ensure that cells which are not
%		  crossed by any chord do not contribute to the local emission [1,1]
%		  (optional,default = 1) 
%		- display: display the convergence in real-time (default = 0) [1,1]
%		- eps1: convergence parameter (optional,default = 0.005)    
%		- eps2: convergence parameter (optional,default = 0.005) 
%		- gamma: convergence parameter (optional,default = 2)
%		- lambda: Lagrange multiplier (optional,default = 1e-9)
%		- Xinit: initial guess for the local emissivity (optional, default = flat 
%		  profile) [1,nlayers]
%
%	Output: 
%
%		- X: local emissivity [1,nlayers]
%		- chi2: global error between Y and T*X (chi2<= 1 for an accurate 
%		  inversion) [1,1]
%
%
%by Y.PEYSSON CEA-DRFC 13/07/1994 <peysson@fedv09.cad.cea.fr>
%
%
if nargin < 4,infoyp(2,'Wrong number of input arguments for maxentyp');return;end
if nargin == 4, zrs = 1;display = 0; end
if nargin == 5, display = 0; end
%
%A fictive line of sight is added to remove the contribution of cells which are not
%crossed by any chord.
%
if zrs == 1,
	S = (T > 0);
	W = sum(S);
	WW = find(W<1);
	T0 = zeros(1,size(T,2));
	T0(WW) = T0(WW)+1;
	T = [T;T0];
	Y = [Y;0];
	dY = [dY;0.0000001];
end
%
maxY = max(Y);
Y = Y/max(Y);
[sT1,sT2] = size(T);
%
%------------------------  Paraboloid  ----------------------------
%
if mode == 1,
	if nargin < 7,
		lambda = 0.000000001;
		eps1 = 0.005;
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 7,
		eps1 = 0.005;
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);		
	elseif nargin == 8,
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 9,
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 10,
		Xinit = ones(sT2,1);
	elseif nargin == 11,
		Xinit = Xinit(:)/max(Xinit(:));
	end
%
	x0 = Xinit;
	W=sparse(diag(1.0./(dY.*dY)));
	A0=sum(x0);
	lambdacur0 = lambda;
	ii=0;
	z=eps*ones(size(x0)); 
%
	S0=-x0'*log(x0/A0)/A0;
	C0=(Y-T*x0)'*W*(Y-T*x0)/sT1;
%
	while C0-1>eps1,
   		ii=ii+1; 
   		g0 = (S0+log(x0/A0))/A0+2*lambdacur0*T'*W*(T*x0-Y);
   		h = (-(2*S0+1)*ones(sT2,1)*ones(sT2,1)'-log((x0/A0)*(x0/A0)'))/(A0^2);
   		hh = h+diag(ones(sT2,1)./x0,0)/A0+2*lambdacur0*T'*W*T;
   		p=hh\g0;
   		x2=x0-p;
   		if (min(x2)<=0),
      		i=find(x2<=0);
      		x2(i)=z(i);
   		end
   		x0=x2;
   		A0=sum(x0);
   		S0=-x0'*log(x0/A0)/A0;
   		C1=(Y-T*x0)'*W*(Y-T*x0)/sT1;
   		if C1 > C0,
     		disp('Unsuccessfull convergence !!!!')
     		break;
   		else
     		C0 = C1;
   		end
  		if display == 1, disp(['log(chi2) = ',num2str(log(C0))]);end
  		gradS=-(S0+log(x0/A0))/A0;
  		gradC=2*T'*W*(T*x0-Y);
  		lambdaopt = (gradS'*gradS)/(gradS'*gradC);
   		lambdacur0 = gamma*abs(lambdaopt);
	end
	if display == 1, infoyp(2,['Number of loops: ',int2str(ii)]);end
	X=x0*maxY;
	chi2 = C0;
%
%-------------------------  Davidon-Fletcher-Powell  ----------------------------
%
elseif mode == 2,
	if nargin < 7,
		lambda = 0.000000001;
		eps1 = 0.005;
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 7,
		eps1 = 0.005;
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);		
	elseif nargin == 8,
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 9,
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 10,
		Xinit = ones(sT2,1);
	elseif nargin == 11,
		Xinit = Xinit(:)/max(Xinit(:));
	end
%
	x0 = Xinit;
	W=diag(1.0./(dY.*dY));
	A0=sum(x0);
	B=eye(sT2,sT2);
	j=1:1:100; 
	q1=ones(100,1)./((2*ones(100,1))'.^j)';
%
	ii = 0;
	lambdacur = lambda;  
%
	S0=-x0'*log(x0/A0)/A0;
	C0=(Y-T*x0)'*W*(Y-T*x0)/sT1;
%
	while C0-1>eps1,
 		ii=ii+1; norm1 = 1;
    	while (norm1 > eps2),
            g0 = (S0+log(x0/A0))/A0+2*lambdacur*T'*W*(T*x0-Y);
            h = (-(2*S0+1)*ones(sT2,1)*ones(sT2,1)'-log((x0/A0)*(x0/A0)'))/(A0^2);
            hh = h+diag(ones(sT2,1)./x0,0)/A0+2*lambdacur*T'*W*T;
            p=B*g0;
            alfa = (p'*g0)/(p'*hh*p);
            x2 = x0-alfa*p;
            if (min(x2)<=0),
            	xm=x0*ones(100,1)'-(alfa*ones(sT2,1)*q1').*(p*ones(100,1)');
              	omega=max(q1(find(min(xm)>0)));
              	x2=x0-omega*alfa*p;
           	end
          	A2=sum(x2); S=-x2'*log(x2/A2)/A2;
          	g= (S+log(x2/A2))/A2+2*lambdacur*T'*W*(T*x2-Y);
          	B=B+((x2-x0)*(x2-x0)')/((x2-x0)'*(g-g0))-((B*(g-g0))*(B*(g-g0))')/((g-g0)'*B*(g-g0));
          	norm1=norm(x0-x2)/norm(x2);
          	x0=x2; A0=sum(x0);
          	S0=-x0'*log(x0/A0)/A0;
    	end
 		C1=(Y-T*x0)'*W*(Y-T*x0)/sT1;
 		if C1 > C0,
     		disp('Unsuccessfull convergence !!!!')
     		break;
 		else
   			C0 = C1;
 		end
		if display == 1, disp(['log(chi2) = ',num2str(log(C0)), 'lambda=',num2str(lambdacur,4)]);end
 		gradS=-(S0+log(x0/A0))/A0;
 		gradC=2*T'*W*(T*x0-Y);
 		lambdaopt = (gradS'*gradS)/(gradS'*gradC);
 		lambdacur = gamma*abs(lambdaopt);
	end
	if display == 1, infoyp(2,['Number of loops: ',int2str(ii)]);end
	X=x0*maxY;
	chi2 = C0;
%
%--------------------------  Gull-Daniell-Delsuc  ------------------------------
%					

elseif mode == 3,
	if nargin < 7,
		lambda = 0.000000001;
		eps1 = 0.005;
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 7,
		eps1 = 0.005;
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);		
	elseif nargin == 8,
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 9,
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 10,
		Xinit = ones(sT2,1);
	elseif nargin == 11,
		Xinit = Xinit(:)/max(Xinit(:));
	end
%
	x0 = Xinit;
	A0=sum(x0);
	W=diag(1.0./(dY.*dY));
%
	S0=-x0'*log(x0/A0)/A0;
	C0=(Y-T*x0)'*W*(Y-T*x0)/sT1;
	alerte=0;
	sbeta=10;betammax=1;
%
	lambdacur0 = lambda;
	ii=0;
%
	lambdacur=lambdacur0;
	while C0-1>eps1,
		ii=ii+1;
 		testproba = 1;
 		while testproba>eps2,
  			x1=A0*exp(-S0-2*lambdacur*A0*T'*W*(T*x0-Y));
%                        keyboard
  			betam=linspace(0,betammax,sbeta);
  			for i=1:sbeta,
   				xm=(1-betam(1,i))*x0+betam(1,i)*x1;Am=sum(xm);
   				Sm(1,i)=-xm'*log(xm/Am)/Am;
   				Cm(1,i)=(Y-T*xm)'*W*(Y-T*xm);
  			end
  			Qm=Sm-lambdacur*Cm;
  			beta=betam(find(Qm==max(Qm)));
  			betammax=beta
  			x2=(1-beta)*x0+beta*x1;A2=sum(x2);
  			testproba1=norm(x2-x0)/norm(x0);
  			if testproba1>testproba,
   				alerte=1;
   				break;
  			else
   				testproba=testproba1;
  			end
  			x0=x2;A0=A2;S0=-x0'*log(x0/A0)/A0;
 		end
 		C1=(Y-T*x0)'*W*(Y-T*x0)/sT1
 		if C1 > C0,
     		disp('Unsuccessfull convergence !!!!')
     		break;  			
		else
  			C0 = C1;
 		end
		if display == 1, disp(['log(chi2) = ',num2str(log(C0))]);end
 		gradS=-(S0+log(x0/A0))/A0;
 		gradC=2*T'*W*(T*x0-Y);
 		lambdaopt = (gradS'*gradS)/(gradS'*gradC);
 		lambdacur = gamma*lambdaopt;
                eval(['save maxent',int2str(ii)])
	end
	if display == 1, infoyp(2,['Number of loops: ',int2str(ii)]);end
	X=x0*maxY;
	chi2 = C0;
%
%----------------------------  Gull-Daniell-Wu  ------------------------------
%		

elseif mode == 4,
	if nargin < 7,
		lambda = 0.000000001;
		eps1 = 0.005;
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 7,
		eps1 = 0.005;
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);		
	elseif nargin == 8,
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 9,
		gamma = 1.5;
		Xinit = ones(sT2,1);
	elseif nargin == 10,
		Xinit = ones(sT2,1);
	elseif nargin == 11,
		Xinit = Xinit(:)/max(Xinit(:));
	end
%
	w0 = sum(2*(T.*T)./((dY.^2)*ones(1,sT2)));
	x0 = Xinit';
	x4 = x0;
	alerte = 1;
	gradC = -2*sum(T.*((-(Y-T*x0')./(dY.^2))*ones(1,sT2)));
	S = -(x0/sum(x0))*log((x0/sum(x0))');
	C0 = sum(((Y-T*x0').^2)./(dY.^2))/sT1;
%
	lambdacur = lambda;
	beta = 0.1;
	ii=0;
%
	while (C0-1 > eps1),
 		ii=ii+1;
 		testproba = 1;
 		while (testproba > eps2),
  			x1 = sum(x0)*exp(-S+lambdacur*gradC*sum(x0));
  			x2 = ((1+lambdacur*(w0.*x0)/sum(x0))./(1+lambdacur*(w0.*x1)/sum(x1))).*x1;
  			x3 = (1-beta)*x0+beta*x2;
  			testproba1 = norm(x3-x0)/norm(x0);
  			if (testproba1 > testproba),
   				alerte=1;
   				break;
  			else
   				testproba=testproba1;
  			end
  			x0=x3;
  			gradC=-2*sum(T.*((-(Y-T*x0')./(dY.^2))*ones(1,sT2)));
  			S=-(x0/sum(x0))*log((x0/sum(x0))');
 		end
 		if alerte==1,
  			beta=beta*0.5;
  			x0=x4;
  			gradC=-2*sum(T.*((-(Y-T*x0')./(dY.^2))*ones(1,sT2)));
  			gradS = -(log(x0/sum(x0))+S)/sum(x0);
  			alerte=0;
 		else
  			x4=x0;
  			C1 = sum(((Y-T*x0').^2)./(dY.^2))/sT1;
  			if C1 > C0,
			    disp('Unsuccessfull convergence !!!!')
   				break;
  			else
   				C0 = C1;
  			end
			if display == 1, disp(['log(chi2) = ',num2str(log(C0))]);end
  			gradS = -(log(x0/sum(x0))+S)/sum(x0);
  			lambdaopt = (gradS*gradS')/(gradS*gradC');
  			lambdacur = gamma*abs(lambdaopt);
 		end
	end
	if display == 1, infoyp(2,['Number of loops: ',int2str(ii)]);end
	X=x0'*maxY;
	chi2 = C0;
else
	infoyp(2,'Maximum entropy algorithm unknown');
	return;
end




