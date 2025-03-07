function [X,chi2,sort] = maxentyp2(mode,Y,dY,T,nevrai)
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
zrs = 1;display = 1; 

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
%-------------------------  Davidon-Fletcher-Powell  ----------------------------
%


		lambda = 0.000000001;
		eps1 = 0.005;
		eps2 = 0.005;
		gamma = 1.5;
		Xinit = ones(sT2,1);
%
%
        lam=logspace(-4,0,100);
%        lam=1;
%	while C0-1>eps1,
	ii = 0;
	j=1:100; 
        for kl=1:length(lam)
	  x0 = Xinit;
	  W=diag(1.0./(dY.*dY));
	  A0=sum(x0);
	  B=eye(sT2,sT2);
	  q1=ones(100,1)./((2*ones(100,1))'.^j)';
%
%
	  S0        = -x0'*log(x0/A0)/A0;
	  C0        = (Y-T*x0)'*W*(Y-T*x0)/sT1;
 	  ii        = ii+1; 
          norm1     = 1;
          lambdacur = lam(kl);
    	  while (norm1 > eps2),
            g0      = (S0+log(x0/A0))/A0+2*lambdacur*T'*W*(T*x0-Y);
            h       = (-(2*S0+1)*ones(sT2,1)*ones(sT2,1)'-log((x0/A0)*(x0/A0)'))/(A0^2);
            hh      = h+diag(ones(sT2,1)./x0,0)/A0+2*lambdacur*T'*W*T;
            p       = B*g0;
            alfa    = (p'*g0)/(p'*hh*p);
            x2      = x0-alfa*p;
            if (min(x2)<=0),
            	xm     = x0*ones(100,1)'-(alfa*ones(sT2,1)*q1').*(p*ones(100,1)');
              	omega  = max(q1(find(min(xm)>0)));
              	x2     = x0-omega*alfa*p;
            end
            A2       = sum(x2); 
            S        = -x2'*log(x2/A2)/A2;
            g        = (S+log(x2/A2))/A2+2*lambdacur*T'*W*(T*x2-Y);
            B        = B+((x2-x0)*(x2-x0)')/((x2-x0)'*(g-g0))-((B*(g-g0))*(B*(g-g0))')/((g-g0)'*B*(g-g0));
            norm1    = norm(x0-x2)/norm(x2);
            x0       = x2; 
            A0       = sum(x0);
            S0       = -x0'*log(x0/A0)/A0;
    	  end
 	  C1         = (Y-T*x0)'*W*(Y-T*x0)/sT1;
 	  if C1 > C0,
     	  disp('Unsuccessfull convergence !!!!')
%     	  break;
 	  else
   	          C0 = C1;
 	  end
	  if display == 1, disp(['log(chi2) = ',num2str(log(C0))]);end
 	  gradS        = -(S0+log(x0/A0))/A0;
 	  gradC        = 2*T'*W*(T*x0-Y);
 	  lambdaopt    = (gradS'*gradS)/(gradS'*gradC);
% 	  lambdacur = gamma*abs(lambdaopt);
          sort.C1(kl) = C1;
          sort.g(kl,:) = g0;
          sort.S(kl) = S;
          sort.A2(kl) = A0;
          sort.f(kl,:) = x0;
          sort.lam(kl)=lambdacur;
          sort.lamopt(kl)=lambdaopt;
          sort.norm(kl)=norm(x2-nevrai)/norm(nevrai);
	end
	if display == 1, infoyp(2,['Number of loops: ',int2str(ii)]);end
	X=x0*maxY;
	chi2 = C0;




