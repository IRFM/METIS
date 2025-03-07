function [X,chi2] = regulyp(order,Y,dY,T,zrs,display,lambda,eps1,Xinit)
%
%       Regularization algorithm
%
%       Input:
%
%               - order: order of the regularization (0,1,2 or 3) [1,1]
%                 (default = 1)
%               - Y: normalized line-integrated emission [nchord,1]
%               - dY: relative error of the line-integrated emission [nchord,1]
%               - T: matrix of transfer between the local emissivity and the line
%                 integrated emission [nchord,nlayers]
%               - zrs: a fictive chord is added to ensure that cells which are not
%                 crossed by any chord do not contribute to the local emission [1,1]
%                 (default = 1) 
%               - display: plot the convergence in real-time (default = 0) [1,1]
%               - lambda: Lagrange multiplier (optional,default = 50)
%               - eps1: convergence parameter (optional,default = 0.02)    
%               - Xinit: initial guess for the local emissivity (optional, default = flat 
%                 profile) [1,nlayers]
%
%       Output: 
%
%               - X: local emissivity [1,nlayers]
%               - chi2: global error between Y and T*X (chi2<= 1 for an accurate 
%                 inversion) [1,1]
%
%
%by Y.PEYSSON CEA-DRFC 12/04/1995 <peysson@fedv09.cad.cea.fr>
%
%
if nargin < 3,infoyp(2,'Wrong number of input arguments for regulyp');return;end
if nargin == 3, order = 1, zrs = 1;display = 0; lambda = 50; eps1 = 0.02;end
if nargin == 4, zrs = 1;display = 0; lambda = 50;  eps1 = 0.02;end
if nargin == 5, display = 0; lambda = 50;  eps1 = 0.02;end
if nargin == 6, lambda = 50; eps1 = 0.02;end
if nargin == 7, eps1 = 0.02;end
%
%A fictive line of sight is added to remove the contribution of cells which are not crossed by any chord.
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
%------------------------  Regularisation --------------------------
%
if nargin < 9,
        Xinit = ones(sT2,1);
else
        Xinit = Xinit(:)/max(Xinit(:));
end             
%
x0 = Xinit;
W = sparse(diag(1.0./(dY.*dY)));
if order == 0,
        B = sparse(eye(sT2));
elseif order ==1,
        B = sparse(-diag(ones(sT2-1,1),1)-diag(ones(sT2-1,1),-1)+2*eye(sT2));
        B(1,1) = 1;
        B(sT2,sT2) = -1;
elseif order == 2,
        B = sparse(diag(ones(sT2-2,1),2) + diag(ones(sT2-2,1),-2)+6*eye(sT2) - 4*diag(ones(sT2-1,1),1) - 4*diag(ones(sT2-1,1),-1));
        B(1,1) = 1;    
        B(1,2) = -2;
        B(2,1) = -2;
        B(2,2) = 5;
        B(sT2-1,sT2-1) = 5;
        B(sT2-1,sT2) = -2;
        B(sT2,sT2-1) = -2;
        B(sT2,sT2) = 1;
elseif order == 3,
        B = sparse(6*diag(ones(sT2-2,1),2) + 6*diag(ones(sT2-2,1),-2)+20*eye(sT2) - 15*diag(ones(sT2-1,1),1) -...
                15*diag(ones(sT2-1,1),-1)-diag(ones(sT2-3,1),3) -diag(ones(sT2-3,1),-3));
        B(1,1) = 1;    
        B(1,2) = -3;
        B(1,3) = 3;
        B(2,1) = -3;
        B(2,2) = 10;
        B(2,3) = -12;
        B(3,1) = 3;
        B(3,2) = -12;
        B(3,3) = 19;
        B(sT2-2,sT2-2) = 19;
        B(sT2-2,sT2-1) = -12;
        B(sT2-2,sT2) = 3;
        B(sT2-1,sT2-2) = -12;
        B(sT2-1,sT2-1) = 10;
        B(sT2-1,sT2) = -3;
        B(sT2,sT2-2) = 3;
        B(sT2,sT2-1) = -3;
        B(sT2,sT2) = 1;
else
        error('Regularization order not valid !')
end
%
x0=inv(T'*W*T+lambda*B)*T'*W*Y;
%x0 = (T'*W*T+lambda*B)'\(T'*W*Y);
C0 = (Y-T*x0)'*W*(Y-T*x0)/sT1;
if display == 1, disp(['log(chi2) = ',num2str(log(C0))]);end
if C0 <= 0.1,
        if display == 1, infoyp(2,['Number of loops: ',int2str(1)]);end
        X = x0*maxY;
        chi2 = C0;
        return;
else
ii=0;
        for i = lambda:-2:1,
                x0 = inv(T'*W*T+i*B)*T'*W*Y;
%               x0 = (T'*W*T+lambda*B)'\(T'*W*Y);
                C0 = (Y-T*x0)'*W*(Y-T*x0)/sT1;
                if display == 1, disp(['log(chi2) = ',num2str(log(C0))]);end
                ii = ii+1;
                if C0-1 < eps1,
                        break
                end 
        end
        if display == 1, infoyp(2,['Number of loops: ',int2str(ii)]);end
        X = x0*maxY;
        chi2 = C0;
end





