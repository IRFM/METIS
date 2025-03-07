function [A,z]=jacobiancsd(fun,x)
% JACOBIANCSD    Complex Step Jacobian
% J = jacobiancsd(f,x) returns the numberical (m x n) Jacobian matrix of a 
% m-vector function, f(x) at the refdrence point, x (n-vector).
% [J,z] = jacobiancsd(f,x) also returns the function value, z=f(x).
%
% Example
% f=@(x)[x(1)*x(2)^2;x(1)*x(3)^2;x(1)^2];
% x=1:3;
% [J,z]=jacobiancsd(f,x)
% J =
%     4     4     0
%     9     0     6
%     2     0     0
% z =
%     4
%     9
%     1
%
% By Yi Cao at Cranfield University, 02/01/2008
%
z=fun(x);                       % function value
n=numel(x);                     % size of independent
m=numel(z);                     % size of dependent
A=zeros(m,n);                   % allocate memory for the Jacobian matrix
h=n*eps;                        % differentiation step size
for k=1:n                       % loop for each independent variable 
    x1=x;                       % reference point
    x1(k)=x1(k)+h*i;            % increment in kth independent variable
    A(:,k)=imag(fun(x1))/h;     % complex step differentiation 
end

% complexe step derivative
%  Copyright (c) 2009, Yi Cao
%  All rights reserved.
%  
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%  
%      * Redistributions of source code must retain the above copyright 
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright 
%        notice, this list of conditions and the following disclaimer in 
%        the documentation and/or other materials provided with the distribution
%        
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
