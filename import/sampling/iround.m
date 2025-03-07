function K = iround(X,X0)
 
% IROUND	index table lookup
% function k = iround(x,x0)
%
% index table lookup
%
% x   monotonic increasing table
% x0  values to look up
% k   index in x of the nearest of x0
% @(#)iround.m	1.4	J.-M. Moret MatLab ToolBox	8/29/90
% version matlab pour CRONOS : ne pas comiler iround.c

if isempty(X) | isempty(X0)
	K = [];
elseif length(X) == 1
	K = 1;
else
	K = interp1(X,1:length(X),X0,'nearest');
end



