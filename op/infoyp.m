function [] = infoyp(mode,a)
%
%	Customized message display
%
%	Input:
%
%		- mode: type of message [1,1]
%		- a: the message (string) [1,m]
%
%	Output: none
%
%
%by Y.PEYSSON CEA-DRFC 24/10/1991 <peysson@drfc.cad.cea.fr>
%revised for MatLab 5.2 (15/04/2000)
%
%
if nargin < 2,
	infoyp(2,'Wrong number of input arguments for infoyp');
	return;
end
%
if mode == 1,
	a1 = [];
	sa = size(a);
	nmax = 68;
	n0 = nmax-sa(2);
	n1 = fix(n0/2);
	n2 = n0-n1;
	for i = 1:n1,a1 = [a1,'-'];end
	if n2==n1,
 		a2 = a1;
	else
 		a2 = [a1,'-'];
	end
	disp([a1,' ',a,' ',a2]);
elseif mode == 2,
	disp(['-----> ',a]);
elseif mode == 3,
	ttt = [];
	sa = size(a);
	disp([a,'...']);
	for i = 1:sa(2)+10,
		ttt= [ttt,'-'];
	end
	disp(ttt);
elseif mode == 4,
	ttt = [];
	sa = size(a);
	disp([' - ',a]);
	disp(ttt);	
end
