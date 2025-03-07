function y = time()
% TIME		Return time in a string.
%function y = time()
%       S = TIME returns a string containing the time
%           in dd-mm-yy  hh:mm:ss format.
% R. MASSET   Aout 1992.

t = round(clock); an=int2str(t(1));
mois = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';...
        'Sep';'Oct';'Nov';'Dec'];

y = [int2str(t(3)),'-',mois(t(2),:),'-',an(3:4),'   ',...
     int2str(t(4)),':',int2str(t(5)),':',int2str(t(6))];
