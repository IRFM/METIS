% JULDAY Calculate the Julian Day Number for a given year, month, and day.
% ------------------------------------------------------------------------
% This is the inverse of the library function CALDAT.

% Adapted from "Numerical Recipes in C",
% by Press, Teukolsky, Vetterling, and Flannery.
% Cambridge University Press, 1992 (2nd Edition).
% ----------------------------------------------------------------------

function julian =julday(year, month, day)


if nargin ~= 3
  error('Usage: [julian]=julday(year, month, day) (integers)');
end

igreg=15+31*(10+12*1582);       % Gregorian Calander adopted Oct. 15, 1582

test=year<=0;
year(test)=year(test)+1;

test=month>2;
jy=year;jm=month;
jm(test)=jm(test)+1;
jy(~test)=jy(~test)-1;
jm(~test)=jm(~test)+13;

julian=fix(floor(365.25*jy)+floor(30.6001*jm)+day+1720995);

% Test whether to change to Gregorian Calandar
test=day+31*(month+12*year)>=igreg;
ja=fix(0.01*jy(test));
julian(test)=julian(test)+2-ja+fix(0.25*ja);
