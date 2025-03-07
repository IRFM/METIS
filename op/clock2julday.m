% retourne la date en secondes juliennes
function julian = clock2julday




% JULDAY Calculate the Julian Day Number for a given year, month, and day.
% ------------------------------------------------------------------------
% This is the inverse of the library function CALDAT.

% Adapted from "Numerical Recipes in C",
% by Press, Teukolsky, Vetterling, and Flannery.
% Cambridge University Press, 1992 (2nd Edition).
% ----------------------------------------------------------------------
%function julian =julday(year, month, day)
% CLOCK = [year month day hour minute seconds]
ckl = clock;
ck = local_date_vector_to_utc(ckl);
ck(6) = ckl(6);
year = ck(1);
month = ck(2);
day = ck(3);

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


julian = julian .* 3600 .* 24  + ck(4) .* 3600 + ck(5) .* 60  + ck(6); 



function dateVectorUtc = local_date_vector_to_utc( dateVectorLocal )
      
  nDates = size(dateVectorLocal,1) ;
  dateVectorUtc = zeros(nDates,6) ;
  
% import the Java classes needed for this process

  import java.text.SimpleDateFormat ;
  import java.util.Date ;
  import java.util.TimeZone ;
  
% instantiate a SimpleDateFormat object with a fixed time/date format and UTC time zone

  utcFormatObject = SimpleDateFormat('yyyy-MM-dd HH:mm:ss') ;
  utcFormatObject.setTimeZone(TimeZone.getTimeZone('UTC')) ;
      
% loop over date strings

  for iDate = 1:nDates

      dateVec = dateVectorLocal(iDate,:) ;
      
%     instantiate a Java Date class object with the local time.  Note that Java year is
%     year since 1900, and Java month is zero-based
      
      localDateObject = Date(dateVec(1)-1900, dateVec(2)-1, dateVec(3), ...
                             dateVec(4), dateVec(5), dateVec(6)) ;                        
                         
%     convert the date object to a string in the correct format and in UTC

      dateStringUtc = char(utcFormatObject.format(localDateObject)) ;
         
%     pick through the resulting string and extract the data we want, converting to
%     numbers as we go

      dateVectorUtc(iDate,1) = str2num(dateStringUtc(1:4)) ;
      dateVectorUtc(iDate,2) = str2num(dateStringUtc(6:7)) ;
      dateVectorUtc(iDate,3) = str2num(dateStringUtc(9:10)) ;
      dateVectorUtc(iDate,4) = str2num(dateStringUtc(12:13)) ;
      dateVectorUtc(iDate,5) = str2num(dateStringUtc(15:16)) ;
      dateVectorUtc(iDate,6) = str2num(dateStringUtc(18:19)) ;
          
  end % loop over dates
