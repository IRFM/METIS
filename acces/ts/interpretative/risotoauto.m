%  RISOTOAUTO  courte description  
%------------------------------------------------------------------------------- 
% fichier :  risotoauto.m  ->  risotoauto 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   [times,ti_h,ti_d,nhnd,va,certout]=risotoauto(shot,E) 
%  
% entrees :  
%  shot = 
%  E    = 
%  
% sorties :  
%  times   = 
%  ti_h    = 
%  ti_d    = 
%  nhnd    = 
%  va      = 
%  certout = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  4.1  du  08/04/2008  
%  
%*@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function [times,ti_h,ti_d,nhnd,va,certout]=risotoauto(shot,E)

times =[];
ti_h =[];
ti_d =[];
nhnd =[];
va =[];
certout=[];

%[times,ti_h,ti_d,nhnd,val]=risotoauto(shot,E)
%
% shot : numero de choc
% E    : energie maximum de DCX en keV
%
tipsh1=[];tipsd1=[];ninh1=[];nind1=[];
tipsh2=[];tipsd2=[];ninh2=[];nind2=[];
tipsh3=[];tipsd3=[];ninh3=[];nind3=[];
tipsh4=[];tipsd4=[];ninh4=[];nind4=[];
tipsh5=[];tipsd5=[];ninh5=[];nind5=[];
tipsh6=[];tipsd6=[];ninh6=[];nind6=[];
times1 = [];
times2 = [];
times3 = [];
times4 = [];
times5 = [];
times6 = [];

nhnd = [];
choc=shot;

[fci,tfci,cert]   = tsbase(choc,'spuiss');

[cxh1,t,eh1,cert] = tsbase(choc,'gh1');
[cxd1,t,ed1,cert] = tsbase(choc,'gd1');
if isempty(cxh1)
  disp('ANR1 HS'); 
else
  times1=t;
  certout =cert;
end
  
[cxh2,t,eh2,cert]=tsbase(choc,'gh2');
[cxd2,t,ed2,cert]=tsbase(choc,'gd2');
  if isempty(cxh2), disp('ANR2 HS'); else, times2=t;certout =cert; end   

[cxh3,t,eh3,cert]=tsbase(choc,'gh3');
[cxd3,t,ed3,cert]=tsbase(choc,'gd3');
  if isempty(cxh3), disp('ANR3 HS'); else, times3=t;certout =cert; end   
  
[cxh4,t,eh4,cert]=tsbase(choc,'gh4');
[cxd4,t,ed4,cert]=tsbase(choc,'gd4');
  if isempty(cxh4), disp('ANR4 HS'); else, times4=t;certout =cert; end   

[cxh5,t,eh5,cert]=tsbase(choc,'gh5');
[cxd5,t,ed5,cert]=tsbase(choc,'gd5');
  if isempty(cxh5), disp('ANR5 HS'); else, times5=t;certout =cert; end   

	
[cxh6,t,eh6,cert]=tsbase(choc,'gh6');
[cxd6,t,ed6,cert]=tsbase(choc,'gd6');
if isempty(cxh6)
    disp('ANR6 HS'); 
else
    times6=t; 
    certout =cert;
end
  
ltimes = [length(times1) length(times2) length(times3) length(times4) length(times5) length(times6)];


if sum(ltimes) == 0
  val = 0;
  return
end
if ltimes(1) > 0
  times = times1(1:min(ltimes(ltimes>0)));
elseif ltimes(2) > 0
  times = times2(1:min(ltimes(ltimes>0)));
elseif ltimes(3) > 0	
  times = times3(1:min(ltimes(ltimes>0)));	
elseif ltimes(4) > 0	
  times = times4(1:min(ltimes(ltimes>0)));  
elseif ltimes(5) > 0	
  times = times5(1:min(ltimes(ltimes>0)));  
elseif ltimes(6) > 0	
  times = times6(1:min(ltimes(ltimes>0)));  
end  
  


val = 1;
temps=5;
stemps=[', t=' num2str(floor(temps*100)./100) 's'];
xx=iround(times,temps);
%
% choix de l'energie maximum pour l'analyse (en keV)
%
emax     = E;

	if ~isempty(cxh1),
	
	cmax     = iround(eh1,emax);
	cinth    = 1:cmax;
	cmax     = iround(ed1,emax);
	cintd    = 1:cmax;	
	cxh11    = cinth(isfinite(cxh1(xx,cinth)./sqrt(eh1(cinth))));
	cxd11    = cintd(isfinite(cxd1(xx,cintd)./sqrt(ed1(cintd))));

	end
	
	
	if ~isempty(cxh2),
	
	cmax     = iround(eh2,emax);
	cinth    = 1:cmax;
	cmax     = iround(ed2,emax);
	cintd    = 1:cmax;	
	cxh22    = cinth(isfinite(cxh2(xx,cinth)./sqrt(eh2(cinth))));
	cxd22    = cintd(isfinite(cxd2(xx,cintd)./sqrt(ed2(cintd))));
	end

	if ~isempty(cxh3),
	
	cmax     = iround(eh3,emax);
	cinth    = 1:cmax;
	cmax     = iround(ed3,emax);
	cintd    = 1:cmax;	
	cxh33    = cinth(isfinite(cxh3(xx,cinth)./sqrt(eh3(cinth))));
	cxd33    = cintd(isfinite(cxd3(xx,cintd)./sqrt(ed3(cintd))));

	end

	if ~isempty(cxh4),
	
	cmax     = iround(eh4,emax);
	cinth    = 1:cmax;
	cmax     = iround(ed4,emax);
	cintd    = 1:cmax;	
	cxh44    = cinth(isfinite(cxh4(xx,cinth)./sqrt(eh4(cinth))));
	cxd44    = cintd(isfinite(cxd4(xx,cintd)./sqrt(ed4(cintd))));
	end

	if ~isempty(cxh5),
	
	cmax     = iround(eh5,emax);
	cinth    = 1:cmax;
	cmax     = iround(ed5,emax);
	cintd    = 1:cmax;	
	cxh55    = cinth(isfinite(cxh5(xx,cinth)./sqrt(eh5(cinth))));
	cxd55    = cintd(isfinite(cxd5(xx,cintd)./sqrt(ed5(cintd))));

	end
	
	if ~isempty(cxh6),
	
	cmax     = iround(eh6,emax);
	cinth    = 1:cmax;
	cmax     = iround(ed6,emax);
	cintd    = 1:cmax;	
	cxh66    = cinth(isfinite(cxh6(xx,cinth)./sqrt(eh6(cinth))));
	cxd66    = cintd(isfinite(cxd6(xx,cintd)./sqrt(ed6(cintd))));
	end
	
% Cette version suppose que atte=1 et sigmavcx = 1	
% masse ions = 0 : sigmavcx=1;
% attenuation par cx, ionisation, atte = exp(-lamdba) suppose =1
% sigmacx supppose =1

mh=0;md=0; 
atte=1;

  for ll=1:length(times),


% ANR1
    if ~isempty(cxh1),
		flux=cxh1(ll,cxh11); ener=eh1(cxh11);
		sigmavcx=sigcx(ener,mh);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energh1=x;
     	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end 
		tipsh1=[tipsh1 -1/cp(1)]; ninh1=[ninh1 (-1/cp(1))^1.5*exp(cp(2))];
     	else
		tipsh1=[tipsh1 0]; ninh1=[ninh1 0];
     	end
    else
        tipsh1=[tipsh1 0]; ninh1=[ninh1 0]; energh1=[];
    end
	
 
    if ~isempty(cxd1),
		flux=cxd1(ll,cxd11); ener=ed1(cxd11);
		sigmavcx=sigcx(ener,md);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energd1=x;
    	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsd1=[tipsd1 -1/cp(1)]; nind1=[nind1 (-1/cp(1))^1.5*exp(cp(2))];
    	else
		tipsd1=[tipsd1 0]; nind1=[nind1 0];
    	end
    else
    	tipsd1=[tipsd1 0]; nind1=[nind1 0]; energd1=[];
    end


% ANR2

    if ~isempty(cxh2),
		flux=cxh2(ll,cxh22);ener=eh2(cxh22);
		sigmavcx=sigcx(ener,mh);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energh2=x;
     	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsh2=[tipsh2 -1/cp(1)]; ninh2=[ninh2 (-1/cp(1))^1.5*exp(cp(2))];
     	else
		tipsh2=[tipsh2 0]; ninh2=[ninh2 0];
     	end
    else
        tipsh2=[tipsh2 0]; ninh2=[ninh2 0];energh2=[];
    end
	

    if ~isempty(cxd2),
		flux=cxd2(ll,cxd22);ener=ed2(cxd22);
		sigmavcx=sigcx(ener,md);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energd2=x;
    	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsd2=[tipsd2 -1/cp(1)]; nind2=[nind2 (-1/cp(1))^1.5*exp(cp(2))];
    	else
		tipsd2=[tipsd2 0]; nind2=[nind2 0]; 
    	end
    else
    	tipsd2=[tipsd2 0]; nind2=[nind2 0]; energd2=[];
    end

% ANR3

    if ~isempty(cxh3),
		flux=cxh3(ll,cxh33);ener=eh3(cxh33);
		sigmavcx=sigcx(ener,mh);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energh3=x;
     	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsh3=[tipsh3 -1/cp(1)]; ninh3=[ninh3 (-1/cp(1))^1.5*exp(cp(2))];
     	else
		tipsh3=[tipsh3 0]; ninh3=[ninh3 0]; 
     	end
    else
        tipsh3=[tipsh3 0]; ninh3=[ninh3 0]; energh3=[];
    end
	
    if ~isempty(cxd3),
		flux=cxd3(ll,cxd33);ener=ed3(cxd33);
		sigmavcx=sigcx(ener,md);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energd3=x;
    	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsd3=[tipsd3 -1/cp(1)]; nind3=[nind3 (-1/cp(1))^1.5*exp(cp(2))];
    	else
		tipsd3=[tipsd3 0]; nind3=[nind3 0]; 
    	end
    else
    	tipsd3=[tipsd3 0]; nind3=[nind3 0]; energd3=[];
    end

% ANR4

    if ~isempty(cxh4),
		flux=cxh4(ll,cxh44);ener=eh4(cxh44);
		sigmavcx=sigcx(ener,mh);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energh4=x;
     	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsh4=[tipsh4 -1/cp(1)]; ninh4=[ninh4 (-1/cp(1))^1.5*exp(cp(2))];
     	else
		tipsh4=[tipsh4 0]; ninh4=[ninh4 0];
     	end
    else
        tipsh4=[tipsh4 0]; ninh4=[ninh4 0]; energh4=[];
    end
	
    if ~isempty(cxd4),
		flux=cxd4(ll,cxd44);ener=ed4(cxd44);
		sigmavcx=sigcx(ener,md);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energd4=x;
    	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsd4=[tipsd4 -1/cp(1)]; nind4=[nind4 (-1/cp(1))^1.5*exp(cp(2))];
    	else
		tipsd4=[tipsd4 0]; nind4=[nind4 0];
    	end
    else
    	tipsd4=[tipsd4 0]; nind4=[nind4 0]; energd4=[];
    end

% ANR5

    if ~isempty(cxh5),
		flux=cxh5(ll,cxh55);ener=eh5(cxh55);
		sigmavcx=sigcx(ener,mh);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energh5=x;
     	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = (y-polyval(cp,x))./y*100;
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsh5=[tipsh5 -1/cp(1)]; ninh5=[ninh5 (-1/cp(1))^1.5*exp(cp(2))];
     	else
		tipsh5=[tipsh5 0]; ninh5=[ninh5 0];
     	end
    else
        tipsh5=[tipsh5 0]; ninh5=[ninh5 0]; energh5=[];
    end
	
    if ~isempty(cxd5),
		flux=cxd5(ll,cxd55);ener=ed5(cxd55);
		sigmavcx=sigcx(ener,md);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energd5=x;
    	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsd5=[tipsd5 -1/cp(1)]; nind5=[nind5 (-1/cp(1))^1.5*exp(cp(2))];
    	else
		tipsd5=[tipsd5 0]; nind5=[nind5 0];
    	end
    else
    	tipsd5=[tipsd5 0]; nind5=[nind5 0]; energd5=[];
    end

% ANR6

    if ~isempty(cxh6),
		flux=cxh6(ll,cxh66);ener=eh6(cxh66);
		sigmavcx=sigcx(ener,mh);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energh6=x;
     	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsh6=[tipsh6 -1/cp(1)]; ninh6=[ninh6 (-1/cp(1))^1.5*exp(cp(2))];
     	else
		tipsh6=[tipsh6 0]; ninh6=[ninh6 0];
     	end
    else
        tipsh6=[tipsh6 0]; ninh6=[ninh6 0]; energh6=[];
    end
	
    if ~isempty(cxd6),
		flux=cxd6(ll,cxd66);ener=ed6(cxd66);
		sigmavcx=sigcx(ener,md);
		kk=find(flux>100); x=ener(kk); y=log(flux(kk)./sqrt(ener(kk))./sigmavcx(kk)./atte);
		energd6=x;
    	if length(y)>3,
		[cp,p]=polyfit(x,y,1);
		xfit=0:0.1:max(x); yfit=cp(1)*xfit+cp(2);
		%
		%  rejet des points a 10 % du fit 
		%  V. Basiuk, 30 mars 2000
		%
		ecart = abs((y-polyval(cp,x))./y*100);
		ind   = ecart>10;
		if sum(~ind) > 1
		  [cp,p]=polyfit(x(~ind),y(~ind),1);
		end
		tipsd6=[tipsd6 -1/cp(1)]; nind6=[nind6 (-1/cp(1))^1.5*exp(cp(2))];
    	else
		tipsd6=[tipsd6 0]; nind6=[nind6 0];
    	end
    else
    	tipsd6=[tipsd6 0]; nind6=[nind6 0]; energd6=[];
    end


  end

ti_h=[tipsh1;tipsh2;tipsh3;tipsh4;tipsh5;tipsh6]';
ti_d=[tipsd1;tipsd2;tipsd3;tipsd4;tipsd5;tipsd6]';
ninh=[ninh1;ninh2;ninh3;ninh4;ninh5;ninh6]';
nind=[nind1;nind2;nind3;nind4;nind5;nind6]';
nhnd=ninh./nind;

% position ANRs

%zcx=ones(size(times))*[0.023 -0.16500 -0.315000 -0.48400 -0.61100 0];
%rcx=2.5*ones(size(zcx));
%[t,a,r0,z0,d0,piqd,e1,ep1,mode_sortie,certif]=tsgetgeo(choc,times,1);
%[rhocx,erreur]=rz2rho(rcx,zcx,a,r0,z0,d0,piqd,e1,ep1);


nhnd(nind==0|ninh==0)=zeros(size(nhnd(nind==0|ninh==0)));
ninh(imag(ninh)~=0)=zeros(size(ninh(imag(ninh)~=0)));
nind(imag(nind)~=0)=zeros(size(nind(imag(nind)~=0)));

