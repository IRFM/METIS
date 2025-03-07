%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to plot relevant ECRH signals in a specific time window
% - without reading the signals first -
% The routine FCE has to be executed before this in order to read the signals
% t1 and t2a are the start and finish times for the plot window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sortie=StCont(t1,t2a,ecemax,stlevel,entree,data)

t      = entree.t;
PA1    = entree.PA1;
PA2    = entree.PA2;
Ptot   = entree.Ptot; 
tece   = entree.tece;
xece   = entree.xece;
rece   = entree.rece;
trece  = entree.trece; 
tIp    = entree.tIp;
Ip     = entree.Ip; 
tne    = entree.tne;
ne     = entree.ne;
tfci   = entree.tfci;
FCI    = entree.FCI; 
tlh    = entree.tlh;
PLH    = entree.PLH;
tvloop = entree.tvloop;
Vloop  = entree.Vloop;
tecenv = entree.tecenv; 
xecenv = entree.xecenv;
shot   = entree.shot;
sonde  = entree.sonde;





% numer of time points for the LH data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
l=length(tlh);

% if LH used - use non validated ECE data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if l>=100
	tece = tecenv ;
	xece = xecenv ;
end

% numer of time points for the ECE data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
lece=length(tece);

%Calculate Vector of Te variation from one point to the next and
% the assosicated time vector which has one point less than the ECE data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diffece = xece([2:lece],:)-xece([1:lece-1],:); 
tdiffece = tece([1:lece-1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialise Vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sttotpos = tdiffece-tdiffece;
sttot=sttotpos;
stperiodPreviousST =sttotpos;
timesinceST =sttotpos;
stperiod=sttotpos;
InvRadNormalOuter=sttotpos;
InvRadNormalInner=sttotpos;
InvRadInvertedOuter=sttotpos;
InvRadInvertedInner=sttotpos;
NormalCount=sttotpos;
InvertedOuterCount=sttotpos;
InvertedInnerCount=sttotpos; 
timelastst=sttotpos; 
TimeCrash=sttotpos;
nCrash=sttotpos;

%Initialise Scalars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tlst=0;
stPlST=0;
NSawtooth=1;
NST=1;

% MIN = minimium numer of timesteps between Sawteeth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MIN=5;


%End Initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Looop over all time points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=2:1:lece-1 

   % set time of last Sawtooth for this time point = time of last Sawtooth for previous time point
   % set Period of last Sawtooth for this time point = Period of last Sawtooth for previous time point
   % set Inversion Radii for this time point = Inversion Radii for previous time point
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   timelastst(n)=tlst;
   stperiodPreviousST(n)=stPlST;
   InvRadNormalOuter(n)=InvRadNormalOuter(n-1);
   InvRadNormalInner(n)=InvRadNormalInner(n-1);
   InvRadInvertedOuter(n)=InvRadInvertedOuter(n-1);
   InvRadInvertedInner(n)=InvRadInvertedInner(n-1);
   
   %initialise values for this time point
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   InvertedOuterCount(n)=0;
   InvertedInnerCount(n)=0;
   NormalFirst=23;
   NormalCount(n)=0;
   NormalLast=1;
   InvertedOuterFirst=23;
   InvertedOuterLast=23;
   InvertedInnerFirst=1;
   InvertedInnerLast=1;
   
   FalseNormalCount=0;
	
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Looop over all relevant channels
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for m=23:-1:1
   
 	%if n-nlast>=5000   
	%   nlast=n;
	%   nb=n
   	%end
	
	
	% if temperature drop is large enought to look like a non inverted sawtooth
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if diffece(n,m) <= -stlevel 
		
		% if at least one Outer Inverted sawtooth and no Inner inverted sawteeth have been found 
		% then this is a real sawtooth in which case - increase the count of channels that sees a normal sawtooth
		% and set the channel number for the last normal Sawtooth = the present channel
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if InvertedOuterCount(n)>= 1 & InvertedInnerCount(n)< 1
			
			% if this is the first Normal sawtooth
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if NormalCount(n)<1
				NormalFirst=m;
			end
			NormalCount(n)=NormalCount(n)+1;
			NormalLast=m;
		% else this a false sawtooth - in which case : Increase the count of false Sawteeth 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		else
		  	FalseNormalCount = FalseNormalCount+1;
		end		
	end
	
	% if temperature increase is large enought to look like an inverted sawtooth
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if diffece(n,m) >= stlevel
		
		% if no Normal sawteeth have been found 
		% then this is an rOuter inverted sawtooth in which case - increase the count of channels 
		% that sees an Outer inverted sawtooth
		% and set the channel number for the last Outer inverted Sawtooth = the present channel
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if NormalCount(n)<1
			
			% if this is the first Outer inverted sawtooth
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if InvertedOuterCount(n)<1
				InvertedOuterFirst=m;
			end
		
			InvertedOuterCount(n)=InvertedOuterCount(n)+1;
			InvertedOuterLast=m;
		
		% else (at least one Normal sawteeth have been found )
		% this is an Inner inverted sawtooth in which case - increase the count of channels 
		% that sees an Inner inverted sawtooth
		% and set the channel number for the last Inner inverted Sawtooth = the present channel
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		else
			% if this is the first Inner inverted sawtooth
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if InvertedInnerCount(n)<1
				InvertedInnerFirst=m;
			end
		
			InvertedInnerCount(n)=InvertedInnerCount(n)+1;
			InvertedInnerLast=m;
		end
	end
			
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %end loop over channels
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % if at least 2 channels sees a normal sawtooth
   % and at least 1 channel sees either an inner or an outer inverted sawtooth
   % and no channels sees this as a false sawtooth
   % and at least MIN time steps has passed since the last detected sawtooth
   % then this is really a sawtooth
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if NormalCount(n)>=2 & ...
      InvertedOuterCount(n)+InvertedInnerCount(n)>=1 & ...
      FalseNormalCount < 1 & ...
      n-NST>MIN
   		
	%set sttot(n) = 1 to indicate that a sawtooth i detected at this time point*
	%duration of this sawtooth = time now - time of last sawtooth
	%set time of last sawtooth = time now 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	sttot(n) = 1;
	stPlST=tdiffece(n)-tlst;
	tlst=tdiffece(n);
	
	%set n of last sawtooth = n now 
	%set nCrash(NSawtooth) which is the n assosicated with sawtooth number NSawtooth = n now 
	%set TimeCrash(NSawtooth) which is the time assosicated with sawtooth number NSawtooth = time now;	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	NST=n;		
	nCrash(NSawtooth)=n;
	TimeCrash(NSawtooth)=tlst;
	
	% increment the sawtooth counter
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	NSawtooth=NSawtooth+1;
	
	%find the sawtooth Inversion radii
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	InvRadNormalOuter(n)=rece(n,NormalFirst);
	InvRadNormalInner(n)=rece(n,NormalLast);
	InvRadInvertedOuter(n)=rece(n,InvertedOuterLast);
	InvRadInvertedInner(n)=rece(n,InvertedInnerFirst);	
   end
   
   
   % for all time points (values of n) 
   % set the the time since last sawtooth = time now - time of last sawtooth
   % set the instantaneous value of the sawtooth period = the maximun of ( the duration of the previous sawtooth and the time since last crash)
   timesinceST(n)=tdiffece(n)-timelastst(n);
   stperiod(n)=max(timesinceST(n),stperiodPreviousST(n));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end of loop over n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stuff for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ldiffece=length(tdiffece);
lsttotpos=length(sttotpos);
lsttot=length(sttot);
lstperiodPreviousST=length(stperiodPreviousST);
ltimesinceST=length(timesinceST);
lstperiod=length(stperiod);

t2=t2a+0.5*(t2a-t1);
numberplots=6;
plothight=0.8/numberplots;
if 1 < 2
h = findobj(0,'type','figure','tag','ddsinv');
if isempty(h)
       h=figure('tag','ddsinv');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',9,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1])




%h=axes('position',[0.1,0.1+2*plothight,0.8,plothight]);
h=axes('position',[0.1,0.1+4*plothight,0.8,2*plothight]);

if l>=10000
	plot(tecenv,xecenv(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 13 17 18 19 20 21 22 23 ]))
else
	plot(tece,xece(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 13 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31]))
end
%plot(tece,xece(:,[15 6 7 8 9  ]))
%plot(tece,xece,tdiffece,10*stcrit,'k-')
legend('ECE electron Temperatures')
%legend(' 5',' 6',' 7',' 8',' 9',' 10',' 11',' 12',' 13',' 14',' 15',' 13',' 17', ...
% 	'18',' 19',' 20',' 21',' 22',' 23')
axis([t1,t2,0,ecemax])
title( shot);
h=axes('position',[0.1,0.1+3*plothight,0.8,1*plothight]);
plot(tdiffece,-0.2*sttot,'k-',...
     tdiffece,stperiodPreviousST,'k-',tdiffece,timesinceST,'b-',tdiffece,stperiod,'r-')
legend('ST detected','Previous ST period', ...
       'Time Since Previous ST','ST period')
axis([t1,t2,-0.0,0.8])

h=axes('position',[0.1,0.1+1.4*plothight,0.8,1.6*plothight]);
plot(tdiffece,InvRadNormalOuter,'b-',tdiffece,InvRadNormalInner,'b-',...
     tdiffece,InvRadInvertedOuter,'g-',tdiffece,InvRadInvertedInner,'g-')
legend('Radius last normal ST','', ...
       'Radius first inverted ST','')
axis([t1,t2,2.0,2.8])

%h=axes('position',[0.1,0.1+1*plothight,0.8,1*plothight]);
%plot(tdiffece, NormalCount,'r-', tdiffece,-InvertedOuterCount-InvertedInnerCount,'b-')
%axis([t1,t2,-20.0,20.0])

h=axes('position',[0.1,0.1+0.7*plothight,0.8,0.7*plothight]);
if l>=100
	plot(tfci,FCI,'r-',tlh,PLH,'b-')
	legend('ICRH power','LH power')
else
	plot(tfci,FCI,'b-')
	legend('ICRH power')
end
axis([t1,t2,-0.05,9])

h=axes('position',[0.1,0.1+0*plothight,0.8,0.7*plothight]);
plot(t,(Ptot)','k-',t,(PA1)','r-',t,(PA2)','b-')
	legend('ECRH total Power','ECRH Power A1','ECRH Power A2')

axis([t1,t2,-50,700])

end



sortie.tdiffece            = tdiffece;
sortie.sttotpos            = sttotpos;
sortie.sttot               = sttot;
sortie.stperiodPreviousST  = stperiodPreviousST; 
sortie.timesinceST         = timesinceST;
sortie.stperiod            = stperiod;   
sortie.InvRadNormalOuter   = InvRadNormalOuter;
sortie.InvRadNormalInner   = InvRadNormalInner;
sortie.InvRadInvertedOuter = InvRadInvertedOuter; 
sortie.InvRadInvertedInner = InvRadInvertedInner;
sortie.NormalCount         = NormalCount;
sortie.InvertedOuterCount  = InvertedOuterCount;
sortie.InvertedInnerCount  = InvertedInnerCount;
sortie.TimeCrash           = TimeCrash;
sortie.nCrash              = nCrash;


%
% passage en rho
%
  val = zeros(length(data.gene.temps),4);
  long = length(data.gene.temps);
  pas = ceil(long/50);
  for k=1:pas:length(data.gene.temps)
    disp(['k =',int2str(k)])
    tana = data.gene.temps(k);
    ind  = min(find(tdiffece>tana));
    if ~isempty(ind)
  
      R = double(squeeze(data.equi.R(k,:,:)));
      Z = double(squeeze(data.equi.Z(k,:,:)));
      x = double(data.equi.rhoRZ(k,:))./data.equi.rhomax(k);
      if ~isnan(R)
        R(:,end) = [];
        Z(:,end) = [];
        R(1,:)   = [];
        Z(1,:)   = [];
        nx  = x' * ones(1,size(R,2));
        nx(1,:)  = [];
        val(k,1) = griddata(R,Z,nx,InvRadNormalOuter(ind),0);
        val(k,2) = griddata(R,Z,nx,InvRadNormalInner(ind),0);
        val(k,3) = griddata(R,Z,nx,InvRadInvertedOuter(ind),0);
        val(k,4) = griddata(R,Z,nx,InvRadInvertedInner(ind),0);
      end
    end
  end

sortie.val = val;
    
    