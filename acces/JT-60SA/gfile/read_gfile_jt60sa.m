function [gdata,ireadok] = read_gfile_funceast(filename, nfcoil, nesum, old_efit)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USAGE:  [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit)
%
%  PURPOSE: Read efit ASCII g-eqdsk file w/ flag iecurr=2 set
%
%  INPUTS: <default>
%     filename = gfile filename [Ex: '/u/leuer/efit/diiid/s87977/g087977.02620']
%     nfcoil = number of F-coils <18 for DIII-D> 
%     nesum =  number of E-coils <6 for DIII-D>
%     old_efit= For iecurr=2 several different flavors of efit output
%               exists. <Optional, default=0>  For old versions of efit 
%               (e.g. NSTX) use old_efit=1 since it uses older efit version.
%
%  OUTPUTS:
%    gdata = data structure containing EFIT G-file variables
%    ireadok = flag to report good read of gfile (0=bad, 1=good)
%
%  RESTRICTIONS:
%   gfile must be produced by EFIT run with flag iecurr=2 in order to 
%         write current density distribution on plasma grid & coil currents.
%        (NEW versions of EFIT use iplcout=2 instead of iecurr)
 
%  METHOD:  

%  WRITTEN BY:  Jim Leuer ON	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULTS:

if(nargin < 1)
    disp('ERROR read_gfile_func: Must specify name of g-file')
    gdata=[];
    ireadok=0;
    help read_gfile_func
    return    
 end
 if ~isstr(filename)
    disp([' %ERROR read_gfile: Filename must be a string',filename]);
    gdata=[];
    ireadok=0;
    return    
 end
 if(nargin < 2)
%    display('read_gfile_fun: Assuming DIII-D 18 PF Coils')
    nfcoil= 18;  
 end
 if(nargin<3)
%    display('read_gfile_fun: Assuming DIII-D 6 E-Coil segments')
    nesum=6;
 end
 if(nargin<4)
    old_efit=0;
 end

% -------------------------- Open File: filename
  filename = deblank(filename);
  fid= fopen(filename, 'r');

   if fid == (-1)
      disp([' %ERROR read_gfile: Couldn''t open file',filename]);
      gdata=[];
      ireadok= 0;
      return
   end

   ireadok= 1;

%Derive shotname from file:
shotname = sprintf('time_%s',filename);

% -------------------------- Start Read of File

  disp('First line:')
%      read (neqdsk,2000,err=30000) (case(i),i=1,6),imfit,nw,nh
  sline= fgetl(fid);
  disp(sline)
  ecase= sline(1,1:48);
  [dum,count]= sscanf(sline(1,49:length(sline)),'%d');
  imfit= dum(1);
  nw=    dum(2);
  nh=    dum(3);

%     read (neqdsk,2020) xdim,zdim,rzero,rgrid(1),zmid
  dum= fscanf(fid,'%f',5);
  xdim=     dum(1);
  zdim=     dum(2);
  rzero=    dum(3);
  rgrid1= dum(4);
  zmid=     dum(5);

% use this data to derive efit grid:
   rgefit= linspace(rgrid1, rgrid1+xdim, nw)';
   zgefit= linspace(zmid-zdim/2,zmid+zdim/2,nh)';
   drefit= (rgefit(end)-rgefit(1))/(nw-1);
   dzefit= (zgefit(end)-zgefit(1))/(nh-1);
   
%      read (neqdsk,2020) rmaxis,zmaxis,ssimag,ssibry,bzero
  dum= fscanf(fid,'%f',5);
  rmaxis= dum(1);
  zmaxis= dum(2);
  ssimag= dum(3);
  ssibry= dum(4);
  bzero=  dum(5);

%      write (neqdsk,2020) cpasma(jtime),ssimag,xdum,rmaxis,xdum
  dum= fscanf(fid,'%f',5);
  cpasma= dum(1);
  ssimag= dum(2);
  xdum=   dum(3);
  rmaxis= dum(4);
  xdum=   dum(5);

%      write (neqdsk,2020) zmaxis,xdum,ssibry,xdum,xdum
  dum= fscanf(fid,'%f',5);
  zmaxis= dum(1);
  xdum=   dum(2);
  ssibry= dum(3);
  xdum=   dum(4);
  xdum=   dum(5);

%      write (neqdsk,2020) (fpol(i),i=1,nw)
  fpol= fscanf(fid,'%f',nw);

%      write (neqdsk,2020) (pres(i),i=1,nw)
  pres= fscanf(fid,'%f',nw);

% ffprim and pprime are input with sign correct for EFIT
% eqdsk is written out as a negative if plasma current is positive

%      write (neqdsk,2020) (workk(i),i=1,nw)= -ffprim(i)
  ffprim= -sign(cpasma)*fscanf(fid,'%f',nw);

%      write (neqdsk,2020) (workk(i),i=1,nw)= -pprime
  pprime= -sign(cpasma)*fscanf(fid,'%f',nw);

%      write (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
  psirz= fscanf(fid,'%f',[nw,nh]);
%  note: nw= plasma width is stored in columns; nh= height in rows

%      write (neqdsk,2020) (qpsi(i),i=1,nw)
  qpsi= fscanf(fid,'%f',nw);

%      write (neqdsk,2022) nbbbs,limitr
  dum= fscanf(fid,'%f',2);
  nbbbs= dum(1);
  limitr= dum(2);

%      write (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
  dum= fscanf(fid,'%f',[2,nbbbs]);
  rbbbs= dum(1,:)';
  zbbbs= dum(2,:)';
  
%      write (neqdsk,2020) (xlim(i),ylim(i),i=1,limitr)
  dum= fscanf(fid,'%f',[2,limitr]);
  xlim= dum(1,:)';
  ylim= dum(2,:)';


% ------------------------end of original eqdsk read

%    fclose(fid);
%    clear dum sline xdum; 
%    return

%-------------------
% Create useful fluxes (even for iecurr not = 2)
% Convert flux objects to REAL units:
  psirz = -psirz'*2*pi;
  psimag = -ssimag*2*pi;
  psibry = -ssibry*2*pi;
  psibnd = -ssibry*2*pi;


% ------------------ 
% Continued read for read_eqdsk_iec2.m below:

%      read next line and display to see where we are
%  sline= fgetl(fid);
%  sline
%  sline= fgetl(fid);
%  sline

%      read (neqdsk,2024,err=30900,end=30900) kvtor,rvtor,nmass
  kvtor= 0;
  dum= fscanf(fid,'%f',3);

%___________________________________________
% 1st end of gfile with out currents
if (~isempty(dum)) 
    kvtor= dum(1);
    rvtor= dum(2);
    nmass= dum(3);

    if (kvtor > 0)
%      read (neqdsk,2020) (presw(i),i=1,mw)
       presw= fscanf(fid,'%f',nw);
%      read (neqdsk,2020) (preswp(i),i=1,mw)
       preswp= fscanf(fid,'%f',nw);
    end

% Read Ion Mass Density and associated Data (if present):
   if nmass>0
%    read (neqdsk,2020) (dmion(i),i=1,nw)
     dmion = fscanf(fid,'%f',nw);
   end

%% Loop to find ishot (if exists):
%   while 1
%     dum = fscanf(fid,'%f',1);
%     ***WORK HERE***
%   end
 
% Read New (sometime > 5/2/97) data in g-file:
  if old_efit==0
    rhovn = fscanf(fid,'%f',nw);    %always written now in g0-files!!
    keecur = fscanf(fid,'%f',1);
    if keecur>0
      workk=fscanf(fid,'%f',nw);
    end        
  end
   
% Read next line of (run-on) integers (and display to see where we are):
%  disp('Line after kvtor, etc...:')
    dum = fscanf(fid,'%f',4);              % The size should be 4 in fscanf(Fid,Format,Size). Added by Lzp(07/07/19)
    disp('ISHOT data line to check:')
    disp(dum)
%___________________________________________
%  2nd end of gfile without currents
  if (~isempty(dum))    
    disp('IECURR=2 data being read...')

%  read (neqdsk) rgrid(1),rgrid(nw),zgrid(1),zgrid(nh)
     dum = fscanf(fid,'%f',4);
  if 0 % now done above using 2nd and 3rd row read
     rgefit=linspace(dum(1),dum(2),nw)';
     zgefit=linspace(dum(3),dum(4),nh)';
     drefit = (dum(2)-dum(1))/(nw-1);
     dzefit = (dum(4)-dum(3))/(nh-1);
  end
%     disp('Grid extreme R,Z values are :')
%     disp(dum)

% read (neqdsk) (brsp(i),i=1,nfcoil)
    brsp = fscanf(fid,'%f',nfcoil);

% read (neqdsk) (ecurrt(i),i=1,nesum)
    ecurrt = fscanf(fid,'%f',nesum);

% read (neqdsk) (pcurrt(i),i=1,nwnh)
    [pcurrt,count] = fscanf(fid,'%f',nw*nh);   %current in grid elements (A)
    disp ('nh, nw, count');
    disp (nh)
    disp (nw)
    disp (count)
    if length(pcurrt) == nh*nw
      pcurrt = reshape(pcurrt,nh,nw);
    else
       disp('No IECURR=2 data.')   
    end
  else
    disp('No IECURR=2 data.')  
  end % if (~isempty(dum))

else
    disp('No IECURR=2 data.')  
end %if (~isempty(dum)) 
%----------------------------------- 
% END of gfile Read
  fclose(fid);
%----------------------------------- 
  if exist('pcurrt')~=1 

%   [pcurrt1a,jphi1a,pcursum1a]= plasma_current(-2*pi*psirz', pprime, ffprim,...
%        psimag,psibnd,rgefit,zgefit,rbbbs,zbbbs,nbbbs,cpasma,iplot,dofast);

    [pcurrt,jphi,pcursum]= plasma_current(psirz, pprime, ffprim,...
         psimag,psibnd,rgefit,zgefit,rbbbs,zbbbs,nbbbs,0,0);

    if(abs(pcursum-cpasma)/cpasma > 0.01)
       disp(['WARNING read_gfile_func: computed plasma current ' ...
		' differs from cpasma by > 1%'])
    end

    brsp = [];
    ecurrt = [];
  else
    pcurrt = [];
    pcursum = [];
    jphi = [];
     
  end
  
%   gdef = gfile_def; % definitions of all variables read in from G-file

  gdata = struct( ...	
    'brsp',brsp, ...
    'bzero',bzero, ...
    'cpasma',cpasma, ...
    'ecase',ecase, ...
    'ecurrt',ecurrt, ...
    'ffprim',ffprim, ...
    'fpol',fpol, ...
    'limitr',limitr, ...
    'nbbbs',nbbbs, ...
    'nh',nh, ...
    'nw',nw, ...
    'pcurrt',pcurrt, ...
    'pcursum',pcursum, ...
    'jphi',jphi, ...
    'pprime',pprime, ...
    'pres',pres, ...
    'psirz',psirz, ...
    'qpsi',qpsi, ...
    'rbbbs',rbbbs, ...
    'rgrid1',rgrid1, ...
    'rmaxis',rmaxis, ...
    'rzero',rzero, ...
    'ssibry',ssibry, ...
    'ssimag',ssimag, ...
    'xdim',xdim, ...
    'xlim',xlim, ...
    'ylim',ylim, ...
    'zbbbs',zbbbs, ...
    'zdim', zdim, ...
    'zmaxis', zmaxis, ...
    'rgefit', rgefit, ...
    'zgefit', zgefit, ...
    'drefit', drefit, ...
    'dzefit', dzefit, ...
    'psimag', psimag, ...
    'psibnd', psibnd, ...
    'psibry', psibry, ...
    'zmid', zmid,...
    'rhovn',rhovn);

  return

%  MODIFICATION HISTORY:
%     DAH  Added new variable definitions and unit converted quantities.
%     DAH  8/8/97  Modified to read new (~>5/2/97) g0-file format *only*
%     DAH  1/13/00 Modified to have nfcoil, nesum defaulted like read_gfile1
%		Now read_gfile *should* be usable for all non-D3D devices 
%		as well...
%     JAL 11-15-02 Fixed so it reads old efit iecurr=2 format for old_efit=1
%     MLW-JAL 1/11/2005 functionalized
%     jal 4/11/06 error found reading no plasma current gfiles
