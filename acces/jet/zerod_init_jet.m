%============================================================
% acces au donnees de JET via MSD+
%============================================================
function z0dinput = zerod_init_jet(mode_exp,shot,gaz,temps,z0dinput)

%  try
%   if nargin('mdsdisconnect') > 0
%        mdsdisconnect( 'mdsplus.jet.efda.org');
%   else 
%        mdsdisconnect;
%   end
%  end

% default = no transp data read through SAL
transp_run_number = NaN;


langue                 =  'anglais';
% cas donnees JET
z0dinput.exp          = [];
% 1 - numero du choc
if isempty(shot) || isempty(gaz) || (~isappadata(0,'MDSPLUS_USERNAME') || ~isappadata(0,'SAL_AUTHORISATION_JET'))
	prompt={'shot number :','charge of main gas:','User name','SecurID key','TRANSP run or revision','Alternative TRANSP user'};
	def={'87737','1','prebut','1234567890','NaN',''};
	dlgTitle='Access to JET data';

    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(answer)
            z0dinput = [];
	    return
    end
    shot  = str2num(answer{1});
    gaz   = str2num(answer{2});
    user  = strtrim(answer{3});
    securidkey = strtrim(answer{4});
    transp_run_number =  str2num(answer{5});
    transp_user =  strtrim(answer{6});
    setappdata(0,'MDSPLUS_USERNAME',user);
    setappdata(0,'TRANSP_RUN',transp_run_number);
    setappdata(0,'TRANSP_USER',transp_user);
    % try SAL authetification if securidkey is defined
    accessok = false;
    if (length(answer{4}) >= 10) 
        try
            authorisation = sal_auth(user,securidkey);
            accessok = true;
        catch
            disp('SAL identification failed !');
        end
    end
    % try MDS+ reading
    liste =  {'ppf/@shot/MAGN/IPLA'};
    try
        data           = cgcgetjet(shot,liste,'','');
    catch
        data = [];
    end
    if ~isempty(data)
        disp('JET data access is available');
        accessok = true;
    else
        rmappdata(0,'MDSPLUS_USERNAME');
    end
    if ~accessok
        error('Access to JET data is no available');
    end
end
% lecture des donnees
%
% detection changement de mur: premier vrai choc avec ILW :80128 
if shot >= 80128
    flag_ilw = 1;
else
    flag_ilw = 0;
end

% preparation des donnees
%  
liste          = {};
liste{end+1}   = 'ppf/@shot/EFIT/BVAC';    % Bvide a 2.96 m 	
liste{end+1}   = 'ppf/@shot/MAGN/IPLA';    % plasma current
liste{end+1}   = 'ppf/@shot/MG2/VSU';      % tension par tour
liste{end+1}   = 'ppf/@shot/MAGN/VL';      % tension par tour
liste{end+1}   = 'ppf/@shot/MAGN/FLX';      % tension par tour
liste{end+1}   = 'ppf/@shot/MG3/TAU';      % taue (magnetique)
liste{end+1}   = 'ppf/@shot/MG2/FLX';      % flux au bord
%liste{end+1}   = 'ppf/@shot/MAGN/FLUX';      % flux au bord
liste{end+1}   = 'ppf/@shot/EFIT/QAX';     % estimation de q0
liste{end+1}   = 'ppf/@shot/EFIT/XLI';     % li
liste{end+1}   = 'ppf/@shot/EFIT/LI3M';     % li
liste{end+1}   = 'ppf/@shot/EFTM/LI3M';     % li
liste{end+1}   = 'ppf/@shot/KS3/ZEFH';     % <zeff> horizontal
liste{end+1}   = 'ppf/@shot/KS3/ZEFV';     % <zeff> vertical
liste{end+1}   = 'ppf/@shot/LIDX/TE0' ;    % Te0 
liste{end+1}   = 'ppf/@shot/LIDX/TEVL';    % <Te>
liste{end+1}   = 'ppf/@shot/LIDX/GAMT';    % piquage de Te
liste{end+1}   = 'ppf/@shot/LIDX/NE0';     % Ne0
liste{end+1}   = 'ppf/@shot/LIDX/NEVL';     % <Ne>
liste{end+1}   = 'ppf/@shot/LIDX/NELA';     % nbar
liste{end+1}   = 'ppf/@shot/LIDX/GAMN';    % piquage de NE
liste{end+1}   = 'ppf/@shot/KG1V/LID3';    % Line integrated (core) (m-2). Look out for fringe jumps (fringe jumps removed when status flag>0). In case of PPF not produced use the JPF signal DF/G1V-MV<3.
liste{end+1}   = 'ppf/@shot/KG1L/LAD3';    % Line averaged (core) m-3. It is created after the fringe jumps have been removed (upon request using ReqCo). It has the temporal resolution of EFIT signals
liste{end+1}   = 'ppf/@shot/HRTX/NELA';    % Line averaged m-3 (when HRTS/HRTX unavailable use LIDR/LIDX equivalents).
liste{end+1}   = 'ppf/@shot/HRTS/NE';      %Profile (m-3) mapped along the line of sight. Similar for LIDAR (LID/NE)
liste{end+1}   = 'ppf/@shot/LHCD/PTOT';    % puissance totale de LH
liste{end+1}   = 'ppf/@shot/MG3/WPD';      % WDIA
liste{end+1}   = 'ppf/@shot/EFIT/WDIA';      % WDIA
liste{end+1}   = 'ppf/@shot/EFIT/QWL';     % q au bord 
liste{end+1}   = 'ppf/@shot/EFIT/Q95';     % q 95
liste{end+1}   = 'ppf/@shot/EFIT/BTPD';    % betap
liste{end+1}   = 'ppf/@shot/EFIT/ELON';    % elongation du plasma
liste{end+1}   = 'ppf/@shot/NBI/PTOT';     % puissance totale de NBI
%liste{end+1}   = 'ppf/@shot/NBI/P080';     % puissance totale de NBI 80 keV
%liste{end+1}   = 'ppf/@shot/NBI/P140';     % puissance totale de NBI 140 keV
liste{end+1}   = 'ppf/@shot/ICRH/PTOT';    % puissance totale de ICRH
liste{end+1}   = 'ppf/@shot/ICRH/FRA';     % frequence ichr Hz antenne A
liste{end+1}   = 'ppf/@shot/ICRH/FRB';     % frequence ichr Hz antenne B
liste{end+1}   = 'ppf/@shot/ICRH/FRC';     % frequence ichr Hz antenne C
liste{end+1}   = 'ppf/@shot/ICRH/FRD';     % frequence ichr Hz antenne D
liste{end+1}   = 'ppf/@shot/ICRH/FRE';     % frequence ichr Hz antenne E (ILA)
liste{end+1}   = 'ppf/@shot/ICRH/PRFA';     % puissance couplee antenne A
liste{end+1}   = 'ppf/@shot/ICRH/PRFB';     % puissance couplee antenne B
liste{end+1}   = 'ppf/@shot/ICRH/PRFC';     % puissance couplee antenne C
liste{end+1}   = 'ppf/@shot/ICRH/PRFD';     % puissance couplee antenne D
liste{end+1}   = 'ppf/@shot/ICRH/PRFE';     % puissance couplee antenne E (ILA)
%  liste{end+1}   = 'ppf/@shot/RFE/FREQ';     % frequence ichr Hz antenne E (ILA)
%  liste{end+1}   = 'ppf/@shot/RFE/PRFE';     % puissance couplee antenne E (ILA)
liste{end+1}   = 'ppf/@shot/EFIT/TRIU';    % triangularite haute
liste{end+1}   = 'ppf/@shot/EFIT/TRIL';    % triangularite basse
liste{end+1}   = 'ppf/@shot/BOLO/TOPI';    % puissance rayonnee totale
liste{end+1}   = 'ppf/@shot/BOLO/TOPO';    % puissance rayonnee totale
liste{end+1}   = 'ppf/@shot/MG2/YOH';      % puissance ohmique
liste{end+1}   = 'ppf/@shot/MG3/YTO';      % pin
liste{end+1}   = 'ppf/@shot/NEX/WAL';      % densite au bord
liste{end+1}   = 'ppf/@shot/NET/WAL';      % densite au bord
liste{end+1}   = 'ppf/@shot/NION/NH';      % Total Number of H Ions  -    
liste{end+1}   = 'ppf/@shot/NION/ND';      % Total Number of D Ions   kf3 ?
liste{end+1}   = 'ppf/@shot/NION/NT';      %  Total Number of T Ions
liste{end+1}   = 'ppf/@shot/NION/NHE3';    %  Total Number of He%3 Ions
liste{end+1}   = 'ppf/@shot/NION/NHE4';    %  Total Number of He%4 Ions
liste{end+1}   = 'ppf/@shot/NION/NBE';     %  Total Number of Be Ions
liste{end+1}   = 'ppf/@shot/NION/NC';      %  Total Number of C Ions
liste{end+1}   = 'ppf/@shot/NION/NN';      %  Total Number of N Ions
liste{end+1}   = 'ppf/@shot/NION/NO';      %  Total Number of O Ions
liste{end+1}   = 'ppf/@shot/NION/NNE';     %  Total Number of Ne Ions
liste{end+1}   = 'ppf/@shot/NION/NNI';     %  Total Number of Ni Ions
liste{end+1}   = 'ppf/@shot/NION/NAR';     %  Total Number of AR Ions
liste{end+1}   = 'ppf/@shot/NION/NXE';     %  Total Number of XE Ions
liste{end+1}   = 'ppf/@shot/NION/NKR';     %  Total Number of KR Ions     
liste{end+1}   = 'ppf/@shot/NION/THE4';    %  temps de remplacement de HE4     
liste{end+1}   = 'ppf/@shot/TIN/RDTB';    %  taux de neutron DT  
liste{end+1}   = 'ppf/@shot/S3AD/AD35';    % D-alpha signal to LH transition
liste{end+1}   = 'ppf/@shot/CXSM/WTOT';    % energie thermique
liste{end+1}   = 'ppf/@shot/CXSM/ANGF';    % rotation
liste{end+1}   = 'ppf/@shot/CXFM/WTOT';    % energie thermique
liste{end+1}   = 'ppf/@shot/CXFM/ANGF';    % rotation
liste{end+1}   = 'ppf/@shot/CXFM/TIMX';    % Ti0
liste{end+1}   = 'ppf/@shot/CXFM/TI';      % profil de Ti
liste{end+1}   = 'ppf/@shot/KG1V/LID3';
liste{end+1}   = 'ppf/@shot/EFIT/FBND';    % flux au bord
liste{end+1}   = 'ppf/@shot/EFIT/RGEO';    % grand rayon
liste{end+1}   = 'ppf/@shot/EFIT/CR0';     % petit rayon
liste{end+1}   = 'ppf/@shot/EFIT/XIP';     % plasma current
liste{end+1}   = 'ppf/@shot/EFIT/AREA';    % surface section
liste{end+1}   = 'ppf/@shot/EFIT/VOLM';    % volume plasma
liste{end+1}   = 'ppf/@shot/EFIT/VJAC';    % vpr
liste{end+1}   = 'ppf/@shot/TION/TIAV';    % <TI>
liste{end+1} = 'ppf/@shot/NBI4/POW1';
liste{end+1} = 'ppf/@shot/NBI4/ENG1';
liste{end+1} = 'ppf/@shot/NBI4/AL1';
liste{end+1} = 'ppf/@shot/NBI4/PFR1';
liste{end+1} = 'ppf/@shot/NBI4/POW2';
liste{end+1} = 'ppf/@shot/NBI4/ENG2';
liste{end+1} = 'ppf/@shot/NBI4/AL2';
liste{end+1} = 'ppf/@shot/NBI4/PFR2';
liste{end+1} = 'ppf/@shot/NBI4/POW3';
liste{end+1} = 'ppf/@shot/NBI4/ENG3';
liste{end+1} = 'ppf/@shot/NBI4/AL3';
liste{end+1} = 'ppf/@shot/NBI4/PFR3';
liste{end+1} = 'ppf/@shot/NBI4/POW4';
liste{end+1} = 'ppf/@shot/NBI4/ENG4';
liste{end+1} = 'ppf/@shot/NBI4/AL4';
liste{end+1} = 'ppf/@shot/NBI4/PFR4';
liste{end+1} = 'ppf/@shot/NBI4/POW5';
liste{end+1} = 'ppf/@shot/NBI4/ENG5';
liste{end+1} = 'ppf/@shot/NBI4/AL5';
liste{end+1} = 'ppf/@shot/NBI4/PFR5';
liste{end+1} = 'ppf/@shot/NBI4/POW6';
liste{end+1} = 'ppf/@shot/NBI4/ENG6';
liste{end+1} = 'ppf/@shot/NBI4/AL6';
liste{end+1} = 'ppf/@shot/NBI4/PFR6';
liste{end+1} = 'ppf/@shot/NBI4/POW7';
liste{end+1} = 'ppf/@shot/NBI4/ENG7';
liste{end+1} = 'ppf/@shot/NBI4/AL7';
liste{end+1} = 'ppf/@shot/NBI4/PFR7';
liste{end+1} = 'ppf/@shot/NBI4/POW8';
liste{end+1} = 'ppf/@shot/NBI4/ENG8';
liste{end+1} = 'ppf/@shot/NBI4/AL8';
liste{end+1} = 'ppf/@shot/NBI4/PFR8';
liste{end+1} = 'ppf/@shot/NBI8/POW1';
liste{end+1} = 'ppf/@shot/NBI8/ENG1';
liste{end+1} = 'ppf/@shot/NBI8/AL1';
liste{end+1} = 'ppf/@shot/NBI8/PFR1';
liste{end+1} = 'ppf/@shot/NBI8/POW2';
liste{end+1} = 'ppf/@shot/NBI8/ENG2';
liste{end+1} = 'ppf/@shot/NBI8/AL2';
liste{end+1} = 'ppf/@shot/NBI8/PFR2';
liste{end+1} = 'ppf/@shot/NBI8/POW3';
liste{end+1} = 'ppf/@shot/NBI8/ENG3';
liste{end+1} = 'ppf/@shot/NBI8/AL3';
liste{end+1} = 'ppf/@shot/NBI8/PFR3';
liste{end+1} = 'ppf/@shot/NBI8/POW4';
liste{end+1} = 'ppf/@shot/NBI8/ENG4';
liste{end+1} = 'ppf/@shot/NBI8/AL4';
liste{end+1} = 'ppf/@shot/NBI8/PFR4';
liste{end+1} = 'ppf/@shot/NBI8/POW5';
liste{end+1} = 'ppf/@shot/NBI8/ENG5';
liste{end+1} = 'ppf/@shot/NBI8/AL5';
liste{end+1} = 'ppf/@shot/NBI8/PFR5';
liste{end+1} = 'ppf/@shot/NBI8/POW6';
liste{end+1} = 'ppf/@shot/NBI8/ENG6';
liste{end+1} = 'ppf/@shot/NBI8/AL6';
liste{end+1} = 'ppf/@shot/NBI8/PFR6';
liste{end+1} = 'ppf/@shot/NBI8/POW7';
liste{end+1} = 'ppf/@shot/NBI8/ENG7';
liste{end+1} = 'ppf/@shot/NBI8/AL7';
liste{end+1} = 'ppf/@shot/NBI8/PFR7';
liste{end+1} = 'ppf/@shot/NBI8/POW8';
liste{end+1} = 'ppf/@shot/NBI8/ENG8';
liste{end+1} = 'ppf/@shot/NBI8/AL8';
liste{end+1} = 'ppf/@shot/NBI8/PFR8';
%  	liste{end+1}   = 'ppf/@shot/NBI4/POW1';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/POW2';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/POW3';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/POW4';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/POW5';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/POW6';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/POW7';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/POW8';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/POW1';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/POW2';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/POW3';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/POW4';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/POW5';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/POW6';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/POW7';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/POW8';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/ENG1';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/ENG2';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/ENG3';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/ENG4';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/ENG5';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/ENG6';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/ENG7';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI4/ENG8';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/ENG1';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/ENG2';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/ENG3';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/ENG4';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/ENG5';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/ENG6';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/ENG7';    %  puissance NBI 1/16
%  	liste{end+1}   = 'ppf/@shot/NBI8/ENG8';    %  puissance NBI 1/16
%  	liste{end+1} = 'ppf/@shot/NBI4/AL1';
%  	liste{end+1} = 'ppf/@shot/NBI4/AL2';
%  	liste{end+1} = 'ppf/@shot/NBI4/AL3';
%  	liste{end+1} = 'ppf/@shot/NBI4/AL4';
%  	liste{end+1} = 'ppf/@shot/NBI4/AL5';
%  	liste{end+1} = 'ppf/@shot/NBI4/AL6';
%  	liste{end+1} = 'ppf/@shot/NBI4/AL7';
%  	liste{end+1} = 'ppf/@shot/NBI4/AL8';
%  	liste{end+1} = 'ppf/@shot/NBI8/AL1';
%  	liste{end+1} = 'ppf/@shot/NBI8/AL2';
%  	liste{end+1} = 'ppf/@shot/NBI8/AL3';
%  	liste{end+1} = 'ppf/@shot/NBI8/AL4';
%  	liste{end+1} = 'ppf/@shot/NBI8/AL5';
%  	liste{end+1} = 'ppf/@shot/NBI8/AL6';
%  	liste{end+1} = 'ppf/@shot/NBI8/AL7';
%  	liste{end+1} = 'ppf/@shot/NBI8/AL8';
liste{end+1}   = 'ppf/@shot/TIN/RNT';      %  neutron DD
liste{end+1}   = 'ppf/@shot/TIN/RDD';      %  neutron DD
liste{end+1}   = 'ppf/@shot/EFIT/RBND';    %  Rsepa
liste{end+1}   = 'ppf/@shot/EFIT/ZBND';    %  Zsepa
liste{end+1}   = 'ppf/@shot/EFIT/NBND';    %  NB sepa
liste{end+1}   = 'ppf/@shot/KK3P/RPED';    %   rayon du haut du piedestal
liste{end+1}   = 'ppf/@shot/KK3P/TPED';    %  temperature en haut du piedestal
%liste{end+1}   = 'ppf/@shot/NBIP/YBNS';    % neutron DD beam-plasma
liste{end+1}   = 'ppf/@shot/NBIP/YNS';    % neutron du au plasma (thermique)
% confinement
liste{end+1}   = 'ppf/@shot/SCAL/H98Y2';    % H_H
liste{end+1}   = 'ppf/@shot/SCAL/H98Y';    % H_H
liste{end+1}   = 'ppf/@shot/SCAL/H89';    % H_L
liste{end+1}   = 'ppf/@shot/SCAL/PT08';    % Martin
liste{end+1}   = 'ppf/@shot/SCAL/PLTH';    % PLOSS

% phase du LH
liste{end+1}   = 'jpf/@shot/LH/A1FPHO';
liste{end+1}   = 'jpf/@shot/LH/C1FPHO';
liste{end+1}   = 'jpf/@shot/LH/E1FPHO';
liste{end+1}   = 'jpf/@shot/LH/A1FPRO';
liste{end+1}   = 'jpf/@shot/LH/C1FPRO';
liste{end+1}   = 'jpf/@shot/LH/E1FPRO';

%
  % lecture globale   
%
data           = cgcgetjet(shot,liste,'','');
jpfdata        = data.jpf;
data           = data.ppf;
if ~isempty(data.EFIT.XIP.data)
	tip   = data.EFIT.XIP.t;
	ip    = data.EFIT.XIP.data;
	ip    = ip(tip >=0);
	temps = tip(tip >=0);
% il faut aussi trouver d'autres donnees pour la position du plasma
%  	elseif ~isempty(data.MAGN.IPLA.data)
%  		tip   = data.MAGN.IPLA.t;
%  		ip    = data.MAGN.IPLA.data;
%  		ip    = ip(tip >=0);
%  		temps = tip(tip >=0);
%  		if length(temps) > 1000
%  			temps = linspace(min(temps),max(temps),1000);
%  		end
else
	disp('No current plasma measurement available')
        z0dinput = [];
	return
end
signe_ip               = sign(mean(ip));
z0dinput.cons.temps    = temps;
z0dinput.cons.ip       = abs(ip);
if ~isempty(data.MAGN.FLX.data)
  z0dinput.cons.flux     = signe_ip .* zechan(data.MAGN.FLX.t,data.MAGN.FLX.data,temps,'nearest')./2./pi;    
elseif isempty(data.MG2.FLX.data) && ~isempty(data.MAGN.VL.data)
  flux = cumtrapz(data.MAGN.VL.t,data.MAGN.VL.data(:,3));
  z0dinput.cons.flux     = signe_ip .* zechan(data.MAGN.VL.t,flux,temps,'nearest')./2./pi; 
else
  z0dinput.cons.flux     = signe_ip .* zechan(data.MG2.FLX.t,data.MG2.FLX.data,temps,'nearest')./2./pi; 
end
if z0dinput.cons.flux(end) > z0dinput.cons.flux(1)
	z0dinput.cons.flux = - z0dinput.cons.flux;
end
z0dinput.geo.a         = interp10d(data.EFIT.CR0.t,data.EFIT.CR0.data,temps,'nearest');
z0dinput.geo.z0         = zeros(size(z0dinput.geo.a));
if isempty(data.EFIT.RGEO.data)
      Rgeo = 0.5*(min(data.EFIT.RBND.data,[],2) + max(data.EFIT.RBND.data,[],2));
      tR   = data.EFIT.RBND.t;
      z0dinput.geo.R         = interp10d(tR,Rgeo,temps,'nearest');
else
    z0dinput.geo.R         = interp10d(data.EFIT.RGEO.t,data.EFIT.RGEO.data,temps,'nearest');
end
z0dinput.geo.K         = interp10d(data.EFIT.ELON.t,data.EFIT.ELON.data,temps,'nearest');
z0dinput.geo.d         = max(interp10d(data.EFIT.TRIL.t,data.EFIT.TRIL.data,temps,'nearest') ,...
			  interp10d(data.EFIT.TRIU.t,data.EFIT.TRIU.data,temps,'nearest'));
rb0                    = abs(interp10d(data.EFIT.BVAC.t,data.EFIT.BVAC.data,temps,'nearest') .* 2.96);
z0dinput.geo.b0        = rb0 ./ z0dinput.geo.R;
z0dinput.geo.vp        = interp10d(data.EFIT.VOLM.t,data.EFIT.VOLM.data,temps,'nearest');
z0dinput.geo.sp        = interp10d(data.EFIT.AREA.t,data.EFIT.AREA.data,temps,'nearest');
if isempty(data.EFIT.VJAC.data)
	z0dinput.geo.sext      = NaN .* z0dinput.cons.temps;
else
	z0dinput.geo.sext      = interp10d(data.EFIT.VJAC.t,data.EFIT.VJAC.data(:,end),temps,'nearest');
end

if ~isempty(data.LIDX.NEVL.data)
	nemoy                  = interp10d(data.LIDX.NEVL.t,data.LIDX.NEVL.data,temps,'nearest');
	nemoy(~isfinite(nemoy))= min(nemoy(isfinite(nemoy)));
else
	nemoy = [];
end
if ~isempty(data.LIDX.NE0.data)
	ne0                    = interp10d(data.LIDX.NE0.t,data.LIDX.NE0.data,temps,'nearest');
	ne0(~isfinite(ne0))    = min(ne0(isfinite(ne0)));
else
	ne0 = [];	
end

if (length(nemoy) > 3) & (length(ne0) > 3)
	ane                    = ne0 ./ nemoy  - 1;
	nbar                   = ne0 ./  (2 .* gamma(ane + 1.5) ./ gamma(ane + 1 ) ./ sqrt(pi)); 	
	z0dinput.cons.nbar     = nbar;
else
	z0dinput.cons.nbar     = interp10d(data.KG1V.LID3.t,data.KG1V.LID3.data,temps,'nearest') ./ z0dinput.geo.a ./ 2 ./ z0dinput.geo.K;
	ane = ones(size(temps)); 
end
% utilisation du vrai nbar si disponible
if ~isempty(data.LIDX.NELA.data)
	nbar     = interp10d(data.LIDX.NELA.t,data.LIDX.NELA.data,temps,'nearest');
	indok    = find(isfinite(nbar));
	z0dinput.cons.nbar(indok) = nbar(indok);
end

% amelioration de la resolution temporelle sur nbar pour l'injection de glacon
if ~isempty(data.KG1V.LID3.data)
	nl      = interp10d(data.KG1V.LID3.t,sgolayfilt(data.KG1V.LID3.data,1,5),temps,'nearest');
	[pp,s,mu]= polyfit(nl,z0dinput.cons.nbar,1);
	xb =(nl - mu(1)) ./ mu(2);
	nbar = polyval(pp,xb);
	indok    = find(isfinite(nbar) & (nbar > 1e18));
	z0dinput.cons.nbar(indok) = nbar(indok);
end
%
if ~isempty(data.ICRH.PTOT.data)
    z0dinput.cons.picrh    = interp10d(data.ICRH.PTOT.t,data.ICRH.PTOT.data,temps,'nearest');
else
    z0dinput.cons.picrh    = zeros(size(temps));
end 
if ~isempty(data.LHCD.PTOT.data)
    z0dinput.cons.plh      = interp10d(data.LHCD.PTOT.t,data.LHCD.PTOT.data,temps,'nearest');
else
    z0dinput.cons.plh    = zeros(size(temps));
end 
if ~isempty(data.NBI.PTOT.data)
    z0dinput.cons.pnbi     = interp10d(data.NBI.PTOT.t,data.NBI.PTOT.data,temps,'nearest');
else
    z0dinput.cons.pnbi    = zeros(size(temps));
end 
z0dinput.cons.pecrh    = zeros(size(temps));
z0dinput.cons.hmore    = ones(size(temps)); 
if ~isempty(data.SCAL.H98Y2)
    data.SCAL.H98Y2 = data.SCAL.H98Y;
end
if ~isempty(data.SCAL.PT08.t) && isnumeric(data.SCAL.PT08.t) %%&& ~isempty(data.SCAL.H98Y2.t) && ~isempty(data.SCAL.H89.t)
    hh = interp10d(data.SCAL.H98Y2.t,sgolayfilt(data.SCAL.H98Y2.data,1,5),temps,'nearest');
    hl = interp10d(data.SCAL.H89.t,sgolayfilt(data.SCAL.H89.data,1,5),temps,'nearest');
    pmartin = interp10d(data.SCAL.PT08.t,sgolayfilt(data.SCAL.PT08.data,1,5),temps,'nearest');
    plcfs = interp10d(data.SCAL.PLTH.t,sgolayfilt(data.SCAL.PLTH.data,1,5),temps,'nearest');
    hmore = hh .* (plcfs > pmartin) + hl .* (plcfs <= pmartin);
    indok = find((hmore > 0.3) & (hmore < 2));
    if length(indok) > 2
        hmore = interp1(temps(indok),hmore(indok),temps,'nearest','extrap');
    end
    hmore(hmore < 0.3) = 1;
    hmore(hmore > 2)   = 1;
    %figure(21);plot(temps,hmore);
    z0dinput.cons.hmore    = hmore;
end
if  ~isempty(data.KS3.ZEFV.t)
	zeffm          = interp10d(data.KS3.ZEFV.t,data.KS3.ZEFV.data,temps,'nearest');
elseif ~isempty(data.KS3.ZEFH.data)
	zeffm          = interp10d(data.KS3.ZEFH.t,data.KS3.ZEFH.data,temps,'nearest');
else
	disp('No Zeff data')
	zeffm = 2 .* ones(size(temps));
end
if gaz == 2
	z0dinput.cons.zeff    = max(2,min(16,zeffm));
else
	z0dinput.cons.zeff    = max(1,min(7,zeffm));
end
	 
if ~isempty(data.EFIT.LI3M.data) & ~isstr(data.EFIT.LI3M.data)
    li                    = interp10d(data.EFIT.LI3M.t,data.EFIT.LI3M.data,temps,'nearest');
    librute               = li;
    li                    = medfilt1(li,3);
    z0dinput.option.li    = li(3);

elseif ~isempty(data.EFTM.LI3M.data) & ~isstr(data.EFTM.LI3M.data)
    li                    = interp10d(data.EFTM.LI3M.t,data.EFTM.LI3M.data,temps,'nearest');
    librute               = li;
    li                    = medfilt1(li,3);
    z0dinput.option.li    = li(3);

elseif ~isempty(data.EFIT.XLI.data) & ~isstr(data.EFIT.XLI.data)
    li                    = interp10d(data.EFIT.XLI.t,data.EFIT.XLI.data,temps,'nearest');
    librute               = li;
    li                    = medfilt1(li,3);
    z0dinput.option.li    = li(3);
else
    z0dinput.option.li    = 1;
    li                    = ones(size(temps));
    librute               = li;
end
        
% donnees experimentale pour comparaison
z0dinput.exp0d.temps  = temps;
z0dinput.exp0d.pin    = zechan(data.MG3.YTO.t,data.MG3.YTO.data,temps,'nearest');
z0dinput.exp0d.zeff   = z0dinput.cons.zeff;
z0dinput.exp0d.vp     = z0dinput.geo.vp;
z0dinput.exp0d.sp     = z0dinput.geo.sp;
z0dinput.exp0d.sext   = z0dinput.geo.sext;
z0dinput.exp0d.ane    = ane;
z0dinput.exp0d.nem    = nemoy;
z0dinput.exp0d.ne0    = ne0;
if ~isempty(data.NET.WAL.data)
	z0dinput.exp0d.nebord = zechan(data.NET.WAL.t,data.NET.WAL.data,temps,'nearest');
elseif  ~isempty(data.NEX.WAL.data)
	z0dinput.exp0d.nebord = zechan(data.NEX.WAL.t,data.NEX.WAL.data,temps,'nearest');
end
z0dinput.exp0d.nhem   = (real(zechan(data.NION.NHE3.t,data.NION.NHE3.data,temps,'nearest')) + ...
		      real(zechan(data.NION.NHE4.t,data.NION.NHE4.data,temps,'nearest'))) ./ z0dinput.exp0d.vp;
z0dinput.exp0d.nimpm  = (real(zechan(data.NION.NBE.t,data.NION.NBE.data,temps,'nearest')) + ...
			real(zechan(data.NION.NC.t,data.NION.NC.data,temps,'nearest')) + ...
			real(zechan(data.NION.NN.t,data.NION.NN.data,temps,'nearest')) + ...
			real(zechan(data.NION.NO.t,data.NION.NO.data,temps,'nearest')) + ...
			real(zechan(data.NION.NNI.t,data.NION.NNI.data,temps,'nearest')) + ...
			real(zechan(data.NION.NAR.t,data.NION.NAR.data,temps,'nearest')) + ...
			real(zechan(data.NION.NKR.t,data.NION.NKR.data,temps,'nearest')) + ...
			real(zechan(data.NION.NXE.t,data.NION.NXE.data,temps,'nearest')) + ...
			real(zechan(data.NION.NNE.t,data.NION.NNE.data,temps,'nearest'))) ./ z0dinput.exp0d.vp;
z0dinput.exp0d.nDm    = real(zechan(data.NION.ND.t,data.NION.ND.data,temps,'nearest')) ./ z0dinput.exp0d.vp;
z0dinput.exp0d.nTm    = real(zechan(data.NION.NT.t,data.NION.NT.data,temps,'nearest')) ./ z0dinput.exp0d.vp;
z0dinput.exp0d.n1m    = z0dinput.exp0d.nTm + z0dinput.exp0d.nDm +  ...
		      real(zechan(data.NION.NH.t,data.NION.NH.data,temps,'nearest')) ./ z0dinput.exp0d.vp;
z0dinput.exp0d.nim    = z0dinput.exp0d.n1m + z0dinput.exp0d.nimpm + z0dinput.exp0d.nhem;
z0dinput.exp0d.ni0    = z0dinput.exp0d.nim .* (1+ane);

z0dinput.exp0d.taue   = zechan(data.MG3.TAU.t,data.MG3.TAU.data,temps,'nearest'); 
z0dinput.exp0d.tauhe  = zechan(data.NION.THE4.t,data.NION.THE4.data,temps,'nearest');
if ~isempty(data.EFIT.WDIA.data)
	z0dinput.exp0d.w      = zechan(data.EFIT.WDIA.t,medfilt1(data.EFIT.WDIA.data,3),temps,'nearest'); 
else
	z0dinput.exp0d.w      = zechan(data.MG3.WPD.t,medfilt1(data.MG3.WPD.data,3),temps,'nearest'); 
end
z0dinput.exp0d.dwdt   = pdederive(temps,z0dinput.exp0d.w,2,2,1,1); 
if ~isempty(data.CXSM.WTOT.data)
    z0dinput.exp0d.wth    = zechan(data.CXSM.WTOT.t,data.CXSM.WTOT.data,temps,'nearest');
elseif ~isempty(data.CXFM.WTOT.data)
    z0dinput.exp0d.wth    = zechan(data.CXFM.WTOT.t,data.CXFM.WTOT.data,temps,'nearest');
else
    z0dinput.exp0d.wth    = z0dinput.exp0d.w; 
end
z0dinput.exp0d.dwthdt = pdederive(temps,z0dinput.exp0d.wth,2,2,1,1); 

z0dinput.exp0d.te0    = zechan(data.LIDX.TE0.t,data.LIDX.TE0.data,temps,'nearest');
z0dinput.exp0d.tem    = zechan(data.LIDX.TEVL.t,data.LIDX.TEVL.data,temps,'nearest');
z0dinput.exp0d.tite   = zechan(data.TION.TIAV.t,data.TION.TIAV.data,temps,'nearest') ./ z0dinput.exp0d.tem;
z0dinput.exp0d.ate    = z0dinput.exp0d.te0 ./ z0dinput.exp0d.tem -1;
z0dinput.exp0d.ape    = ane + z0dinput.exp0d.ate; 
  
z0dinput.exp0d.pohm   = zechan(data.MG2.YOH.t,data.MG2.YOH.data,temps,'nearest'); 
if ~isempty(data.MAGN.VL.data)
  z0dinput.exp0d.vloop  = signe_ip .* zechan(data.MAGN.VL.t,data.MAGN.VL.data(:,3),temps,'nearest'); 
else
  z0dinput.exp0d.vloop  = signe_ip .* zechan(data.MG2.VSU.t,data.MG2.VSU.data,temps,'nearest'); 
end
z0dinput.exp0d.qa     = zechan(data.EFIT.QWL.t,data.EFIT.QWL.data,temps,'nearest'); 
z0dinput.exp0d.q95    = zechan(data.EFIT.Q95.t,data.EFIT.Q95.data,temps,'nearest'); 
z0dinput.exp0d.q0     = zechan(data.EFIT.QAX.t,data.EFIT.QAX.data,temps,'nearest'); 
z0dinput.exp0d.betap  = zechan(data.EFIT.BTPD.t,data.EFIT.BTPD.data,temps,'nearest'); 
z0dinput.exp0d.ip     = z0dinput.cons.ip;
z0dinput.exp0d.pfus   = 3.56e6 .* 1.602176462e-19 .* zechan(data.TIN.RDTB.t,data.TIN.RDTB.data,temps,'nearest');
if isempty(data.BOLO.TOPI.data)
  z0dinput.exp0d.prad   = zechan(data.BOLO.TOPO.t,data.BOLO.TOPO.data,temps,'nearest'); 
else
  z0dinput.exp0d.prad   = zechan(data.BOLO.TOPI.t,data.BOLO.TOPI.data,temps,'nearest'); 
end
z0dinput.exp0d.ploss  = z0dinput.exp0d.pin - z0dinput.exp0d.prad ./ 3;

z0dinput.exp0d.nbar  = z0dinput.cons.nbar;
z0dinput.exp0d.li    = librute;
z0dinput.exp0d.picrh = z0dinput.cons.picrh;
z0dinput.exp0d.pnbi  = z0dinput.cons.pnbi;
z0dinput.exp0d.plh   = z0dinput.cons.plh;
z0dinput.exp0d.pecrh = z0dinput.cons.pecrh;
z0dinput.exp0d.edgeflux  = z0dinput.cons.flux;
try 
	modeha = sgolayfilt(data.S3AD.AD35.data,1,51);
catch
	modeha = medfilt1(data.S3AD.AD35.data,1,3);
end

modeha = modeha - std(modeha);
modeha(modeha <0) =0;
%modeha = zechan(data.S3AD.AD35.t,modeha,temps,'nearest') ./ z0dinput.exp0d.nem;
if isempty( z0dinput.exp0d.nem)
    modeha = zechan(data.S3AD.AD35.t,modeha,temps,'nearest') ./ z0dinput.cons.nbar;
else
    modeha = zechan(data.S3AD.AD35.t,modeha,temps,'nearest') ./ z0dinput.exp0d.nem;
end
modeha = modeha - min(modeha);
modeha = modeha ./ max(modeha);
z0dinput.exp0d.modeh =modeha; 

z0dinput.exp0d.ndd           = zechan(data.TIN.RNT.t,data.TIN.RNT.data,temps,'nearest');
%z0dinput.exp0d.ndd_th        = zechan(data.NBIP.YPNS.t,data.NBIP.YPNS.data,temps,'nearest');
%z0dinput.exp0d.ndd_nbi_th    = zechan(data.NBIP.YBNS.t,data.NBIP.YBNS.data,temps,'nearest');


% rotation 
if ~isempty(data.CXSM.ANGF.data)
      z0dinput.exp0d.wrad =  zechan(data.CXSM.ANGF.t,max(data.CXSM.ANGF.data,[],2),temps,'nearest');
elseif ~isempty(data.CXFM.ANGF.data)
    z0dinput.exp0d.wrad =  zechan(data.CXFM.ANGF.t,max(data.CXFM.ANGF.data,[],2),temps,'nearest');
else
    z0dinput.exp0d.wrad = NaN .* ones(size(temps));
end
switch flag_ilw
case 1
    z0dinput.option.zmax       = 6;
    z0dinput.option.zimp       = 4;
    z0dinput.option.rimp       = 0.1;
    z0dinput.option.zeff       = 0;
otherwise
    z0dinput.option.zmax       = 4;
    z0dinput.option.zimp       = 6;
    z0dinput.option.rimp       = 0.1;
    z0dinput.option.zeff       = 7;
end
z0dinput.option.nphi      = 15;
z0dinput.option.angle_nbi = 90;
z0dinput.option.angle_nbi2 = 90;
z0dinput.option.einj      = 100e3;
z0dinput.option.modeh     = 1;
iso                       = real(z0dinput.exp0d.nTm ./max(eps,z0dinput.exp0d.nDm));
ind                       = find(isfinite(iso) & iso > 0);
if length(ind) <= 5
	z0dinput.cons.iso       =  zeros(size(temps));
	if gaz  == 1
		z0dinput.option.gaz       =  2;
	else
		z0dinput.option.gaz       =  4;
	end
else
	indnk                   = find(~isfinite(iso) | iso < 0);
	iso(indnk)              = mean(iso(ind));
	z0dinput.option.gaz       =  3;
	fprintf('INFORMATION: Tritium detected !\n');
	z0dinput.cons.iso       = iso;
end

%calcul de l'energie des neutres	
nbi4 = data.NBI4;
nbi8 = data.NBI8;

% boucle sur les pinis
%  po = 0;
%  eo = 0;
%  ro = 0;
%  zo = 0;
%  zs = 0.15;
%  zu = 0.05;
%  pnbi_s = zeros(size(z0dinput.cons.pnbi));
%  for kk = 1:8
%      el = getfield(nbi4,sprintf('ENG%d',kk),'data');
%      pl = getfield(nbi4,sprintf('POW%d',kk),'data');
%      pfr = getfield(nbi4,sprintf('PFR%d',kk),'data');
%      al = getfield(nbi4,sprintf('AL%d',kk),'data');
%      tl = getfield(nbi4,sprintf('POW%d',kk),'t');
%      if ~isempty(el) & ~isempty(pl)
%  	po = po + sum(pl);
%  	if ~isempty(pfr)
%  		%eo = eo + sum(pl .* el .* (pfr(1) + pfr(2) ./ 2  + pfr(3) ./ 3) ./ max(1,sum(pfr)));
%  		eo = eo + sum(pl .* el .*  sum(pfr) ./ max(1,pfr(1) + pfr(2) .* 2  + pfr(3) .* 3));
%  	else
%  		eo = eo + sum(pl .* el);
%  	end
%  	switch kk
%  	case {1,2,7,8}
%  		ro = ro + sum(pl .* 1.85);
%  	otherwise
%  		ro = ro + sum(pl .* 1.3);
%  	end
%  	if isempty(al)
%  		zo = zo + sum(pl .* zs);
%  	elseif al == 1
%  		zo = zo + sum(pl .* zu);
%  	else
%  		zo = zo + sum(pl .* zs);
%  	end
%  	p1 = interp10d(tl,pl,temps,'nearest');
%  	pnbi_s = pnbi_s + p1;
%      end
%      el = getfield(nbi8,sprintf('ENG%d',kk),'data');
%      pl = getfield(nbi8,sprintf('POW%d',kk),'data');
%      pfr = getfield(nbi8,sprintf('PFR%d',kk),'data');
%      al = getfield(nbi8,sprintf('AL%d',kk),'data');
%      tl = getfield(nbi8,sprintf('POW%d',kk),'t');
%      if ~isempty(el) & ~isempty(pl)
%  	po = po + sum(pl);
%  	if ~isempty(pfr)
%  		%eo = eo + sum(pl .* el .* (pfr(1) + pfr(2) ./ 2  + pfr(3) ./ 3) ./ max(1,sum(pfr)));
%  		eo = eo + sum(pl .* el .*  sum(pfr) ./ max(1,pfr(1) + pfr(2) .* 2  + pfr(3) .* 3));
%  	else
%  		eo = eo + sum(pl .* el);
%  	end
%  	switch kk
%  	case {1,2,7,8}
%  		ro = ro + sum(pl .* 1.85);
%  	otherwise
%  		ro = ro + sum(pl .* 1.3);
%  	end
%  	if isempty(al)
%  		zo = zo + sum(pl .* zs);
%  	elseif al == 1
%  		zo = zo + sum(pl .* zu);
%  	else
%  		zo = zo + sum(pl .* zs);
%  	end
%  	p1 = interp10d(tl,pl,temps,'nearest');
%  	pnbi_s = pnbi_s + p1;
%      end
%  end
%  if  (eo >0) & (po >0)
%  	z0dinput.option.einj = eo ./ po;
%  	z0dinput.option.rtang      = ro ./ po;		
%  	z0dinput.option.zext       = zo ./ po;		
%  else
%  	z0dinput.option.einj = 8e4;
%  	z0dinput.option.rtang      = (1.85 + 1.31) ./ 2;
%  	z0dinput.option.zext       =  0.1;		
%  end
%         if any(pnbi_s)
%  		    z0dinput.cons.pnbi = pnbi_s;
%         end
% boucle sur les pinis
po_1 = 0;
eo_1 = 0;
ro_1 = 0;
zo_1 = 0;
zs_1 = 0.15;
zu_1 = 0.05;
po_2 = 0;
eo_2 = 0;
ro_2 = 0;
zo_2 = 0;
zs_2 = 0.15;
zu_2 = 0.05;
pnbi_s1 = zeros(size(z0dinput.cons.pnbi));
pnbi_s2 = zeros(size(z0dinput.cons.pnbi));
for kk = 1:8
    el = getfield(nbi4,sprintf('ENG%d',kk),'data');
    pl = getfield(nbi4,sprintf('POW%d',kk),'data');
    pfr = getfield(nbi4,sprintf('PFR%d',kk),'data');
    al = getfield(nbi4,sprintf('AL%d',kk),'data');
    tl = getfield(nbi4,sprintf('POW%d',kk),'t');
    if ~isempty(el) & ~isempty(pl)
	switch kk
	case {1,2,7,8}
		po_1 = po_1 + sum(pl);
		ro_1 = ro_1 + sum(pl .* 1.85);
		if ~isempty(pfr)
			eo_1 = eo_1 + sum(pl .* el .*  sum(pfr) ./ max(1,pfr(1) + pfr(2) .* 2  + pfr(3) .* 3));
		else
			eo_1 = eo_1 + sum(pl .* el);
		end
		if isempty(al)
			zo_1 = zo_1 + sum(pl .* zs_1);
		elseif al == 1
			zo_1 = zo_1 + sum(pl .* zu_1);
		else
			zo_1 = zo_1 + sum(pl .* zs_1);
		end
		p1 = interp10d(tl,pl,temps,'nearest');
	        pnbi_s1 = pnbi_s1 + p1;
	otherwise
		po_2 = po_2 + sum(pl);
		ro_2 = ro_2 + sum(pl .* 1.3);
		if ~isempty(pfr)
			eo_2 = eo_2 + sum(pl .* el .*  sum(pfr) ./ max(1,pfr(1) + pfr(2) .* 2  + pfr(3) .* 3));
		else
			eo_2 = eo_2 + sum(pl .* el);
		end
		if isempty(al)
			zo_2 = zo_2 + sum(pl .* zs_2);
		elseif al == 1
			zo_2 = zo_2 + sum(pl .* zu_2);
		else
			zo_2 = zo_2 + sum(pl .* zs_2);
		end
		p2 = interp10d(tl,pl,temps,'nearest');
	        pnbi_s2 = pnbi_s2 + p2;
	end
    end
    el = getfield(nbi8,sprintf('ENG%d',kk),'data');
    pl = getfield(nbi8,sprintf('POW%d',kk),'data');
    pfr = getfield(nbi8,sprintf('PFR%d',kk),'data');
    al = getfield(nbi8,sprintf('AL%d',kk),'data');
    tl = getfield(nbi8,sprintf('POW%d',kk),'t');
    if ~isempty(el) & ~isempty(pl)
	switch kk
	case {1,2,7,8}
		po_1 = po_1 + sum(pl);
		ro_1 = ro_1 + sum(pl .* 1.85);
		if ~isempty(pfr)
			eo_1 = eo_1 + sum(pl .* el .*  sum(pfr) ./ max(1,pfr(1) + pfr(2) .* 2  + pfr(3) .* 3));
		else
			eo_1 = eo_1 + sum(pl .* el);
		end
		if isempty(al)
			zo_1 = zo_1 + sum(pl .* zs_1);
		elseif al == 1
			zo_1 = zo_1 + sum(pl .* zu_1);
		else
			zo_1 = zo_1 + sum(pl .* zs_1);
		end
		p1 = interp10d(tl,pl,temps,'nearest');
	        pnbi_s1 = pnbi_s1 + p1;
	otherwise
		po_2 = po_2 + sum(pl);
		ro_2 = ro_2 + sum(pl .* 1.3);
		if ~isempty(pfr)
			eo_2 = eo_2 + sum(pl .* el .*  sum(pfr) ./ max(1,pfr(1) + pfr(2) .* 2  + pfr(3) .* 3));
		else
			eo_2 = eo_2 + sum(pl .* el);
		end
		if isempty(al)
			zo_2 = zo_2 + sum(pl .* zs_2);
		elseif al == 1
			zo_2 = zo_2 + sum(pl .* zu_2);
		else
			zo_2 = zo_2 + sum(pl .* zs_2);
		end
		p2 = interp10d(tl,pl,temps,'nearest');
	        pnbi_s2 = pnbi_s2 + p2;
	end
    end
end

if  (eo_1 >0) & (po_1 > 0)
	z0dinput.option.einj       = eo_1 ./ po_1;
	z0dinput.option.rtang      = ro_1 ./ po_1;		
	z0dinput.option.zext       = zo_1 ./ po_1;
else
	z0dinput.option.einj       = 8e4;
	z0dinput.option.rtang      = 1.85;
	z0dinput.option.zext       = 0.1;		
end

if  (eo_2 >0) & (po_2 >0)
    if  (po_1  == 0)
	z0dinput.option.einj       = eo_2 ./ po_2;
	z0dinput.option.rtang      = ro_2 ./ po_2;		
	z0dinput.option.zext       = zo_2 ./ po_2;	
 	z0dinput.option.nb_nbi     = 1;
 	pnbi_s1                    = pnbi_s2;
 	pnbi_s2(:)                 = 0;  
	z0dinput.option.einj2      = 8e4;
	z0dinput.option.rtang2     = 1.31;
	z0dinput.option.zext2      = 0.1;		
    else
	z0dinput.option.nb_nbi      = 2;
	z0dinput.option.einj2       = eo_2 ./ po_2;
	z0dinput.option.rtang2      = ro_2 ./ po_2;		
	z0dinput.option.zext2       = zo_2 ./ po_2;	
    end
else
	z0dinput.option.einj2       = 8e4;
	z0dinput.option.rtang2      = 1.31;
	z0dinput.option.zext2       = 0.1;		
end
pnbi_s = pnbi_s1 + sqrt(-1) .* pnbi_s2;
z0dinput.cons.pnbi     = pnbi_s .* z0dinput.cons.pnbi ./ max(1,real(pnbi_s) + imag(pnbi_s));


z0dinput.option.etalh     = 0.8;
z0dinput.option.lhmode    = 3;
z0dinput.option.xlh        = 0.4;
z0dinput.option.dlh        = 0.3;
z0dinput.option.wlh        = 11e-3 * 32;
z0dinput.option.npar0      = 2;
% z0dinput.option.zext       = 0.1;
z0dinput.option.ane       = 0;
z0dinput.option.vane       = mean(1 + ane);
% il faut affiner ce parametre
z0dinput.machine     = 'JET';
z0dinput.shot        = shot;
z0dinput.option.configuration = 2;

% activation des dents de scies	
z0dinput.option.qdds         = 1;
z0dinput.option.kidds        = 3;
z0dinput.option.xiioxie      = 0;
z0dinput.option.xieorkie     = 0;
	
	
% calcul de la frequence ICRH la plus probable (pour H)
% calcul detailler de la frequence
ppi = 0;
ffi = 0;
if ~isempty(data.ICRH.PRFA.data) & ~isempty(data.ICRH.FRA.data)
      ppi = sum(data.ICRH.PRFA.data) + ppi;
      ffi = sum(data.ICRH.FRA.data .* data.ICRH.PRFA.data) + ffi;
end
if ~isempty(data.ICRH.PRFB.data) & ~isempty(data.ICRH.FRB.data)
      ppi = sum(data.ICRH.PRFB.data) + ppi;
      ffi = sum(data.ICRH.FRB.data .* data.ICRH.PRFB.data) + ffi;
end
if ~isempty(data.ICRH.PRFC.data) & ~isempty(data.ICRH.FRC.data)
      ppi = sum(data.ICRH.PRFC.data) + ppi;
      ffi = sum(data.ICRH.FRC.data .* data.ICRH.PRFC.data) + ffi;
end
if ~isempty(data.ICRH.PRFD.data) & ~isempty(data.ICRH.FRD.data)
      ppi = sum(data.ICRH.PRFD.data) + ppi;
      ffi = sum(data.ICRH.FRD.data .* data.ICRH.PRFD.data) + ffi;
end
if ~isempty(data.ICRH.PRFE.data) & ~isempty(data.ICRH.FRE.data)
      ppi = sum(data.ICRH.PRFE.data) + ppi;
      ffi = sum(data.ICRH.FRE.data .* data.ICRH.PRFE.data) + ffi;
%  elseif ~isempty(data.RFE.PRFE.data) & ~isempty(data.RFE.FREQ.data)
%        ppi = sum(data.RFE.PRFE.data) + ppi;
%        ffi = sum(data.RFE.FREQ.data .* data.RFE.PRFE.data) + ffi;
end

if ppi > 0
      z0dinput.option.freq = ffi ./ ppi ./ 1e6;
else
      bres   = z0dinput.geo.R .*  z0dinput.geo.b0 ./ (z0dinput.geo.R - 0.15 .* z0dinput.geo.a);
      %bres   =  2 .* pi .* option.freq .* 1e6 ./ (95.5e6 .*zg ./ ag);
      freq = bres .* 95.5e6 ./ 1e6 ./ 2 ./ pi;
      z0dinput.option.freq = trapz(z0dinput.cons.temps,z0dinput.cons.picrh .*  freq) ./ ...
			      max(1,trapz(z0dinput.cons.temps,z0dinput.cons.picrh));
      if z0dinput.option.freq <= 0
	      z0dinput.option.freq = max(1,mean(freq));
      end
end
z0dinput.option.signe = sign(mean(data.EFIT.XIP.data .* data.EFIT.BVAC.data));	


% Phase et N// pour LH
if ~isempty(data.LHCD.PTOT.data)
    philh = zeros(size(data.LHCD.PTOT.data));
    prlh  = zeros(size(data.LHCD.PTOT.data));
    if ~isempty(jpfdata.LH.E1FPHO.data) && ~isempty(jpfdata.LH.E1FPRO.data)
        ploc  =  interp1(jpfdata.LH.E1FPRO.t,jpfdata.LH.E1FPRO.data,data.LHCD.PTOT.t,'nearest');
        philh =  philh + interp1(jpfdata.LH.E1FPHO.t,jpfdata.LH.E1FPHO.data,data.LHCD.PTOT.t,'nearest')  .* ploc;
        prlh  =  prlh + ploc;
    end
    if ~isempty(jpfdata.LH.C1FPHO.data) && ~isempty(jpfdata.LH.C1FPRO.data)
        ploc  =  interp1(jpfdata.LH.C1FPRO.t,jpfdata.LH.C1FPRO.data,data.LHCD.PTOT.t,'nearest');
        philh =  philh + interp1(jpfdata.LH.C1FPHO.t,jpfdata.LH.C1FPHO.data,data.LHCD.PTOT.t,'nearest')  .* ploc;
        prlh  =  prlh + ploc;
    end
    if ~isempty(jpfdata.LH.A1FPHO.data) && ~isempty(jpfdata.LH.A1FPRO.data)
        ploc  =  interp1(jpfdata.LH.A1FPRO.t,jpfdata.LH.A1FPRO.data,data.LHCD.PTOT.t,'nearest');
        philh =  philh + interp1(jpfdata.LH.A1FPHO.t,jpfdata.LH.A1FPHO.data,data.LHCD.PTOT.t,'nearest')  .* ploc;
        prlh  =  prlh + ploc;
    end
    indok = find(data.LHCD.PTOT.data >1e5);
    if length(indok)>=3
        philhav = trapz(data.LHCD.PTOT.t(indok),philh(indok)) ./ max(1,trapz(data.LHCD.PTOT.t(indok),prlh(indok)));
        phi_in  = [-90 -45 0 45 90];
        npar_in = [1.4 1.6 1.8 2.1 2.3];
        nparlh  = interp1(phi_in,npar_in,philhav,'pchip','extrap');
        Dn    = 2.01 - 0.63 .*  (abs(nparlh - 1.8) + 1.8);
        if isfinite(nparlh)
            z0dinput.option.etalh     = Dn;
            z0dinput.option.npar0     = nparlh;
        end
    end
end
  
	
% utilisation de la separatrice si possible
if ~isempty(data.EFIT.RBND.data) && ~isempty(data.EFIT.ZBND.data) && ~isempty(data.EFIT.NBND.data)
	
	vt             = ones(size(z0dinput.cons.temps));
	teta           = linspace(0,2*pi,201);
	ve             = ones(size(teta));
	z0dinput.exp0d.Rsepa = NaN .* vt * ve;
	z0dinput.exp0d.Zsepa = NaN .* vt * ve;		
	for k = 1:length(vt)
		indn   = find(data.EFIT.NBND.t >= z0dinput.cons.temps(k),1);
		if isempty(indn) 
			indn  = min(length(data.EFIT.NBND.t),length(data.EFIT.RBND.t));
		end
		rr     = data.EFIT.RBND.data(indn,1:data.EFIT.NBND.data(indn));
		zz     = data.EFIT.ZBND.data(indn,1:data.EFIT.NBND.data(indn));
		rr(:,end+1)    = rr(:,1);
		zz(:,end+1)    = zz(:,1);
		rmin  = min(rr);
		rmax  = max(rr);
		ra    = max(1,0.5 .* (rmin + rmax));
		a     = max(0.01,0.5 .* (rmax - rmin));
		% SECURITE TEMPS ETRANGE
		ra = max(a .* 1.01 , ra);
		zmin  = min(zz);
		zmax  = max(zz);
		za    = (zmin + zmax) ./ 2;
		b     = 0.5 .* (zmax - zmin);
		elon     = max(0.5,b ./ a);
		mask1 = (zz == zmax);
		mask2 = (zz == zmin);	
		rzmax = max(rr .* mask1,[],2);
		rzmin = max(rr .* mask2,[],2);
		cl    = ra - rzmin;
		cu    = ra - rzmax;
		d     = (cl+cu) ./2 ./ a;
		z0dinput.geo.a(k)      = a;
		z0dinput.geo.R(k)      = ra;
		z0dinput.geo.z0(k)     = za;
		z0dinput.geo.K(k)      = elon;
		z0dinput.geo.d(k)      = d;	
					
		% iso angle
		cl   = (rr - ra) + sqrt(-1) .* (zz -za);
		thl  = unwrap(angle(cl),[],2);
		thl(thl<0) = thl(thl<0) + 2 .* pi;
		rhol = abs(cl);
		[thl,indl] = sort(thl);
		rhol    =   rhol(indl);
		rhol = cat(2,rhol,rhol,rhol);
		thl = cat(2,thl -2.*pi,thl,thl+2.*pi);
		indnok = find(any(diff(thl,1,2)<=0,1));
		thl(:,indnok) =[];
		rhol(:,indnok)  = [];
		rho    = pchip(thl,rhol,teta);
		Rext   = ra + rho .* cos(teta);
		Zext   = za + rho .* sin(teta);			
		z0dinput.exp0d.Rsepa(k,:) = Rext;
		z0dinput.exp0d.Zsepa(k,:) = Zext - za;

		% securite
		rrmin = min(z0dinput.exp0d.Rsepa(:));
		if rrmin <= 0;
			z0dinput.exp0d.Rsepa = z0dinput.exp0d.Rsepa - rrmin + 0.01;
		end

		
	end
	
	rb0                    = abs(interp10d(data.EFIT.BVAC.t,data.EFIT.BVAC.data,temps,'nearest') .* 2.96);
	b0mem 		       = z0dinput.geo.b0;
	z0dinput.geo.b0        = rb0 ./ z0dinput.geo.R;
	fprintf('relative variation b0 : %g\n', std((b0mem - z0dinput.geo.b0) ./(rb0 ./ 2.96))); 
	
end

% reactivate ionisation of neutral beam by fast ions
z0dinput.option.fast_ion_sbp   = 1;

% parametrage optimal pour ILW
if flag_ilw == 1
     z0dinput.option.ane            = 11;
     z0dinput.option.Recycling      = 0.7;
     z0dinput.option.scaling        = 12;
     z0dinput.option.dilution       = 1;
     z0dinput.option.tau_limitation = 'On';
     z0dinput.option.ploss_exp      = 'max_power';
     z0dinput.option.l2hscaling     = 2;
     z0dinput.option.plhthr         = 'pel+pion';
     z0dinput.option.fpl2h_lim      = 2;
     z0dinput.option.usepped_scl    = 1;
     z0dinput.option.xiioxie        = -4.5;
     z0dinput.option.grad_ped       = 3;
     z0dinput.option.qdds           = 1;
     z0dinput.option.kidds          = 3;
     z0dinput.option.runaway        = 3;
     z0dinput.option.zeff           = 0;
     z0dinput.option.faccu          = 0;
     z0dinput.option.W_effect       = 1;
     z0dinput.option.density_model  = 'minconv';
     z0dinput.option.frad           = 1;
     z0dinput.option. matthews      = -1;
     z0dinput.option.z_prad         = 'Stangeby';
     z0dinput.option.gaunt          = 1;
     z0dinput.option.noncoronal     = 0;
     z0dinput.option.configuration  = 2;
     z0dinput.option.lambda_scale   = 3;
     z0dinput.option.factor_scale   = 1;
     z0dinput.option.eioniz         = 0;
     z0dinput.option.sol_model      = 'scaling';
     z0dinput.option.sol_rad        = 'decoupled';
     z0dinput.option.fcond          = -1;
     z0dinput.option.fmom           = 0;
     z0dinput.option.mach_corr      = 1;
     z0dinput.option.cw_factor      = 0;
     z0dinput.option.cw_offset      = 5e-5;
     z0dinput.option.e_shielding    = 'Honda-Sauter';
     z0dinput.option.sitb           = 0;
     z0dinput.option.mode_expo_inte = 1;
     z0dinput.option.cronos_regul   = 2;
     z0dinput.option.impur_rot      = 'max';
     z0dinput.option.available_flux = 20;
     z0dinput.option.hmore_pped     = 1;
     
end

if z0dinput.option.gaz == 3
    % securite fausse information provenant de JET
    z0dinput.cons.ftnbi  = 0 .* z0dinput.cons.iso;
end

% % close old MDS+ 
% % do nothing for recent version 
% status = mdsdisconnect();
% if mod( status, 2) == 0
% 		disp( 'MDS disconnect failed');
% end



function yy = zechan(x,y,xx,mode)

if isempty(x) | isempty(y)
	yy = sqrt(-1) .* ones(size(xx));
else
	yy =interp10d(x,y(:,end),xx,mode);
end


function yi = interp10d(x,y,xi,methode)

if (size(x,1)>1) & (size(x,2)==1)
	indnok = 0;
	nb     = 100;
	while (~isempty(indnok)) & (nb >0)
		indnok = find(diff(x)<=0);
		if ~isempty(indnok)
			x(indnok) = [];
			y(indnok,:) = [];
	
		end
		nb = nb - 1;
	end
end
yi = interp1(x,y,xi,methode);



% fonction d'interpolation pour consigene TS
function  yi = interp1c(x,y,xi)

	[x,ind] = sort(x);
	y       = y(ind);
	while (any(diff(x) <=0))
		indbad = find(diff(x) <=0);
		x(indbad+1) = x(indbad+1) + 1e-6;
	end
	yi = interp1(x,y,xi,'linear');




