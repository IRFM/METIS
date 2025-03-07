% selection of method depending on platform
function varargout = zmdsplusjet(varargin)

global watchdog_mdsplus

if isappdata(0,'CACHE_MDSPLUS_JET')
  [ok,out] = zcachemdsplusjet(varargin);
  if ok == 1
    if ~isempty(out)
         varargout =out;
         return
    end
  end
end

% put MDS+ in trouble sometime
%  if ~ispc
%      % time to close the session before crash of the connection
%      if isempty(watchdog_mdsplus)
%          watchdog_mdsplus = timer;
%          watchdog_mdsplus.StartDelay = 1000;
%          watchdog_mdsplus.BusyMode = 'queue';
%          watchdog_mdsplus.TimerFcn = @(myTimerObj, thisEvent)watchdog_mdsplus_callback(myTimerObj, thisEvent);
%          start(watchdog_mdsplus);
%      end
%  end
% if isappdata(0,'MDSPLUS_USERNAME')
%     if ~isempty(strfind(varargin{2},'jpf')) || ispc
%         [varargout{1:nargout}] =zmdsplusjet_old(varargin{:});
%     else
%         [varargout{1:nargout}] =zmdsplusjet_fast(varargin{:});
%     end
% else
    [varargout{1:nargout}] =zsaljet(varargin{:});
%end
% shot    = numero du choc
% seq     = sequence (vide =  derniere)
% uid     = user identifier (defaut = ppf)
% name    = nom  du signal
% connect = 0 -> connection et deconnection, 1 -> connection seule, 2 -> deconnection seule, 3 -> rien
% Example : [data,time,space,units,cr,flag_lec] = zmdsplusjet(53521,'PPF/MAGN/IPLA');
% Update record:
% YYMMDD Who Comments
% 080521 AJC Complete rewrite. Note still assumes PPF/... Does not put client or server MDSPlus layers
%            into error handlers (unless real error) to get data.
zcachemdsplusjet(varargin,varargout);

function [data,time,space,units,cr,flag_lec] = zmdsplusjet_fast(shot,name,seq,uid,connect)

% reservation des sorties

data     =    [];
time     =    [];
space    =    [];
units    =    [];
cr       =    0;
flag_lec =    0;

%
% gestion des entrees
%

if nargin < 2
   	error('syntaxe :  [data,time,space,units,cr,flag_lec] = zmdsplusjet(shot,name,seq,uid,connect)');
end

if isempty( shot)
   	error('il faut donner un numero de choc');
end	

if isempty( name)
   	error('il faut donner le nom du signal');
end	

%
% Missing sequence, use 0.
%

if nargin < 3
  	seq = 0;
end
if isempty( seq)	
	seq = 0;
end

%
% Missing or empty uid, use JETPPF.
%

if nargin < 4
  	uid = '';
end
if isempty( uid)
	uid = 'jetppf';
end

%
% Missing or empty connect, use 0.
%

if nargin < 5 
  	connect = 0;
end
if isempty( connect)
  	connect = 0;
end

%
% Ensure connected. Note in original mdsconnect it never returns an error, but aborts.
%

if (connect == 0) || (connect == 1) || (isempty(strmatch('mdsipmex',inmem,'exact')) &&  ~isempty(which('mdsipmex')))
    % bypass error appening on some system where the first try return an error 
    try
	status = mdsconnect(sprintf('ssh://%s@mdsplus.jetdata.eu',getappdata(0,'MDSPLUS_USERNAME')));
    catch
	status = mdsconnect(sprintf('ssh://%s@mdsplus.jetdata.eu',getappdata(0,'MDSPLUS_USERNAME')));    
    end
    if ~isempty( status) && (status ~= 1)
           error( 'Connection failed to mdsplus.jet.efda.org');
            cr = -1;
            return
    end
elseif isunix
            [s,t] = unix('lsof -i | grep mdsplus.jet.efda.org | grep -i ESTABLISHED | wc -l');
            if s ~= 0 
		% bypass error appening on some system where the first try return an error 
		try                
		    status = mdsconnect(sprintf('ssh://%s@mdsplus.jetdata.eu',getappdata(0,'MDSPLUS_USERNAME')));
		catch
		    status = mdsconnect(sprintf('ssh://%s@mdsplus.jetdata.eu',getappdata(0,'MDSPLUS_USERNAME')));		
		end
                if ~isempty( status) && (status ~= 1)
                    error( 'Connection failed to mdsplus.jet.efda.org');
                    cr = -1;
                    return
               end

            end 
end

%
% Set ppfuid, cache the signal (may be none), geting ppf error code, return this error code. Note
% use substitution in tdi expresion and its a multi-part one.
%

[ier status] = mdsvalue( 'ppfuid($),_s=jet($,long($),_ier),_ier', uid, sprintf( '%s/%d', name, fix( seq)), fix( shot));

%
% Status is odd if the mds level worked. Normally ier will be an MDS message if error. Note this
% is used below - when we read data it may contain an error message if fault occured.
%

if mod( status, 2) == 0
	disp( sprintf( '%d/%s/%d (%s) MDS error: %d - no data', fix( shot), name, fix( seq), uid, fix( status)));
	disp( ier);
	cr = -3;
	return
end

%
% Check ppf error code. There is only data if its zero. The ppf interface returns
% {data,t} or {data,x,t}. Attempt to get an error message assuming its ppfget related.
%

if ier ~= 0
	[msg status] = mdsvalue( 'ppfemsg("ppfget",long($))', fix( ier));
	disp( sprintf( '%d/%s/%d (%s) PPF error: %d "%s" - no data', fix( shot), name, fix( seq), uid, ier, msg));
	cr = -10001;
	return
end

%
% Get data from saved signal.
%

[data status] = mdsvalue( 'data(_s)');

if mod( status, 2) == 0
	disp( sprintf( '%d/%s/%d (%s) MDS error: %d tdi: data(_s)', fix( shot), name, fix( seq), uid, fix( status)));
	disp( data);
	data = [];
	cr   = -4;
	return
end

%
% Get units from saved signal.
%

[units status] = mdsvalue( 'units(_s)');

if mod( status, 2) == 0
	disp( sprintf( '%d/%s/%d (%s) MDS error: %d tdi: units(_s)', fix( shot), name, fix( seq), uid, fix( status)));
	disp( units);
	data  = [];
	units = [];
	cr    = -5;
	return
end

%
% If we use ppf function and it has data then its a signal with n descriptors, get this value as its 
% needed to avoid mdsplus errors so that we get only the bits present.
%

[ndesc status] = mdsvalue( 'ndesc(_s)');

if mod( status, 2) == 0
	disp( sprintf( '%d/%s/%d (%s) MDS error: %d tdi: ndesc(_s)', fix( shot), name, fix( seq), uid, fix( status)));
	disp( ndesc);
	data  = [];
	units = [];
	cr    = -6;
	return
end

%
% Get x or t from saved signal in dimension 0.
%

if ndesc > 3
	[space status] = mdsvalue( 'dim_of(_s,0)');
else
	[time status] = mdsvalue( 'dim_of(_s,0)');
end

if mod( status, 2) == 0
	disp( sprintf( '%d/%s/%d (%s) MDS error: %d tdi: dim_of(_s,0)', fix( shot), name, fix( seq), uid, fix( status)));
	if ndesc > 3
		disp( space);
	else
		disp( time);
	end
	data  = [];
	units = [];
	space = [];
	time  = [];
	cr    = -7;
	return
end

%
% Get t from saved signal in dimension 1 if had x.
%

if ndesc > 3
	[time status] = mdsvalue( 'dim_of(_s,1)');
	if mod( status, 2) == 0
		disp( sprintf( '%d/%s/%d (%s) MDS error: %d tdi: dim_of(_s,1)', fix( shot), name, fix( seq), uid, fix( status)));
		disp( time);
		data  = [];
		units = [];
		space = [];
		time  = [];
		cr    = -8;
		return
	end
end

%
% Final adjustments, this logic was in original.
%

if ~isempty( space)
	if size( space, 2) == 1
		space = space';
	end
end
if ~isempty( time)
	if size( time, 1) == 1
		time = time';
	end
	if ~isempty( data)
		if size( time, 1) ~= size( data, 1)
			data = data';
		end
	end
end

flag_lec = 1;

disp( sprintf( '%d/%s/%d (%s) data read', fix( shot), name, fix( seq), uid));

%
% This does nothing, but is here in case run against old mdsdisconnect.
%
%  if ~isdeployed
%    if (connect == 0) || (connect == 2)
%  	  status = mdsdisconnect();
%  	  if mod( status, 2) == 0
%  		  disp( 'MDS disconnect failed');
%  	  end
%    end
%  end
% shot    = numero du choc
% seq     = sequence (vide =  derniere)
% uid     = user identifier (defaut = ppf)
% name    = nom  du signal
% connect = 0 -> connection et deconnection, 1 -> connection seule, 2 -> deconnection seule, 3 -> rien
function [data,time,space,units,cr,flag_lec] = zmdsplusjet_old(shot,name,seq,uid,connect)

% reservation des sorties
data     =    [];
time     =    [];
space    =    [];
units    =    [];
cr       =    0;
flag_lec =    0;
% gestion des entrees
if nargin < 2
   error('syntaxe :  [data,time,space,units,cr] = zmdsplujet(shot,name,seq,uid,connect)');
end
if isempty(shot)
   error('il faut donner un numero de choc');
end	
if isempty(name)
   error('il faut donner le nom du signal');
end	
if nargin < 3
  seq = [];
elseif seq == 0
  seq = [];
end
if nargin < 4
    uid = '';
end
if nargin < 5
    connect = 0;
elseif isempty(connect)
    connect = 0;
end

% debut de la transaction
if (connect == 0) | (connect == 1) || (isempty(strmatch('mdsipmex',inmem,'exact')) &&  ~isempty(which('mdsipmex')))
    try
        % bypass error appening on some system where the first try return an error
        try
            cr = mdsconnect(sprintf('ssh://%s@mdsplus.jetdata.eu',getappdata(0,'MDSPLUS_USERNAME')));
        catch
            cr = mdsconnect(sprintf('ssh://%s@mdsplus.jetdata.eu',getappdata(0,'MDSPLUS_USERNAME')));
        end
        if isempty(cr)
            cr =0;
        end
    catch
        cr = -1;
    end
    if cr < 0
        
        error('Probleme de connexion a mdsplus.jet.efda.org');
        return
    end
elseif isunix
    [s,t] = unix('lsof -i | grep mdsplus.jet.efda.org | grep -i ESTABLISHED | wc -l');
    if s ~= 0
        try
            % bypass error appening on some system where the first try return an error
            try
                cr = mdsconnect(sprintf('ssh://%s@mdsplus.jetdata.eu',getappdata(0,'MDSPLUS_USERNAME')));
            catch
                cr = mdsconnect(sprintf('ssh://%s@mdsplus.jetdata.eu',getappdata(0,'MDSPLUS_USERNAME')));
            end
            if isempty(cr)
                cr =0;
            end
        catch
            cr = -1;
		end
		if cr < 0
		      
			error('Probleme de connexion a mdsplus.jet.efda.org');
			return
		end
            end 
end


% choix du user
if ~isempty(uid)
   try
	   cr = mdsvalue(sprintf('_sig = ppfuid("%s")',upper(uid)));
           if isempty(cr) 
		cr =0;
	   end
		if isstr(cr)
		   disp(cr)
			cr = -2;
		end
	catch
			cr = -2;
	end 
	if cr ~= 0
		disp('Probleme dans le uid');
		return
	end
else
   try
	   cr = mdsvalue('_sig = ppfuid("JETPPF")');
           if isempty(cr) 
		cr =0;
	   end
		if isstr(cr)
		   disp(cr)
			cr = -2;
		end
	catch
			cr = -2;
	end 
         if cr ~= 0
		disp('Probleme dans le uid');
		return
	end
end


% lecture des donnees
if isempty(seq)
	namec = name;
else
   namec = sprintf('%s/%d',name,fix(seq));
end
try
	[data,crd]  = mdsvalue(sprintf('_sig = jet("%s",%d)',namec,shot));
	if (crd == 1) & (all(size(data) == 1)) & (data(1,1) == 0)
	   data = [];
		if isempty(uid)
		    if isempty(seq)
	    		fprintf('pas de donnees pour %s @ %d\n',name,shot);
			 else
	    		fprintf('pas de donnees pour %s @ %d.%d\n',namec,shot,seq);
			 end
		else
		    if isempty(seq)
    			fprintf('pas de donnees pour %s: %s @ %d\n',uid,name,shot);
			 else
    			fprintf('pas de donnees pour %s: %s @ %d.%d\n',uid,namec,shot,seq);
			 end
		end
	else
		time       = mdsvalue('dim_of(_sig,0)');
		if ~isempty(time)
			space   = mdsvalue('dim_of(_sig,1)');
	   else
			time    = mdsvalue('dim_of(_sig)');
		end
		units    = mdsvalue('UNITS(_sig)');
		if ~isempty(space)
		   if ~isstr(space)
		   	void  = time;
				time  = space;
 	     	 	data  = data;
				space = void;
		   	if size(space,2) == 1
			        space = space';
				end
			else
			   space = [];
			end
		end
		if size(time,1) == 1
			     time =time';
	   end
		if size(time,1) ~= size(data,1)
				     data =data';
		end
		if isempty(uid)
	    		fprintf('Lecture de %s @ %d\n',name,shot);
                        flag_lec = 1;
		else
    		        fprintf('Lecture de %s: %s @ %d\n',uid,name,shot);
                        flag_lec = 1;
		end
	end
catch
	if isempty(uid)
    		fprintf('Probleme de lecture de %s@%d\n',name,shot);
	else
    		fprintf('Probleme de lecture de %s:%s@%d\n',uid,name,shot);
	end
	 cr = -10001
end

% on remet le user par defaut
if ~isempty(uid)
   try
	   cr = mdsvalue('_sig = ppfuid("JETPPF")');
           if isempty(cr) 
		cr =0;
	   end
		if isstr(cr)
		   disp(cr)
			cr = -3;
		end
	catch
			cr = -3;
	end 
	if cr ~= 0
		disp('Probleme dans le uid PPF');
		return
	end
end

%  % fin de la transaction
%  if (connect == 0) | (connect == 2)
%   	try
%   			cr = mdsdisconnect;
%     catch
%  			cr = -5;
%  	end
%             if isempty(cr) 
%  		cr =0;
%  	   end
%  	if cr ~= 1
%  		disp('Probleme de la deconnexion connexion de mdsplus.jet.efda.org');
%  		return
%  	end
%  end

% si on arrive ici 
if crd == 3
   cr = 0;
else
	cr = crd;
end

function  [ok,out] = zcachemdsplusjet(in,data)

out = {};
ok = 0;
if nargin == 2
    if isappdata(0,'CACHE_MDSPLUS_JET')
        cache = getappdata(0,'CACHE_MDSPLUS_JET');
    else
        cache = {};
    end
    cache{end+1}.in = in;
    cache{end}.out  = data;
    setappdata(0,'CACHE_MDSPLUS_JET',cache);
    if isappdata(0,'PATH2CACHE_MDSPLUS_JET')
        inname = getinname(getappdata(0,'PATH2CACHE_MDSPLUS_JET'),in);
        save(inname,'in','out');
    end
    
else
    cache = getappdata(0,'CACHE_MDSPLUS_JET');
    for k=1:length(cache)
        inc = cache{k}.in;
        if length(inc)  == length(in)
            ok = 1;
            for l = 1:min(4,length(inc))
                if ~same(inc{l},in{l})
                    ok = 0;
                end
            end
            if ok == 1
                out = cache{k}.out;
                if ~isempty(out)
		  fprintf('%d/%s retrieved from cache\n',cache{k}.in{1},cache{k}.in{2});
		else 
		  ok = 0;
		end
                break;
            end
        end
    end
    if (ok == 1) && isappdata(0,'PATH2CACHE_MDSPLUS_JET')
        inname = getinname(getappdata(0,'PATH2CACHE_MDSPLUS_JET'),in);
        if ~exist(sprintf('%s.mat',inname),'file')
            save(inname,'in','out');
        end
    end
    
    if (ok == 0) && isappdata(0,'PATH2CACHE_MDSPLUS_JET')
        inname = getinname(getappdata(0,'PATH2CACHE_MDSPLUS_JET'),in);
        if exist(sprintf('%s.mat',inname),'file')
            data = load(inname);
            out = data.out;
            if ~isempty(out)
		fprintf('%d/%s retrieved from cache\n',in{1},in{2});
                ok = 1;		
		if isappdata(0,'CACHE_MDSPLUS_JET')
		    cache = getappdata(0,'CACHE_MDSPLUS_JET');
		else
		    cache = {};
		end
		cache{end+1}.in = in;
		cache{end}.out  = out;
		setappdata(0,'CACHE_MDSPLUS_JET',cache);
	    end
	end
            
     end
     
    
end

function inname = getinname(path,in)

inname = 'CACHE_MDSPLUS_JET';
for k=1:min(4,length(in))
  if isempty(in{k})
    % rien
  elseif isnumeric(in{k})
    inname = sprintf('%s_%g',inname,in{k});
  else
     s = in{k};
     s = strrep(s,' ','_');
     s = strrep(s,'+','_plus_');
     s = strrep(s,'?','_ptq_');
     s = strrep(s,'*','_star_');
     s = strrep(s,'=','_equal_');
     s = strrep(s,':','_2pt_');
     s = strrep(s,'/','_slash_');
     s = strrep(s,'\','_anti_');   
     s = strrep(s,'(','_pg_');
     s = strrep(s,')','_pd_');   
     s = strrep(s,'[','_cg_');
     s = strrep(s,']','_cd_');   
     s = strrep(s,'{','_cg_');
     s = strrep(s,'}','_cd_');   
     inname = sprintf('%s_%s',inname,s); 
  end

end
inname = fullfile(path,inname);


%[data,time,space,units,cr,flag_lec] = zsaljet(shot,name,seq,uid,connect)
function [data,time,space,units,cr,flag_lec] =zsaljet(shot,name,seq,uid,connect)

% initialisation
data     = [];
time     = [];
space    = [];
units    = '';
cr       = 0;
flag_lec = false;

if (nargin < 3) || isempty(seq)
    seq = 0;
end
data_path = sprintf('data/pulse/%d/ppf/signal/jet%s',shot,name);
if ~isempty(uid)
  data_path = strrep(data_path,'jetppf',uid);
end
[rep,raw_text] = sal_get(data_path,'JET',seq);
if isempty(rep) && ~isempty(strfind(name,'jpf'))
    data_path = sprintf('data/pulse/%d/%s',shot,name);
    [rep,raw_text] = sal_get(data_path,'JET',seq);
end
if ~isempty(strfind(upper(name),'PPF/CX')) && isempty(rep)
    data_path = strrep(data_path,'jetppf','cxsbatch');
    [rep,raw_text] = sal_get(data_path,'JET',seq);
end

fprintf('Reading via SAL: %s ',data_path);

raw_text(raw_text < ' ') = [];

if isempty(rep)
    cr = -1;
    fprintf('%s\n',raw_text);
    return
elseif ~isfield(rep.object,'data')
    fprintf('%s\n',raw_text);
    return
end

flag_lec = true;
data = rep.object.data.value.data;
if isfield(rep.object,'dimensions') && isfield(rep.object.dimensions.value,'x0')
    time = rep.object.dimensions.value.x0.value.data.value.data;
end
if isfield(rep.object,'dimensions') && isfield(rep.object.dimensions.value,'x1')
    space = rep.object.dimensions.value.x1.value.data.value.data;
end
if isfield(rep.object,'units')
    units = rep.object.units.value;
end
ss = size(data);
fprintf('[%d %d] \n',ss(1),ss(2));



