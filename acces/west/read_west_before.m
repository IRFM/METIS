% lecture des donnees WEST avant choc
function wf = read_west_before(shot)

wf = [];
todo = 1;
if nargin == 0
  shot = 0;
elseif isempty(shot)
  shot = 0;
end

dir_mem = pwd;
dir_tmp = tempname;
[s,m] = mkdir(dir_tmp);
if s ~= 1
  disp(m);
  return
end
cd(dir_tmp);
if (shot == 0) && isdir('/Home/cgc/WEST_PCS')
	mkdir('DCS/Operation/DocumentRoot/Params')
	[s,t] = unix('cp -rp /Home/cgc/WEST_PCS/* ./DCS/Operation/DocumentRoot/Params');
        if s == 0
            [s,t] = unix('cat ./DCS/Operation/DocumentRoot/Params/shot_number.txt');
            if  s == 0
	      shot = str2num(t);
	      dsup = dir('./DCS/Operation/DocumentRoot/Params/Sup.xml');
	      ddp  = dir('./DCS/Operation/DocumentRoot/Params/DP.xml');
	      if (length(dsup) > 0) && (length(ddp) > 0)
		  disp('DCS data before shot recovered');
		  todo = 0;
	      end
	  end
       end
       if todo == 1
	  shot = 0;
       end
end
if (shot == 0) && isdir('/donnees/WEST_PCS/')
	mkdir('DCS/Operation/DocumentRoot/Params')
	[s,t] = unix('cp -rp /donnees/WEST_PCS/* ./DCS/Operation/DocumentRoot/Params');
        if s == 0
            [s,t] = unix('cat ./DCS/Operation/DocumentRoot/Params/shot_number.txt');
            if  s == 0
	      shot = str2num(t);
	      dsup = dir('./DCS/Operation/DocumentRoot/Params/Sup.xml');
	      ddp  = dir('./DCS/Operation/DocumentRoot/Params/DP.xml');
	      if (length(dsup) > 0) && (length(ddp) > 0)
		  disp('DCS data before shot recovered');
		  todo = 0;
	      end
	  end
       end
end
if todo == 1
    if shot == 0
	shot = tsdernier_choc;
    end
    if exist('tsrfile') && (todo == 1) 
	sr=tsrfile(shot, 'FPCSPARAM','param.file'); 
    else
	return;
    end
    if sr ~= 0
      cd(dir_mem);
      rmdir(dir_tmp,'s');
      disp('unable to access to pre shot data');
      return
    end
    % compression option change depending on shot
    [s,t] = unix('tar xzvf param.file');
    if s~= 0
	[s,t] = unix('tar xjvf param.file');
    end
    if s~= 0
	[s,t] = unix('tar xvf param.file');
    end
end
if s~= 0
  cd(dir_mem);
  rmdir(dir_tmp,'s');
  disp(t);
  return
end

% copy python script to local directory
[s,t] = unix(sprintf('cp -rp %s .',fullfile(fileparts(which('read_west_before')),'py','*.py')));
if s ~= 0
  cd(dir_mem);
  rmdir(dir_tmp,'s');
  disp(t);
  return
end

%sup = xml_read('DCS/Operation/DocumentRoot/Params/Sup.xml');
%dp  = xml_read('DCS/Operation/DocumentRoot/Params/DP.xml');
% call python script for data decoding
[s,t] = unix('python3 west_before.py');
if s~= 0
   [s,t] = unix('python west_before.py');
end
if s~= 0
   [s,t] = unix('module purge;  module load anaconda/python27; python west_before.py');
end

if s ~= 0
  cd(dir_mem);
  rmdir(dir_tmp,'s');
  disp(t);
  return
end
% retrieve data to Matlab
tt = strrep(sprintf('%s;',strrep(t,sprintf('\n'),',')),'],','];');
index_start = min(findstr(tt,'wf.'));
index_end   = max(findstr(tt,'];')) + 1;
try
    eval(tt(index_start:index_end));
catch
    indstart = strfind(t,'(''wf.');
    indend   = strfind(t,sprintf(')\n'));
    for k=1:min(length(indstart),length(indend))
        tt=t(indstart(k)+2:indend(k)-1);
        tt = strrep(tt,'='', array(','=');
        tt = strrep(tt,'])','];');
        tt = strrep(tt,sprintf('\n'),'');
        tt = strrep(tt,', dtype=float64)','');
        eval(sprintf('%s;',tt));
    end
end
  
% end of the function  
cd(dir_mem);
rmdir(dir_tmp,'s');

% assign shot number
wf.shot = shot;