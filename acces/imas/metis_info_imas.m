% return information about version in repository
function [repository,commit,version] = metis_info_imas

pwd_mem = pwd;
cd(fileparts(which('metis')));

if isdir(fullfile(fileparts(which('metis')),'.svn')) || isdir(fullfile(fileparts(which('metis')),'..','.svn'))
  % this is a svn version
  [s,t] = unix('svn info --show-item url');
  if s ~= 0
      [s,t] = unix('svn info | grep -i URL');
      if s == 0
	[v,t]= strtok(t,':');
	if ~isempty(t)
	  t = t(3:end);
	end
      end
  end
  if s == 0
        tt = tseparec(t);
        repository = deblank(tt(end,:));
  else
        repository = '??????';
  end
  [s,t] = unix('svn info --show-item revision');
  if s ~= 0
      [s,t] = unix('svn info | grep -i Revision');
      if s == 0
	[v,t]= strtok(t,':');
	if ~isempty(t)
	  t = t(3:end);
	end
      end
  end
  if s == 0
        tt = tseparec(t);
        commit = deblank(tt(end,:));
  else
        commit = '??????';
  end
  [s,t] = unix('svn info --show-item kind');
  if s ~= 0
      [s,t] = unix('svn info | grep -i "Node Kind"');
      if s == 0
	[v,t]= strtok(t,':');
	if ~isempty(t)
	  t = t(3:end);
	end
      end
  end
  if s == 0
        tt = tseparec(t);
        version = sprintf('%d @ %s',zinebversion,deblank(tt(end,:)));
  else
        version = sprintf('%d @ %s',zinebversion,'??????');  
  end
elseif isdir(fullfile(fileparts(which('metis')),'.git')) || isdir(fullfile(fileparts(which('metis')),'..','.git'))
  %git version
  [s,t] = unix('git remote -v');
  if s == 0
        tt = tseparec(t);
        repository = deblank(tt(end,:));
  else
        repository = '??????';
  end
  [s,t] = unix('git rev-parse  HEAD');
  if s == 0
        tt = tseparec(t);
        commit = deblank(tt(end,:));
  else
        commit = '??????';
  end
  [s,t] = unix('git describe --all');
  if s == 0
        tt = tseparec(t);
        version = sprintf('%d @ %s',zinebversion,deblank(tt(end,:)));
  else
        version = sprintf('%d @ %s',zinebversion,'??????');  
  end
else
   repository = '??????'
   commit  = '??????'
   version = '??????'
end

cd(pwd_mem);