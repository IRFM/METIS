% recompute UAL data and put it in matlab file
% usefull to rebuilt 2d equilibrium data
% 
[error_flag,cpos_data_like] = metis4itm(post.z0dinput.shot,1,'',post,[]);
if isappdata(0,'METIS_FILENAME')
	[p,name_root]= fileparts(getappdata(0,'METIS_FILENAME'));
else
	name_root = 'noname';
end
filename = sprintf('%s_cpos_data_like4%s',name_root,'unknown_tokamak');
if isfield(post.z0dinput,'machine')
	if ~isempty(post.z0dinput.machine)
		filename = sprintf('%s_cpos_data_like4%s',name_root,post.z0dinput.machine);
	end
end
if isfield(post.z0dinput,'shot')
	if ~isempty(post.z0dinput.shot)
		filename = sprintf('%s@%d',filename,post.z0dinput.shot);
	end
end
save(filename,'error_flag','cpos_data_like');
fprintf('CPOs data like structure have been save in %s.mat file\n',filename);
