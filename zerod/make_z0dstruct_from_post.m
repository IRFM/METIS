% initialize z0dstruct for zerodevolution from METIS simulation post data
% structure. 
% Syntax:
%   z0dstruct = make_z0dstruct_from_post(post,time);
%
% inputs:
%     post = METIS output data structure;
%     time = strating point for the simulation (s)
%
% output: 
%     zs        = zerod data for time
%     profil    = profil for time
%     z0dstruct = initialized data structure for z0devolution
%
function [zs,profil,z0dstruct] = make_z0dstruct_from_post(post,time)

% extract data from post
time0 = post.z0dinput.cons.temps;
ind0 = find(post.z0dinput.cons.temps >= time,1);
if isempty(ind0)
    ind0  = length(post.z0dinput.cons.temps) -1;
elseif ind0 == length(post.z0dinput.cons.temps)
    ind0  = length(post.z0dinput.cons.temps) -1;
end    
ind4t = (ind0-2):(ind0+1);
ind4t(ind4t < 1) = 1;
noms = fieldnames(post.zerod);
for k=1:length(noms)
    if length(post.zerod.(noms{k})) == length(time0)
        z0dstruct.zs.(noms{k})      = post.zerod.(noms{k})(ind4t);
    else
        z0dstruct.zs.(noms{k})      = post.zerod.(noms{k});
    end
end
noms = fieldnames(post.z0dinput.cons);
for k=1:length(noms)
    if length(post.z0dinput.cons.(noms{k})) == length(time0)
        z0dstruct.z0dinput.cons.(noms{k})      = post.z0dinput.cons.(noms{k})(ind4t);
    else
        z0dstruct.z0dinput.cons.(noms{k})      = post.z0dinput.cons.(noms{k});
    end
end
noms = fieldnames(post.z0dinput.geo);
for k=1:length(noms)
    if length(post.z0dinput.geo.(noms{k})) == length(time0)
        z0dstruct.z0dinput.geo.(noms{k})      = post.z0dinput.geo.(noms{k})(ind4t);
    else
        z0dstruct.z0dinput.geo.(noms{k})      = post.z0dinput.geo.(noms{k});
    end
end
noms = fieldnames(post.z0dinput.exp0d);
for k=1:length(noms)
    switch noms{k}
        case {'Rsepa','Zsepa'}
            %skip
        otherwise
            if length(post.z0dinput.exp0d.(noms{k})) == length(time0)
                try
                    z0dstruct.z0dinput.exp0d.(noms{k})      = post.z0dinput.exp0d.(noms{k})(ind4t);
                catchzs
                    keyboard
                end
            else
                z0dstruct.z0dinput.exp0d.(noms{k})      = post.z0dinput.exp0d.(noms{k});
            end
    end
end
%  SEPA
if isfield(post.z0dinput.exp0d,'Rsepa') || isfield(post.z0dinput.exp0d,'Zsepa')
    % coherence test
    if isfield(post.z0dinput.exp0d,'Rsepa') && ~isfield(post.z0dinput.exp0d,'Zsepa')
        post.z0dinput.exp0d = rmfield(post.z0dinput.exp0d,'Rsepa');
    end
    if ~isfield(post.z0dinput.exp0d,'Rsepa') && isfield(post.z0dinput.exp0d,'Zsepa')
        post.z0dinput.exp0d = rmfield(post.z0dinput.exp0d,'Zsepa');
    end
    if all(size(post.z0dinput.exp0d.Rsepa) ~= size(post.z0dinput.exp0d.Zsepa))
        post.z0dinput.exp0d = rmfield(post.z0dinput.exp0d,'Rsepa');
        post.z0dinput.exp0d = rmfield(post.z0dinput.exp0d,'Zsepa');
    end
    if (size(post.z0dinput.exp0d.Rsepa) ~= length(time0))
        post.z0dinput.exp0d = rmfield(post.z0dinput.exp0d,'Rsepa');
        post.z0dinput.exp0d = rmfield(post.z0dinput.exp0d,'Zsepa');
    end
    post.z0dinput.exp0d.Rsepa = post.z0dinput.exp0d.Rsepa(ind4t,:);
    post.z0dinput.exp0d.Zsepa = post.z0dinput.exp0d.Zsepa(ind4t,:);

end
%
z0dstruct.z0dinput.option = post.z0dinput.option;
post.z0dinput.option = 1;
if isfield(post.z0dinput,'info')
    z0dstruct.z0dinput.zsinfo = post.z0dinput.zsinfo;
    z0dstruct.info            = post.z0dinput.zsinfo;
else
    z0dstruct.z0dinput.zsinfo = zerod_0dinfo;
    z0dstruct.info     = zerod_0dinfo;
end
z0dstruct.z0dinput.langue = post.z0dinput.langue;
if isfield(post.z0dinput,'profinfo')
    z0dstruct.z0dinput.profinfo = post.z0dinput.profinfo;
else
    z0dstruct.z0dinput.profinfo = z0dprofinfo;
end
z0dstruct.z0dinput.mode_exp = -1;
z0dstruct.z0dinput.machine  = post.z0dinput.machine;
z0dstruct.z0dinput.shot     = post.z0dinput.shot;
z0dstruct.info              = zerod_0dinfo;

%  profiles
noms = fieldnames(post.profil0d);
time1d = post.profil0d.temps;
for k=1:length(noms)    
    if size(post.profil0d.(noms{k}),1) == length(time1d)
        z0dstruct.profil.(noms{k}) = interp1(time1d,post.profil0d.(noms{k}),z0dstruct.z0dinput.cons.temps,'linear','extrap');
    else
        z0dstruct.profil.(noms{k}) = post.profil0d.(noms{k});
    end
end

% rebased time
new_time = z0dstruct.z0dinput.cons.temps + time - z0dstruct.z0dinput.cons.temps(end);
dtime = diff(new_time);
dtime(dtime <= 0) = mean(dtime(dtime > 0));
time_inter = cat(1,0,cumsum(dtime(:)));
new_time = time_inter + new_time(end) - time_inter(end);
z0dstruct.z0dinput.cons.temps = new_time;
z0dstruct.zs.temps = new_time;
z0dstruct.profil.temps = new_time;

% zs
noms = fieldnames(z0dstruct.zs);
zs   = [];
for k = 1:length(noms)
    zs.(noms{k})  = z0dstruct.zs.(noms{k})(end);
end

% profil 
noms   = fieldnames(z0dstruct.profil);
profil = [];
for k = 1:length(noms)
    profil.(noms{k})  = z0dstruct.profil.(noms{k})(end,:);
end


