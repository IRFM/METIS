function [sepa1t,z0dstruct,nbpts] = manage_LCFS_nbpoints(sepa1t,z0dstruct,nbpts_metis,nbpts_equi)

% test is something should be done 
if isempty(sepa1t)
    nbpts_sepa1t = NaN;
elseif isempty(sepa1t.Rsepa)
    nbpts_sepa1t = NaN;
else
    nbpts_sepa1t = size(sepa1t.Rsepa,2);
end
if isempty(z0dstruct) ||  ...
   ~isfield(z0dstruct,'z0dinput') || ...
   ~isfield(z0dstruct.z0dinput,'exp0d')
    nbpts_exp0d = NaN;
elseif ~isfield(z0dstruct.z0dinput.exp0d,'Rsepa')
    nbpts_exp0d = NaN;
elseif isempty(z0dstruct.z0dinput.exp0d.Rsepa)
    nbpts_exp0d = NaN;  
else
    nbpts_exp0d = size(z0dstruct.z0dinput.exp0d.Rsepa,2);
end
if isempty(z0dstruct) ||  ...
   ~isfield(z0dstruct,'profil') 
    nbpts_prof = NaN;
elseif ~isfield(z0dstruct.profil,'Rsepa')
    nbpts_prof = NaN;
elseif isempty(z0dstruct.profil.Rsepa)
    nbpts_prof = NaN;  
else
    nbpts_prof = size(z0dstruct.profil.Rsepa,2);
end
if isempty(z0dstruct) ||  ...
   ~isfield(z0dstruct,'profil0d')
    nbpts_prof0d = NaN;
elseif ~isfield(z0dstruct.profil0d,'Rsepa')
    nbpts_prof0d = NaN;
elseif isempty(z0dstruct.profil0d.Rsepa)
    nbpts_prof0d = NaN;  
else
    nbpts_prof0d = size(z0dstruct.profil0d.Rsepa,2);
end

% all lengths
nbpts = cat(1,nbpts_sepa1t,nbpts_exp0d,nbpts_metis,nbpts_equi,nbpts_prof,nbpts_prof0d);
nbpts = max(nbpts(isfinite(nbpts)));
% For METIS internal resampling use with external data
setappdata(0,'LCFS_nbpts',nbpts);


% use the longest size
if isfinite(nbpts_sepa1t) && (size(sepa1t.Rsepa,2) ~= nbpts)
    % resize sepa1t
    [sepa1t.Rsepa,sepa1t.Zsepa] = metis4muscle3_resize_LCFS(sepa1t.Rsepa,sepa1t.Zsepa,nbpts);
end
if isfinite(nbpts_exp0d) && (size(z0dstruct.z0dinput.exp0d.Rsepa,2) ~= nbpts)
    % resize previous LCFS
    [z0dstruct.z0dinput.exp0d.Rsepa,z0dstruct.z0dinput.exp0d.Zsepa] = metis4muscle3_resize_LCFS(z0dstruct.z0dinput.exp0d.Rsepa,z0dstruct.z0dinput.exp0d.Zsepa,nbpts);
end
if isfinite(nbpts_prof) && (size(z0dstruct.profil.Rsepa,2) ~= nbpts)
     [z0dstruct.profil.Rsepa,z0dstruct.profil.Zsepa] = metis4muscle3_resize_LCFS(z0dstruct.profil.Rsepa,z0dstruct.profil.Zsepa,nbpts);   
end
if isfinite(nbpts_prof0d) && (size(z0dstruct.profil0d.Rsepa,2) ~= nbpts)
     [z0dstruct.profil0d.Rsepa,z0dstruct.profil0d.Zsepa] = metis4muscle3_resize_LCFS(z0dstruct.profil0d.Rsepa,z0dstruct.profil0d.Zsepa,nbpts);   
end

% figure(21);clf
% plot(sepa1t.Rsepa',sepa1t.Zsepa','.',exp0d.Rsepa',exp0d.Zsepa','-')
% drawnow



