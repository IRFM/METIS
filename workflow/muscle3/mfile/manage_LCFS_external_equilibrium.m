function  data = manage_LCFS_external_equilibrium(data,LCFS_nbpts)

% init 
equi_loc = [];
to_be_updated = false;

if isfield(data,'equilibrium')
    equi_loc = data.equilibrium;
elseif isfield(data,'current_diffusion')
    equi_loc = data.current_diffusion;
end

if ~isempty(equi_loc)
    for k=1:length(equi_loc.time_slice)
        Rsepa = equi_loc.time_slice{k}.boundary.outline.r;
        Zsepa = equi_loc.time_slice{k}.boundary.outline.z;
        if ~isempty(Rsepa) && length(Rsepa) ~= LCFS_nbpts
            [Rsepa_out,Zsepa_out] = metis4muscle3_resize_LCFS(Rsepa(:)',Zsepa(:)',LCFS_nbpts);
            equi_loc.time_slice{k}.boundary.outline.r = Rsepa_out;
            equi_loc.time_slice{k}.boundary.outline.z = Zsepa_out;
            to_be_updated = true;
        end
    end
    if to_be_updated
        if isappdata(0,'EQUILIBRIUM_EXP') || isappdata(0,'CURDIFF_EXP')
            create_external_equilibrium_from_imas([],0,equi_loc);
        end        
        if isfield(data,'equilibrium')
            data.equilibrium = equi_loc;
        end
        if isfield(data,'current_diffusion')
            data.current_diffusion = equi_loc;
        end
        % For METIS internal resampling use with external data
        setappdata(0,'LCFS_nbpts',LCFS_nbpts);
    end
end


