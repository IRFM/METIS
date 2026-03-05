function metis4muscle3_set_external_data(data,settings,state,sawtooth)

%  loop on fieldnames
noms = fieldnames(state);
equi_done = false;
equi_ext  = false;
curdif_ext = false;
for k=1:length(noms)
    if isfield(data,noms{k}) && ~isempty(data.(noms{k}))
        
        switch noms{k}
            case 'equilibrium'
                % separate cases
                if ~equi_done
                    create_external_equilibrium_from_imas([],0,data.equilibrium);
                    equi_done = true;
                end
               equi_ext = true;
            case 'current_diffusion'
                % separate cases
                if ~equi_done
                    create_external_equilibrium_from_imas([],0,data.equilibrium);
                    equi_done = true;
                end
                curdif_ext = true;
            otherwise
                setappdata(0,settings.(noms{k}).ext_name,data.(noms{k}));
        end
    end
end

% remove the non activated external data
if ~equi_ext && isappdata(0,'EQUILIBRIUM_EXP')
    rmappdata(0,'EQUILIBRIUM_EXP');
end
if ~curdif_ext && isappdata(0,'CURDIFF_EXP')
    rmappdata(0,'CURDIFF_EXP');
end

% to be sure but not used in this workflow
setappdata(0,'METIS_EXTERNAL_CURDIF_EQUI_TRANSITION_FACTOR',1);

% for mechanism mimicking sawtooth
switch sawtooth.on_off
    case 'on'
        setappdata(0,'LAST_TIME_ST',sawtooth.LAST_ST_TIME);
        setappdata(0,'INDICE_INV_ST',sawtooth.INDICE_INV_ST);
        setappdata(0,'MINIMAL_TAU_ST',sawtooth.MINIMAL_TAU_ST);
    otherwise
        % clear external data for mechanism mimicking ST
        if isappdata(0,'LAST_TIME_ST')
            rmappdata(0,'LAST_TIME_ST');
        end
        if isappdata(0,'INDICE_INV_ST')
            rmappdata(0,'INDICE_INV_ST');
        end
        if isappdata(0,'MINIMAL_TAU_ST')
            rmappdata(0,'MINIMAL_TAU_ST');
        end        
end



