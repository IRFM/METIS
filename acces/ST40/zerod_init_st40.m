function z0dinput = zerod_init_st40(mode_exp,shot,gaz,temps,z0dinput)

% selection of ASTRA file
[f,p] = uigetfile('*.mat', 'Select a ASTRA simulation for ST40 (@ matfile format)');
if isempty(f) || isnumeric(f)
    z0dinput = [];
    return
end
z0dinput = zerod_init_st40_from_astra(fullfile(p,f),z0dinput);

    