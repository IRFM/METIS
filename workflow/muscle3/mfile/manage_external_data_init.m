function data_out = manage_external_data_init(data_1t,clock)

% make 3 time sclices data structure to use always the same method for
% resamping
vt       = ones(3,1);
data_out = [];
%dt       = abs(clock.t_next - clock.t_cur);
% loop on field names of input data structure
noms = fieldnames(data_1t);
for k=1:length(noms)
    switch noms{k}
        case 'x'
            data_out.x = data_1t.x;
        case 'temps'
            data_out.temps = data_1t.temps + cat(1,-abs(clock.dt),0,abs(clock.dt));
        otherwise
            data_out.(noms{k}) = vt * data_1t.(noms{k})(:)';
    end
end
