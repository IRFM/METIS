function msg = make_checkpoint_message(internal_state)

% rm python fields that can't be saved in Matlab and will be recreate when
% workflow resume
internal_state.actor    = rmfield(internal_state.actor,'instance');
internal_state.actor    = rmfield(internal_state.actor,'logger');
internal_state.actor.io = rmfield(internal_state.actor.io,'ports');
% loop on raw data
if isfield(internal_state.actor.io,'data')
    noms = fieldnames(internal_state.actor.io.data);
    for k=1:length(noms)
        internal_state.actor.io.data.(noms{k}) = rmfield(internal_state.actor.io.data.(noms{k}),'data_raw');
    end
end

tmpf = sprintf('%s.mat',tempname);
save(tmpf,'internal_state');
fid = fopen(tmpf,'r');
data = fread(fid,Inf,'uint8');
fclose(fid);
delete(tmpf);
% make a Python bytes stream
var = make_python_bytes(data);
%var = py.bytes(uint8(data));
msg = py.libmuscle.Message(internal_state.workflow.clock.t_cur,internal_state.workflow.clock.t_next,var);

% for testing
%save snapshot_message

