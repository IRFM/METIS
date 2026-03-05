function mstruct = get_matlab_serialized_structure(data)

tmpf = sprintf('%s.mat',tempname);
fid = fopen(tmpf,'w');
fwrite(fid,uint8(data),'uint8');
fclose(fid);
mstruct = load(tmpf,'internal_state');
mstruct = mstruct.internal_state;
delete(tmpf);
