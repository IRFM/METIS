function provenance = make_provenance

% allocate provenance
ps = ids_gen('dataset_fair');
provenance = ps.ids_properties.provenance;
% create also provenance for version 4 of the dictionary
provenance.node{1}.path                   = which('metis'); 
provenance.node{1}.reference{1}.name      ='METIS'; 
provenance.node{1}.reference{1}.timestamp = datestr(now,30); 
[matlab_version,matlab_date] = version;
[s,t] = unix('which matlab');
provenance.node{end+1} = provenance.node{1};
provenance.node{end}.path                   = strtrim(t);
provenance.node{end}.reference{1}.name      = matlab_version;
provenance.node{end}.reference{1}.timestamp = matlab_date;
[s,t] = unix('uname -n -m -r -v -p -i -o');
if s == 0
    provenance.node{end+1} = provenance.node{1};
    provenance.node{end}.path                   = strtrim(t); 
    provenance.node{end}.reference{1}.name      = strtrim(strtok(t,'.')); 
    provenance.node{end}.reference{1}.timestamp = datestr(now,30); 
end
provenance.node{end+1} = provenance.node{1};
provenance.node{end}.path                   = pwd;
provenance.node{end}.reference{1}.name      = getenv('USER');
provenance.node{end}.reference{1}.timestamp = datestr(now,30);
[zver,zdate,git_hash] = zinebversion;
provenance.node{end+1} = provenance.node{1};
provenance.node{end}.path                   = git_hash;
provenance.node{end}.reference{1}.name      = zver;
provenance.node{end}.reference{1}.timestamp = zdate;
