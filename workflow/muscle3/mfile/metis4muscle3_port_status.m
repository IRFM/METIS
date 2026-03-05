function status = metis4muscle3_port_status(instance,io)

% look if ports are connected
name = fieldnames(io.list);
for k=1:length(name)
    try
        status.(name{k}) = instance.is_connected(name{k});
    catch
        status.(name{k}) = false;
    end
end
