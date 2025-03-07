% script de gestion du nom de la machine pour METIS
if ~isfield(z0dinput.option,'machine')
    z0dinput.option.machine = z0dinput.machine;
elseif isempty(z0dinput.option.machine)
    z0dinput.option.machine = z0dinput.machine;
else
    z0dinput.machine = z0dinput.option.machine;
end
if ~isfield(z0dinput.option,'shot')
    z0dinput.option.shot = z0dinput.shot;
elseif isempty(z0dinput.option.shot)
    z0dinput.option.shot = z0dinput.shot;
else
    z0dinput.shot = z0dinput.option.shot;
end
