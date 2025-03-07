function post = metis_func(z0dinput, fast_run)

if nargin < 2
    fast_run = false;
end

evalin('base','clear  z0dstruct');

% overwrite parameters whith reference parameters if defined
z0dinput.option = z0doverwriteparam(z0dinput.option);

% sort du mode evolution
z0dinput.option.evolution = 0;

% recopie machine dans option
zerod_machine_name;

if fast_run
    [post.zerod,void,post.profil0d] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
else
    [post.zerod,void,post.profil0d] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
end

post.z0dinput = z0dinput;

