% script pour le calcul de metis avec zerodevolution
if ~exist('z0dstruct','var') && exist('post','var') && isfield(post,'profil0d')
    % create z0dstruct from metis data if possible
    [z0dstruct,z0dinput_void,zerod,profil0d] = metis2z0dstruct(post);
end
% securite
if isfield(z0dinput.exp0d,'Rsepa') && isempty(z0dinput.exp0d.Rsepa)
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
end
if isfield(z0dinput.exp0d,'XDURx') && isempty(z0dinput.exp0d.XDURx)
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURx');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURt');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURv');
end
% initialisation
% choix de l'index
z0plotsc;
drawnow
tmin = min(z0dstruct.zs.temps);
index_restart = length(z0dstruct.zs.temps) + 1;
if index_restart > length(z0dinput.cons.temps)
	tmax =  max(z0dinput.cons.temps);
else
	tmax = z0dinput.cons.temps(index_restart);
end
prompt={sprintf('Time for restart (%g <= t <= %g s)',tmin,tmax)};
name='METIS evolution restart';
numlines=1;
defaultanswer={sprintf('%g',z0dstruct.zs.temps(end))};
answer=inputdlg(prompt,name,numlines,defaultanswer);
if ~isempty(answer)
	trestart = min(tmax,max(tmin,str2num(answer{1})));
	index_restart = find(z0dinput.cons.temps >= trestart,1)
end
if index_restart > length(z0dinput.cons.temps);
	error('METIS simulation is already finished');
end
hdlg = z0dpatience('evolution');
% boucle sur les temps :
for k_index = index_restart:length(z0dinput.cons.temps)
        z0dpatience(k_index/(length(z0dinput.cons.temps)+1));
	cons1t  = zerod_get1t(z0dinput.cons,k_index);
	geo1t   = zerod_get1t(z0dinput.geo,k_index);
	exp0d1t = zerod_get1t(z0dinput.exp0d,k_index);
	[zs,profil,z0dstruct] = zerodevolution(z0dstruct,z0dinput.option,z0dinput.cons.temps(k_index),cons1t,geo1t,z0dinput.exp0d,exp0d1t,exp0d1t);
end
z0dpatience(1);

post.z0dinput = z0dstruct.z0dinput;
post.zerod    = z0dstruct.zs;
post.profil0d = z0dstruct.profil;
delete(hdlg)
