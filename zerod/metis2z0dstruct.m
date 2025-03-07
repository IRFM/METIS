% create z0dstruct structure from post structure
function [z0dstruct,z0dinput,zerod,profil0d] = z0post2evolution(post,temps)

z0dstruct.zs          = post.zerod;
z0dstruct.profil      = post.profil0d;
z0dinput              = post.z0dinput;
% suppression de la separatrice et Xdur
if isfield(z0dinput.exp0d,'Rsepa')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
end
if isfield(z0dinput.exp0d,'Zsepa')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
end

if isfield(z0dstruct.profil,'Rsepa')
	z0dstruct.profil = rmfield(z0dstruct.profil,'Rsepa');
end
if isfield(z0dstruct.profil,'Zsepa')
	z0dstruct.profil = rmfield(z0dstruct.profil,'Zsepa');
end
		
if isfield(z0dinput.exp0d,'XDURt')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURt');
end
if isfield(z0dinput.exp0d,'XDURx')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURx');
end
if isfield(z0dinput.exp0d,'XDURv')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURv');
end
% backward compatibilty (bug shortcut)
if isfield(z0dinput.exp0d,'flux')
    if ~isfield(z0dinput.exp0d,'edgeflux')
	  z0dinput.exp0d.edgeflux = z0dinput.exp0d.flux;
    elseif all(~isfinite(z0dinput.exp0d.edgeflux))
 	  z0dinput.exp0d.edgeflux = z0dinput.exp0d.flux;   
    end
    z0dinput.exp0d = rmfield(z0dinput.exp0d,'flux');
end
z0dstruct.z0dinput = z0dinput;
z0dstruct.info     = post.z0dinput.zsinfo;
zerod		= [];
profil0d	= [];
if (nargin > 1) && ~isempty(temps)
	indt = find(temps <= z0dstruct.zs.temps,1);
	if isempty(indt)
		indt = length(z0dstruct.zs.temps);
	end
	noms = fieldnames(z0dstruct.zs);
	for k = 1:length(noms)
		nomc = noms{k};
		if size(z0dstruct.zs.(nomc),1) == length(z0dstruct.zs.temps)
			zerod.(nomc) = z0dstruct.zs.(nomc)(indt);
		else
			zerod.(nomc) = z0dstruct.zs.(nomc);
		end
	end
	indt = find(temps >= z0dstruct.profil.temps,1);
	if isempty(indt)
		indt = length(z0dstruct.profil.temps)
	end
	noms = fieldnames(z0dstruct.profil);
	for k = 1:length(noms)
		nomc = noms{k};
		if size(z0dstruct.profil.(nomc),1) == length(z0dstruct.profil.temps)
			profil0d.(nomc) = z0dstruct.profil.(nomc)(indt,:);
		else
			profil0d.(nomc) = z0dstruct.profil.(nomc);
		end
	end

end 
