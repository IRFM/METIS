% importation d'une separtrice
function cr = zimport_sepa(fichier)

% cr 
cr = 0;

% test des entrees
if nargin <1 
	fichier ='';
end

if ishandle(fichier)
    hfc     = fichier;
    fichier = '';
    set(hfc,'visible','off');
else
    hfc = [];
end


% selection du fichier
if isempty(fichier)
	[file,path]=uigetfile('*.*','Nom du fichier a charger ?');
	if ~ischar(file)
		% annulation
                if ishandle(hfc)
		     set(hfc,'visible','on');
		end
		return
	end
	fichier = fullfile(path,file);
end

% chargement
try
   sepas = load(fichier);
catch
   disp('erreur de lecture :')
   disp(lasterr)
   cr = -1;
   if ishandle(hfc)
	 set(hfc,'visible','on');
   end
   return
end

if isempty(sepas)
	disp('pas de donnees ...')
	cr = 1;
        if ishandle(hfc)
	    set(hfc,'visible','on');
        end
	return
end

% matlab ou ascii
if isstruct(sepas) 
	% c'est un fichier matlab
	noms =fieldnames(sepas)	;
	if length(noms) == 1
		% tout dans la meme variable
		RZ  = getfield(sepas,noms{1});
		si  = size(RZ);
		if si(1) == 2
			R = RZ(1,:);
			Z = RZ(2,:);
			T = 0;
		elseif si(2) == 2
			R = RZ(:,1)';
			Z = RZ(:,2)';
			T = 0;
		elseif si(1) == 3
			T = RZ(1,:);
			R = RZ(2,:);
			Z = RZ(3,:);
		elseif si(2) == 3
			T = RZ(:,1)';
			R = RZ(:,2)';
			Z = RZ(:,3)';
		else
			disp('Erreur : la matrice doit avoir une dimension 2 * N ou N * 2')
			cr = 2;
                        if ishandle(hfc)
	                   set(hfc,'visible','on');
                        end
			return
		end
	elseif length(noms) == 2
		% 2 variables seulements
                R   = getfield(sepas,noms{1});
                % si plus d'une dimension, la premiere est le temps
                if size(R,2) == 1
      	           R = R';
                end
                Z   = getfield(sepas,noms{2});
                if size(Z,2) == 1
      	           Z = Z';
                end
	else 
		% recherche des noms utiles
		% par ordre conventionnel
		if strmatch('R',noms,'exact')
			R = sepas.R;
		elseif strmatch('R',noms)
			R   = getfield(sepas,noms{min(strmatch('R',noms))});
		elseif strmatch('r',noms,'exact')
			R = sepas.r;
		elseif strmatch('r',noms)
			R   = getfield(sepas,noms{min(strmatch('r',noms))});
		elseif strmatch('X',noms,'exact')
			R = sepas.X;
		elseif strmatch('X',noms)
			R   = getfield(sepas,noms{min(strmatch('X',noms))});
		elseif strmatch('x',noms,'exact')
			R = sepas.x;
		elseif strmatch('x',noms)
			R   = getfield(sepas,noms{min(strmatch('x',noms))});
		end
                % si plus d'une dimension, la premiere est le temps
                if size(R,2) == 1
      	           R = R';
                end

		if strmatch('Z',noms,'exact')
			Z = sepas.Z;
		elseif strmatch('Z',noms)
			Z   = getfield(sepas,noms{min(strmatch('Z',noms))});
		elseif strmatch('z',noms,'exact')
			Z = sepas.z;
		elseif strmatch('z',noms)
			Z   = getfield(sepas,noms{min(strmatch('z',noms))});
		elseif strmatch('X',noms,'exact')
			Z = sepas.Y;
		elseif strmatch('Y',noms)
			Z   = getfield(sepas,noms{min(strmatch('Y',noms))});
		elseif strmatch('y',noms,'exact')
			Z = sepas.y;
		elseif strmatch('y',noms)
			Z   = getfield(sepas,noms{min(strmatch('y',noms))});
		end
                % si plus d'une dimension, la premiere est le temps
                if size(Z,2) == 1
      	           Z = Z';
                end
      
		if strmatch('temps',noms,'exact')
			T =sepas.temps;
		elseif strmatch('TEMPS',noms,'exact')
			T =sepas.TEMPS;
		elseif strmatch('time',noms,'exact')
			T =sepas.time;
		elseif strmatch('times',noms,'exact')
			T =sepas.times;
		elseif strmatch('TIME',noms,'exact')
			T =sepas.TIME;
		elseif strmatch('TIMES',noms,'exact')
			T =sepas.TIMES;
		elseif strmatch('T',noms)
			T =  getfield(sepas,noms{min(strmatch('T',noms))});
		elseif strmatch('t',noms)
			T =  getfield(sepas,noms{min(strmatch('t',noms))});
		else
			T = 0;
		end
                % si plus d'une dimension, la premiere est le temps
                if size(T,1) == 1
      	            T = T';
                end
   
       end
else
	% c'est un fichier ascii
	% tout dans la meme variable
	RZ  = sepas;
	si  = size(RZ);
	if si(1) == 2
		R = RZ(1,:);
		Z = RZ(2,:);
		T = 0;
	elseif si(2) == 2
		R = RZ(:,1)';
		Z = RZ(:,2)';
		T = 0;
	elseif si(1) == 3
		T = RZ(1,:);
		R = RZ(2,:);
		Z = RZ(3,:);
	elseif si(2) == 3
		T = RZ(:,1)';
		R = RZ(:,2)';
		Z = RZ(:,3)';
	else
		disp('Erreur : la matrice doit avoir une dimension 2 * N ou N * 2')
		cr = 2;
                if ishandle(hfc)
	             set(hfc,'visible','on');
                end
		return
	end	
end

% test si pas vide
if isempty(R) | isempty(Z)
	disp('Pas de donnee ....');
	cr = 4;
        if ishandle(hfc)
	     set(hfc,'visible','on');
        end
	return
end
if any(size(R)~=size(Z))
	disp('dimension incompatible');
	cr = 5;
        if ishandle(hfc)
	    set(hfc,'visible','on');
        end
	return
end

% type
R = double(R);
Z = double(Z);

% swap de R et Z si R change de signe et z non
if any(R <0) & all(Z >0)
	S = Z;
	Z = R;
	R = S;
end

% mise en forme de la separatrice
for k = 1:size(R,1)
	% ordre 0 -> 2*pi
	ra  = 0.5 .* (min(R(k,:)) + max(R(k,:)));
	za  = 0.5 .* (min(Z(k,:)) + max(Z(k,:)));
	a = angle((R(k,:)-ra) + i .* (Z(k,:)-za)); 
	[a,ind] = sort(a);
	R(k,:)  = R(k,ind);% manque ra,za
	Z(k,:)  = Z(k,ind);% manque ra,za
	% angle distinct
	ind     = find(diff(a) ==0);
	if ~isempty(ind)
		for l =1:length(ind)
			indlp = ind(l) + 1;
			if indlp > length(a)
				indlp = 1;
			end
			indlm = ind(l) - 1;
			if indlm < 1
				indlm = length(a);
			end
			R(k,ind(l)) = 0.5 .* (R(k,indlp) + R(k,indlm));
			Z(k,ind(l)) = 0.5 .* (Z(k,indlp) + Z(k,indlm));
		end
	end
end

% recuperation des donnees du workspace
t      = evalin('base','data.gene.temps');
nbmax  = evalin('base','param.gene.nbsepa');
vt     = ones(size(t));
r0      = evalin('base','data.geo.r0');
b0      = evalin('base','data.geo.b0');

% rechantillonage de la separatrice 
if length(T) > 1
	RR      = tsample(R,T,t,'fen');
	ZZ      = tsample(Z,T,t,'fen');
   % on repporte au extremite la dreniere valeur definie	
	ind      = find(t<=min(T));
	if ~isempty(ind)
		RR(ind,:) = ones(length(ind),1) * R(1,:); 
		ZZ(ind,:) = ones(length(ind),1) * Z(1,:); 
	end
	ind      = find(t>=max(T));
	if ~isempty(ind)
		RR(ind,:) = ones(length(ind),1) * R(end,:); 
		ZZ(ind,:) = ones(length(ind),1) * Z(end,:); 
	end
	R = RR;
	Z = ZZ;
else

   R = ones(size(t)) * R;
   Z = ones(size(t)) * Z;

end





if (size(R,2) > 35) & (size(R,2) < (nbmax -1)) &  ...
	((all(R(:,end) == R(:,1)) & all(Z(:,end) == Z(:,1))) | ...
	(all(R(:,end) ~= R(:,1)) & all(Z(:,end) ~= Z(:,1))))
	% pas de reechantillonage spatiale
	% fermeture
	if (all(R(:,end) ~= R(:,1)) & all(Z(:,end) ~= Z(:,1)))
		% la separtrice doit se refermer
		R(:,end+1) = R(:,1);
		Z(:,end+1) = Z(:,1);
	end
else
	% rechantillonage spatiale
   vee   = ones(1,size(R,2));
	teta  = vt * linspace(0,2*pi,nbmax - 1);
	rmin  = min(R,[],2) * vee;
	rmax  = max(R,[],2) * vee;
	ra    = 0.5 .* (rmin + rmax);   
	zmin  = min(Z,[],2) * vee;
	zmax  = max(Z,[],2) * vee;
	za    = (zmin + zmax) ./ 2;
	alpha = angle( (R-ra) + i .* (Z-za)) + pi;
   ind   = find(alpha(:,1) > pi);
   if ~isempty(ind)
      alpha(ind,:) = alpha(ind,:) - pi;
   end
	R     = tsplinet(alpha,R,teta);
	Z     = tsplinet(alpha,Z,teta);
	% la separtrice doit se refermer
	R(:,end) = R(:,1);
	Z(:,end) = Z(:,1);
end


% extraction des infos de la geometrie
ve    = ones(1,size(R,2));
rmin  = min(R,[],2);
rmax  = max(R,[],2);
ra    = 0.5 .* (rmin + rmax);   
a     = 0.5 .* (rmax - rmin);
zmin  = min(Z,[],2);
zmax  = max(Z,[],2);
za    = (zmin + zmax) ./ 2;
b     = 0.5 .* (zmax - zmin);
k     = b ./ a;
mask1 = (Z == (max(Z,[],2)*ve));
mask2 = (Z == (min(Z,[],2)*ve));
   
rzmax = max(R .* mask1,[],2);
rzmin = max(R .* mask2,[],2);
cl    = ra - rzmin;
cu    = ra - rzmax;
d     = (cl+cu) ./2 ./ a;
    
% calcul des parametres pour helena
hr0     = ra;
hz0     = za;
ha      = a;
he1     = k;
htrl    = - asin(cl ./ a); 
htrh    =   asin(cu ./ a); 

% pour b0

prompt={'Valeur de B0 dans le vide:','Rayon de mesure de B0 :'};
def={sprintf('%g',mean(b0 .* r0 ./ ra)),sprintf('%g',mean(ra))};
dlgTitle='Donnees pour le champ toroidal';
lineNo=1;
answer=zinputdlg(prompt,dlgTitle,lineNo,def);

if isempty(answer)
	b0 = b0 .* r0 ./ ra;
else
	b0 = str2num(answer{1});
	r0 = str2num(answer{2});
	b0 = b0 .* r0 ./ ra;
end

% assignation des donnees dans le workspace
zassignin('base','data.geo.r0',hr0);  
zassignin('base','data.geo.z0',hz0),  
zassignin('base','data.geo.a',ha);   
zassignin('base','data.geo.e1',he1);  
zassignin('base','data.geo.trh1',htrh);
zassignin('base','data.geo.trb1',htrl);
zassignin('base','data.geo.ind1',0 .* vt);
zassignin('base','data.geo.b0',b0);
zassignin('base','data.geo.mode',2 .* vt);
zassignin('base','data.geo.R',single(R));   
zassignin('base','data.geo.Z',single(Z));

if ishandle(hfc)
	 set(hfc,'visible','on');
         zuisavenonok
end

