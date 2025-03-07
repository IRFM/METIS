% METIS open metis simulator GUI
% Authors : J.F. Artaud, F. Imbeaux, J. Garcia, G. Giruzzi, 
%           T. Aniel, V. Basiuk, A. Bécoulet, C. Bourdelle, Y. Buravand, J. Decker, R. Dumont,
%           L.G. Eriksson, X. Garbet, R. Guirlet, G.T. Hoang, P. Huynh, E. Joffrin, X. Litaudon,
%           P. Maget, D. Moreau, R. Nouailletas, B. Pégourié, Y. Peysson, M. Schneider and J. Urban
% Copyright holder : Commissariat à l’Energie Atomique et aux Energies Alternatives (CEA), France
% CEA authorize the use of the METIS software under the CeCILL-C open source license https://cecill.info/licences/Licence_CeCILL-C_V1-en.html  
% The terms and conditions of the CeCILL-C license are deemed to be accepted upon downloading the software and/or exercising any of the rights granted under the CeCILL-C license.
%------------------------------------------
function metis

if isdeployed 
    dir
    if exist('Metis_splash_screen3.png')
        try
             splash('Metis_splash_screen3','png');
             drawnow
        catch
          f = errordlg(lasterr, 'splash error');
          waitfor(f);
        end
    elseif exist('Metis_splash_screen.png')
        try
             splash('Metis_splash_screen','png');
             drawnow
        catch
          f = errordlg(lasterr, 'splash error');
          waitfor(f);
        end
    else
         f = errordlg('unable to find image', 'splash error');
         waitfor(f);
    end
    warning off
    maxNumCompThreads('automatic');
    warning on
    
    if  ~verLessThan('matlab','R2014a') && ispc
	   opengl software
    end
end

% adaptation for MAC
if strcmp(computer,'MACI64')
    %------------------------------------------------------------------------
    % Retrieve screen dimensions (in correct units)
    %------------------------------------------------------------------------
    if exist('groot')
	set(groot,'units','pixels');
    else
	set(0,'units','pixels');    
    end
    %screensize = get(0,'ScreenSize');  % record screensize
end

% root du programme et initilisation du path
root = getappdata(0,'root');
if isempty(root)
        zineb_path;
        root = getappdata(0,'root');
end

% display data source for Lz
try
    display_lz_source
catch
    warning('Unable to display Lz data source');
end


% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('zeroda');
if ishandle(hform)
        zuiformvisible(hform);
	return
end

% version
[zver,zdate]        = zinebversion;

if ~isdeployed 
    warning off
    clear functions
    warning on
end

setappdata(0,'langue_cronos','anglais');
setappdata(0,'uicrossref',[]) ; % securite
rmappdata(0,'uicrossref') ;

if exist('import_ual')
	try 
		import_ual;
	end
end

if evalin('base','exist(''z0dinput'')')
	zcall0d('deja');
else
	zcall0d;
end

if ~isempty(getenv('METIS_INITIAl_DIR'))
    cd(getenv('METIS_INITIAl_DIR'))
end



