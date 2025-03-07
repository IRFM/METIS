% METIS open metis simulator GUI
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
    else
         f = errordlg('unable to find image', 'splash error');
         waitfor(f);
    end
    warning off
    maxNumCompThreads('automatic');
    warning on
end

% root du programme et initilisation du path
root = getappdata(0,'root');
if isempty(root)
        zineb_path;
        root = getappdata(0,'root');
end

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('zeroda');
if ishandle(hform)
        zuiformvisible(hform);
	return
end

% version
[zver,zdate]        = zinebversion;

clear functions

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



