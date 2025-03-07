% ZVERIFEXECUTION  checks that the simulation (in case of a result file) has run up to the end
%------------------------------------------------------------------------------- 
% fichier : zverifexecution.m  ->  zverifexecution
% 
% fonction Matlab 5 : 
%	checks that the simulation (in case of a result file) has run up to the end
%
% syntaxe :  
% 	zverifexecution
%  
% entrees :  
%  
% sorties :  
%  
% fonction �crite par F. Imbeaux
% version  4.0 du 05/11/2008 
%  
% liste des modifications :  
%  
%-------------------------------------------------------------------------------  
%  
function zverifexecution

% file type (Input or resultat)
filetype = evalin('base','param.gene.filetype');

if strcmp(filetype,'resultat')   % the test is done only for a result file
   kmax = evalin('base','param.gene.kmax');  % final time slice
   k    = evalin('base','param.gene.k');     % time slice where the calculation has stopped
   rebuilt = evalin('base','param.gene.rebuilt');     % if 1, there is an automatic rebuilt when the crash occurs
   time = evalin('base','data.gene.temps');  % Cronos time base
   if k ~= kmax
    if kmax > length(time)
    	h = warndlg('Warning : kmax is not consistent with the size of data.gene.temps','Warning');
    else
        if rebuilt   % if the rebuild has been done automatically already
	       h = warndlg(sprintf('The simulation has stopped at t =  %g s \nwhile it should have run until t = %g s', ...
	             time(k),time(kmax)), ...
		     'Simulation has stopped before its end');
        else
	       h = warndlg('The simulation has crashed, you should try to Recover Results','Simulation has crashed');
        end
    end
   else
      disp('The simulation has run successfully until its end')  % disp in the Matlab command window only
   end
end
