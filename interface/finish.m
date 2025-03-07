%  FINISH  courte description  
%------------------------------------------------------------------------------- 
% fichier :  finish.m  ->  finish 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   finish 
%  
% entrees :  
%  
%  
% sorties :  
%  
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  18/11/2005  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function finish

% save z0dinput for reuse 
filename = sprintf('METIS_z0dinput@%s',datestr(now,'yyyy-mm-dd-THH-MM-SS'));
filename2 = sprintf('METIS_data_backup@%s',datestr(now,'yyyy-mm-dd-THH-MM-SS'));
if evalin('base','exist(''z0dinput'')')
   if  evalin('base','exist(''post'')')
       z0dinput_new =  evalin('base','z0dinput');
       z0dinput_new = rmfield(z0dinput_new,'exp0d');
       z0dinput_new = rmfield(z0dinput_new,'mode_exp');
       if isfield(z0dinput_new,'exp')
	  z0dinput_new = rmfield(z0dinput_new,'exp');
       end
       if isfield(z0dinput_new.option,'reference_parameters') && isempty(z0dinput_new.option.reference_parameters)
 	  z0dinput_new = rmfield(z0dinput_new.option,'reference_parameters');      
       end

       z0dinput_post =  evalin('base','post.z0dinput');
       z0dinput_post = rmfield(z0dinput_post,'exp0d');
       z0dinput_post = rmfield(z0dinput_post,'mode_exp');
       if isfield(z0dinput_new,'exp')
	  z0dinput_post = rmfield(z0dinput_post,'exp');
       end
       if isfield(z0dinput_post.option,'reference_parameters') && isempty(z0dinput_post.option.reference_parameters)
 	  z0dinput_post = rmfield(z0dinput_post.option,'reference_parameters');      
       end
       
       warning off
       rep=zcompstruct(z0dinput_post,z0dinput_new,1e-3,0,'',0,0);
       warning on
       if rep ~= 0
                evalin('base',sprintf('save(''%s'',''z0dinput'',''post'')',filename));           
                fprintf('METIS input (z0dinput) have been changed but not used for a simulation: z0dinput have been saved in %s\n',filename);
       elseif ~isappdata(0,'METIS_FILENAME') ||  ~isappdata(0,'METIS_INTERFACE_TITLE')
                evalin('base',sprintf('save(''%s'',''z0dinput'',''post'')',filename2));           
                fprintf('METIS data have been changed but not saved : data have been saved in %s\n',filename2);     
       end
   else
       evalin('base',sprintf('save(''%s'',''z0dinput'')',filename));
       fprintf('METIS input (z0dinput) have been changed but not used for a simulation: z0dinput have been saved in %s\n',filename);
   end
end
