%  ZNE  courte description
%------------------------------------------------------------------------------- 
% fichier :  zne.m  ->  zne 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   [filer,liner] = zne(file,line) 
%  
% entrees :  
%  file = 
%  line = 
%  
% sorties :  
%  filer = 
%  liner = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  1.9  du  11/03/2002  
%  
% #auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function [filer,liner] = zne(file,line)

mess = lasterr;
if nargin < 1
	file = '';
end
if nargin < 2
	line =[];
end
mode  =0;

if ~isempty(file)
   % cas help
   ind = findstr(file,'file:');
   if ~isempty(ind)
      file = file((ind+5):end);
   end
end


if isempty(file) & ~isempty(mess)
	tok = 'in ==>';
	indf = findstr(mess,tok);
	if ~isempty(indf)
		messf = mess((indf+length(tok)):end);
		indnl = find(messf==sprintf('\n')); 
		if ~isempty(indnl)
			file  = messf(1:(indnl-1));
		end
	end
	mode = 1;
else
	filem = which(file);
	if ~isempty(filem)
		if isempty(findstr(filem,'built-in'))
			file =filem;
		end
	end
end

if isempty(line) & ~isempty(mess) & (mode == 1)
	tok = 'On line ';
	indf = findstr(mess,tok);
	if ~isempty(indf)
		messf = mess((indf+length(tok)):end);
		indend = findstr(messf,'==>'); 
		if ~isempty(indend)
			line  = str2num(messf(1:(indend-1)));
		end
	end
	
end


if strcmp(computer,'GLNX86')
	if ~isempty(file)
		if ~isempty(line)
			[s,t]= unix(sprintf('unset LD_LIBRARY_PATH; %s -line %d %s &', ...
                                     getappdata(0,'editeur'),line,file));
		else
			[s,t]= unix(sprintf('unset LD_LIBRARY_PATH; %s  %s &', ...
                                     getappdata(0,'editeur'),file));
		end
	end
else
	if ~isempty(file)
		if ~isempty(line)
			[s,t]= unix(sprintf('%s -line %d %s &', ...
                                     getappdata(0,'editeur'),line,file));
		else
			[s,t]= unix(sprintf('%s  %s &', ...
                                     getappdata(0,'editeur'),file));
		end
	end
end


if nargout > 0
	filer = file; 
	liner = line;
end
