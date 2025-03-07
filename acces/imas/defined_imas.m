%% --------------------------------------------------------------
%% FUNCTIONS RETURNS 1 IF VARIABLE IS INITIALIZED (FLOAT/INTEGER) 
%% OR FILLED (STRING), RETURN 0 OTHERWISE 
%% --------------------------------------------------------------

function [answer] = defined_imas(var)

if nargin < 1

  disp(['-----------------------------------------------------------------'])
  disp(sprintf('To be used with 1 argument:'))
  disp(sprintf('undefine_imas(variable)'))
  disp(sprintf('=> Returns 1 if variable defined, 0 otherwise'))
  disp(['-----------------------------------------------------------------'])
  answer = [];
   
else
  
  %% STRING VARIABLE
  if isstr(var) 
    answer = ~isempty(var);
    
  %% FLOAT OR INTEGER VARIABLE
  else
    
    if length(var)<1
      answer = 0;
    else
    
      if var==-999999999
	answer = 0;
      else
	answer = 1;
      end    
    end
  end

end

