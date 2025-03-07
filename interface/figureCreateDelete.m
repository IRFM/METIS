% FIGURECREATEDELETE  Handles the positioning of figure windows
%                     so they do not overlap upon creation.
%
% Douglas L. Harriman
% Hewlett-Packard VCD
% doug_harriman@hp.com
% 9/30/98



function [] = figureCreateDelete(varargin)

% Parse on message
if nargin == 0
   msg = 'Create';
else   
   msg = varargin{1};
end
msg = lower(msg);
% Persistent variables
persistent positionList    % Definition of all figure positions
persistent figureHandles   % List of all user figure handles
persistent defaultPosition % Default figure position
% See if this is the first run for this session
if isempty(defaultPosition)
   % First run, initialize things
   figureHandles = [];
   % Get the root properties
   defaultPosition = get(0,'DefaultFigurePosition');
   % Attempt to load the default positions from file
   if exist('figures.ini','file')
      load('figures.ini','-mat');
   else
      %disp('No figure position definition file exists.');    
      return;
   end
   % Lock the mfile   
   mlock;
end



% Message defined behavior

switch msg
case 'create'

   

   % Get handle
   han = gcf;
   % Check to see if this is a user figure
   % Change this IF statement if you want to change the definition of 
   % a user figure window.
   if isempty(get(han,'Name')) & strcmp(get(han,'NumberTitle'),'on') ... 
         & prod(get(han,'position') == defaultPosition),

      % We have a user figure, log handle
      figureHandles = [figureHandles han];

      

      % Update close fcn.
      set(han,'DeleteFcn',['figureCreateDelete(''Delete'',' num2str(han) ');']);

      

      % Position the figures
      % See how many definitions we have
      numDefs = size(positionList);
      numDefs = numDefs(1);
      numFigures = length(figureHandles);
      if numFigures <= numDefs

         % Cycle through all figure handles placing figures
         for i = 1:numFigures

            % Set positions
           set(figureHandles(i),'Position',positionList(i,:,numFigures));

         end

      end

      

   end

   

   

case 'delete'
   % Store the handle

   if nargin > 1,
     % We have legal handle
      han = varargin{2};
		% Remove the handle from the list
		figureHandles = figureHandles(figureHandles ~= han);      
      % Position the figures
      % See how many definitions we have
      numDefs = size(positionList);
      numDefs = numDefs(1);
      numFigures = length(figureHandles);
      if numFigures <= numDefs,
         % Cycle through all figure handles placing figures
         for i = 1:numFigures,
            % Set positions
            set(figureHandles(i),'Position',positionList(i,:,numFigures));
         end

      end

      

   end

end



% EOF
