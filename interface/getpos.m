%GETPOS  Takes figure positions from open Figure windows.

%        These positions are used by AUTOPLACE.

%

% Usage:

%  (1)   CD to the directory in which you want to store definitions.

%  (2)   Close all open Figure windows.

%  (3)   Open and position Figure windows for which

%        you would like position definitions.

%  (4)   Run 'getpos'.

%

%   See also: autoplace, figureCreateDelete

%



% Douglas L. Harriman

% Hewlett-Packard VCD

% doug_harriman@hp.com

% 9/30/98



function [] = getpos()



% Load existing definition file if possible

% Attempt to load the default positions from file

if exist('figures.ini','file'),

   load('figures.ini','-mat');

else,

   positionList = [];

end



% Get handles of open windows

figureHandles = sort(get(0,'Children'));

numFigures    = length(figureHandles);



% Initialize positions

figurePositions = zeros(numFigures,4);



% Get and store positions of figures

for i = 1:numFigures,

	figurePositions(i,:) = get(figureHandles(i),'Position');	   

end



% Make sure dimensions match

oldDim = size(positionList);

oldDim = oldDim(1);



if numFigures > oldDim,

   % Adding definitions

   % Create an empty matrix of the new size

   newPositionList = zeros(numFigures,4,numFigures);

   

   % Copy old data into the matrix

   for i = 1:oldDim,

      for j = 1:i,

      	newPositionList(j,:,i) = positionList(j,:,i);   

      end

   end

   

   % Copy new data into the matrix

   newPositionList(:,:,numFigures) = figurePositions;

	positionList = newPositionList;   

   

elseif numFigures == oldDim,

   % Just overwrite the last entry

   positionList(:,:,numFigures) = figurePositions;

   

elseif numFigures < oldDim,

   % Build the pad we need to slide the new definition in

   pad = zeros(oldDim-numFigures,4);

   

   % Replace the given entry

   positionList(:,:,numFigures) = [figurePositions; pad];

   

end



% Save the ini file

save('figures.ini','positionList');
