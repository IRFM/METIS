%AUTOPLACE  Automatic placement of figure windows
%           so that they do not all appear in the
%           same position.
%
% AUTOPLACE     by itself toggles the autoplacement state.
% AUTOPLACE ON  turns the autoplacment on.
% AUTOPLACE OFF turns the autoplacment off.
%
%  See also: getpos, figureCreateDelete
%
% Douglas L. Harriman
% Hewlett-Packard VCD
% doug_harriman@hp.com
% 9/30/98



function [] = autoplace(msg)



% If no message, assume empty string

if nargin == 0, msg = ''; end



% Switch on message

switch lower(msg),

case 'on'

   % Turn it on by setting the default fcns.

   set(0,'DefaultFigureCreateFcn','figureCreateDelete(''Create'');');

	set(0,'DefaultFigureDeleteFcn','figureCreateDelete(''Delete'');');

   disp('Figure autoplacement activated');drawnow

   

case 'off'

   % Turn it on by setting the default fcns.

   set(0,'DefaultFigureCreateFcn','');

	set(0,'DefaultFigureDeleteFcn','');

   disp('Figure autoplacement deactivated');drawnow

   

case ''

   % Check current default fcn

   fcn = get(0,'DefaultFigureCreateFcn');

   if strcmp(fcn,'figureCreateDelete(''Create'');')

      % It's on, so turn it off

      autoplace off;

   else,

      % It's off, so turn it on

      autoplace on;

   end

end

