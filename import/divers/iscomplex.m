function a = iscomplex(x)

% ISCOMPLEX	test whether a number is complex
% function a = iscomplex(x)
%
% test whether x is complex

if ~isnumeric(x)
  % non-numeric input (cell, struct, char): not complex. Needed under GNU Octave,
  % whose own library functions call iscomplex on cell arrays.
  a = false;
  return
end
a = any(any(imag(x) ~= 0));
