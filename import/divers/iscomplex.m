function a = iscomplex(x)

% ISCOMPLEX	test whether a number is complex
% function a = iscomplex(x)
%
% test whether x is complex

a = any(any(imag(x) ~= 0));
