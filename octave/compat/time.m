function t = time ()
% Octave only: METIS's trait/time.m (a date-string helper) shadows Octave's
% built-in time(), which Octave's own clock() relies on, causing infinite
% recursion. Restore the built-in. Nothing in METIS calls time() as a function.
  t = builtin ('time');
end
