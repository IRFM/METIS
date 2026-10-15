function root = metis_octave_setup ()
% METIS_OCTAVE_SETUP  Set up the METIS path for GNU Octave.
%
%   root = metis_octave_setup ();
%
% Uses METIS's own zineb_path, so functions resolve in the same order as in
% MATLAB, then puts octave/compat first on the path. octave/compat holds small
% shims for MATLAB functions that Octave lacks (griddedInterpolant, pwch,
% chckxy) and restores Octave's built-in time(). It is never added to the
% MATLAB path.
%
% Compile the C mex files once beforehand with octave/build_mex.sh (optional
% but about 5x faster).

  if ! exist ("OCTAVE_VERSION", "builtin")
    error ("metis_octave_setup is for GNU Octave only; use metis_path in MATLAB");
  end

  root = fileparts (fileparts (mfilename ("fullpath")));
  addpath (root);
  zineb_path ();
  addpath (fullfile (root, "octave", "compat"), "-begin");

  % METIS's compiled PDE solver needs MPI/MUMPS; use the pure-M solver.
  setappdata (0, "MEXSOLVER_IN_METIS", 0);
end
