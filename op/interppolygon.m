function varargout = interppolygon(X,N,method)
% INTERPPOLYGON  Interpolates a polygon.
%
%   Y = INTERPPOLYGON(X,N) interporpolates a polygon whose coordinates are
%   given in X over N points, and returns it in Y.
%
%   The polygon X can be of any dimension (2D, 3D, ...), and coordinates
%   are exected to be along columns. For instance, if X is 10x2 matrix (2D
%   polygon of 10 points), then Y = INTERPPOLYGON(X,50) will return a
%   50x2 matrix.
%
%   INTERPPOLYGON uses interp1 for interpolation. The interpolation method
%   can be specified as 3rd argument, as in Y=INTERPPOLYGON(X,N,'method').
%   Allowed methods are 'nearest', 'linear', 'spline' (the default),
%   'pchip', 'cubic' (see interp1).
%
%   ALGORITHM
%   The point to point distance array is used to build a metric, which will
%   be interpolated by interp1 over the given number of points.
%   Interpolated points are thus equally spaced only in the linear case,
%   not in other cases where interpolated points do not lie on the initial
%   polygon. If this is an issue, try caalling INTERPPOLYGON twice, as in:
%   >> Y = INTERPPOLYGON(X,50); % Make a spline interpolation%   
%   >> Z = INTERPPOLYGON(Y,50); % Will correct for uneven space between points
%
%   OUTPUT
%   On top of the interpolated polygon, [Y M] = INTERPPOLYGON(X,N) will
%   return in M the metric of the original polygon.
%
%   EXAMPLE
%   X = rand(5,2);
%   Y = interppolygon(X,50);
%   figure, hold on
%   plot(X(:,1),X(:,2),'ko-')
%   plot(Y(:,1),Y(:,2),'r.-')

% ------------------ INFO ------------------
%   Authors: Jean-Yves Tinevez 
%   Work address: Max-Plank Institute for Cell Biology and Genetics, 
%   Dresden,  Germany.
%   Email: tinevez AT mpi-cbg DOT de
%   January 2009
%   Permission is given to distribute and modify this file as long as this
%   notice remains in it. Permission is also given to write to the author
%   for any suggestion, comment, modification or usage.
% ------------------ BEGIN CODE ------------------

    if nargin < 3
        method = 'spline';
    end
    if nargin < 2 || N < 2
        N = 2;
    end

    % 1. Build metric
    [pl bl]         = polygonlength(X);
    orig_metric     = [   0 ;  cumsum(bl/pl) ];
    
    % 2. Interpolate
    interp_metric   = ( 0 : 1/(N-1) : 1)';
    Y               = interp1( ...
        orig_metric,...
        X,...
        interp_metric,...
        method);
    
    % 3. Ouputs
    varargout{1} = Y;
    if nargout > 1
        varargout{2} = orig_metric;
    end
    
    %% Subfunction
    % Calculate the point to point distanceof each polygon point
    function varargout = polygonlength(X)
        
        n_dim       = size(X,2);
        delta_X     = 0;
        for dim = 1 : n_dim
            delta_X = delta_X + ...
                diff(X(:,dim)).^2 ;
        end
        
        branch_lengths  = sqrt( delta_X );
        pl              = sum( branch_lengths );
        
        varargout{1}    = pl;
        if nargout > 1
            varargout{2} = branch_lengths;
        end
        
    end

end


%  Copyright (c) 2009, Jean-Yves Tinevez
%  All rights reserved.
%  
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%  
%      * Redistributions of source code must retain the above copyright 
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright 
%        notice, this list of conditions and the following disclaimer in 
%        the documentation and/or other materials provided with the distribution
%        
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
