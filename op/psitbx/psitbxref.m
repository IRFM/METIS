%% The TCV Psi-Toolbox reference guide

% Insert_Contents_Here

%<HR>
%%% "psitbxdcd"
%<HR>

%%%% Purpose

% Creates a "psitbxdcd" object containing the description of the geometry of the chords of a diagnostic.

%%%% Syntax

% <PRE>dcd = psitbxdcd(rd,zd,phid,pvd,tvd,nd,t)</PRE>

%%%% Description

% The function "psitbxdcd" is used to create a "psitbxdcd" object that contains the description of the geometry of the chords of a diagnostic. The fields of a "psitbxdcd" object and the corresponding input arguments have the following meaning:

% <DL>
% <DT>"rd"<DD>The radial position in cylindrical coordinate of a reference point along the chord, very often the pin-hole
% <DT>"zd"<DD>The vertical position of the same point
% <DT>"phid"<DD>The toroidal angle of the same point
% <DT>"pvd"<DD>The poloidal viewing angle, measured between the horizontal plane and the chord, positive if the chord is pointing above the plane 
% <DT>"tvd"<DD>The toroidal viewing angle, measured between the vertical plane containing the machine axis and the chord reference point and the vertical plane containing the chord itself, positive if the chord is looking in the direction of increasing toroidal angle
% <DT>"nd"<DD>The number of points used to discretise the chord when building a "psitbxgrid" object of type "'Diagnostic-Chords'"
% <DT>"t"<DD>The time base for moving chords
% </DL>

% The size of "rd","zd","phid","pvd","tvd" follows the following rules:

% <UL>
% <LI>The size can retain the organisation of the diagnostic, for example 20&times;10 for 10 cameras with 20 detectors each, or 1028&times;1028 for a CCD.
% <LI>All these arguments must have the same size unless
% <UL><LI>It is the same for all chords and can be summarised in a scalar
%     <LI>It is empty and will be replaced by a default value</UL>
% <LI>For chords moving in time, the last dimension must correspond to the the length of "t".
% </UL>

%%%% Example

%%%%% Example 1

% Here is an example with two cameras with 20 detectors each, the first one with a poloidal coverage, the second one with a toroidal coverage:

rd = 1.2;
zd = 0;
phid = pi/8;
pvd = [linspace(-pi/3,pi/3,20);zeros(1,20)];
tvd = [zeros(1,20);linspace(-pi/3,pi/3,20)];
dcd = psitbxdcd(rd,zd,phid,pvd,tvd)
plot(dcd)
waitforbuttonpress

%%%%% Example 2: moving chord

% The possibility to have real-time scaning diagnostic is included. In that case the last dimension of the fields of the "psitbxdcd" object represents the time base. Given as an example <CITE>Le cauchemar du gyrotroniste</CITE>:

td = linspace(0,1,73);
rd = 1.2;
zd = 0;
phid = pi/2;
pvd = pi/6*cos(pi*td);
tvd = pi/6*sin(2*pi*td);
LeCauchemarDuGyrotroniste = psitbxdcd(rd,zd,phid,pvd,tvd,[],td);
plot(LeCauchemarDuGyrotroniste)
waitforbuttonpress

%<HR>
%%% "psitbxdcd/intersect"
%<HR>

%%%% Purpose

% Computes the intersection between diagnostic chords described by a "psitbxdcd" object and a cylindrical coordinate grid contained in a "psitbxgrid" object.

%%%% Syntax

% <PRE>s = intersect(dcd,g)</PRE>

%%%% Description

% The method "intersect" applies to a "psitbxdcd" object and computes the intersection between the chords and a grid in cylindrical coordinates. It returns an array of linear coordinates along the chords of size "[dcd.nd,size(dcd)]". The main purpose of intersect is to provide the coordinate of the points for building a "psitbxgrid" object of type "'Diagnostic-Chords'".
 
% For those chords who do not cross the domain of the grid, "NaN"'s are returned. For those chords who cross the domain twice, the first segment is considered.

%%%% Example

dcd = psitbxdcd('Demo');
s = intersect(dcd,psitbxgrid('Cylindrical','Grid','Default'));
gdc = psitbxgrid('Diagnostic-Chords','Points',{s},dcd)

%%%% Algorithm

% Only awfull calculation of intersections between lines, planes and cylinders using analytical geometry.

%<HR>
%%% "psitbxdcd/psitbxdxp"
%<HR>

%%%% Purpose

% Computes the intersection between diagnostic chords described by a "psitbxdcd" object and flux surfaces inferred fromm a flux distribution "psitbxpsi" object.

%%%% Syntax

% <PRE>s = psitbxdxp(dcd,psi,ps)</PRE>

%%%% Description

%%%% Example

dcd = psitbxdcd('Fan');
s = intersect(dcd,psitbxgrid('Cylindrical','Grid','Default'));
gdc = psitbxgrid('Diagnostic-Chords','Points',{s},dcd)

%%%% Algorithm

% Only awfull calculation of intersections between lines, planes and cylinders using analytical geometry.

%<HR>
%%% "psitbxfun"
%<HR>

%%%% Purpose

% Creates a "psitbxfun" object containing a function defined on a spatial grid.

%%%% Syntax

% <PRE>
% f = psitbxfun(x,g)
% f = psitbxfun(x,g,t)
% </PRE>

%%%% Description

% The function "psitbxfun" is used to create a "psitbxfun" object that contains a function defined on a spatial grid. If this function does not depend on time, "psitbxfun(x,g)" is used with "x" containing the function value(s) defined on the "psitbxgrid" g. If the function varies in time, use the form "psitbxfun(x,g,t)".

% The size of the matrix "x" must be "[size(g)[,n1,...][,nt]]". The optional dimension "n1,..." can be used for multi-value functions (for example the wave length for spectroscopic measurements). If the time argument is given, then the last dimension of "x", "nt" must correspond to "length(t)".
 
% For degenerated coordinates, that have a size of one, the following rules apply:

% <UL>
% <LI>if this coordinate is a scalar, the function is defined only on the corresponding subspace and extrapolation outside this subspace is not possible
% <LI>if this coordinate is a "nAN", the function is assumed to be independant of this spatial coordinate
% </UL>
% Singleton dimension in "x" arising from degenerated coordinates can be omitted.
 
%%%% Example

%%%%% Example 1

% A fairly simple function defined on a cylindrical coordinate grid:

rzphi = {0.62:.02:0.7,-.05:.02:.05,pi/180*(-3:2:3)};
g = psitbxgrid('Orthogonal','Grid',rzphi);
[rzphi{:}] = ndgrid(rzphi{:});
x = exp(-300*(rzphi{1} - g.x{1}(1)).^2) .* exp(-3000*rzphi{2}.^2);
f = psitbxfun(x,g)
plot(f)
waitforbuttonpress

%%%%% Example 2

% A very common concept in Tokamak physics is the profile, that is a quantity constant on a flux surface and therefore depending only on the flux coordinate. Such a profile is trivialy represented by a function on a 1-D grid, here the canonical parabolic profile:

g = psitbxgrid('Flux','Grid',{linspace(0,1,41),NaN,NaN});
f = psitbxfun(1-linspace(0,1,41)'.^2,g)
plot(f)
waitforbuttonpress

%<HR>
%%% "psitbxfun/cumsum"
%<HR>

%%%% Purpose

% Method to compute the indefinite integral of a "psitbxfun" object over spatial coordinates

%%%% Syntax

% <PRE>
% f = cumsum(f,dim)
% f = cumsum(f,dim,'metric')
% </PRE>

%%%% Description

% The method "cumsum" applies on a "psitbxfun" object defined on grid with "'Grid'" storage and returns a "psitbxfun" containing its indefinite integral  along the specified dimensions "dim". For example if "dim = 1", the integral is performed along the first coordinate, if "dim = [2,3]", the surface integral along coordinates 2 and 3 is returned, if "dim = [1,2,3]", the volume integral is returned.
 
% By specifying the option "'metric'", the appropriate metric will be taken into account; simply the function "f.*metric(f.grid,dim)" is integrated.

%%%% Algorithm

% Uses a trapeze integration.

%<HR>
%%% "psitbxfun/mean"
%<HR>

%%%% Purpose

% Method to average a "psitbxfun" object over spatial coordinates

%%%% Syntax

% <PRE>
% f = mean(f,dim)
% f = mean(f,dim,w)
% f = mean(f,dim,'metric')
% f = mean(f,dim,w,'metric')
% </PRE>

%%%% Description

% The method "mean" applies on a "psitbxfun" object defined on grid with "'Grid'" storage and returns a "psitbxfun" containing the function averaged along the specified dimensions "dim". For example if "dim = 1", the average is performed along the first coordinate, if "dim = [2,3]", the surface average along coordinates 2 and 3 is returned, if "dim = [1,2,3]", the volume average is returned.
 
% Note that the grid on which the returned function is defined will become degenerated in the dimension along which average is carried out.

% In the call "mean(f,dim,w)", "w" is a "psitbxfun" object playing the role of a weight in the average.

% By specifying the option "'metric'", the appropriate metric will be taken into account; simply the function "f.*metric(f.grid,dim)" is average.

%%%% Example

% Average plasma current density

jphi = psitbxtcv(10000,'JPHI');
jphiavg = mean(jphi,[1,2],'metric');
plot(jphiavg)

%%%% Algorithm

% <DL>
% <DT>"mean(f,dim,w)"<DD>"sum(f.*w,dim)"./sum(w,dim)"
% <DT>"mean(f,dim,w,'metric')"<DD>"sum(f.*w.*metric(f.grid,dim),dim)./sum(w.*metric(f.grid,dim),dim)"
% </DL>

%<HR>
%%% "psitbxfun/psitbxf2f"
%<HR>

%%%% Purpose

% Method to interpolate a "psitbxfun" object on a specified grid.

%%%% Syntax

% <PRE>
% f2 = psitbxf2f(f,g2)
% f2 = psitbxf2f(f,g2,der)
% f2 = psitbxf2f(f,der)
% </PRE>

%%%% Description

% The method "psitbxf2f" applies on a "psitbxfun" object to interpolate the embended function on new space points. In the call "psitbxf2f(f,g2)" the "psitbxfun" object represents a function defined on a "psitbxgrid" which must have "'Grid'" storage. "g2" is a "psitbxgrid" object specifying the new space points on which the function "f" must be interpolated. Those can be given in any coordinate systems. Points outside the definition domain of the function will be assigned a "NaN". If any coordinate of "g2" is a "NaN", which means that any function defined on "g2" is independant of this coordinate, then the function "f" must also be independant along the corresponding coordinate.

% The optional argument "der" is a three element vector that specifies the spatial derivative order to apply along each spatial coordinate. Note that since functions are interpolated with cubic spline, any order bigger than 3 will results in a zero function.

%%%% Example

%%%%% Example 1

% How to compute the normalised flux coordinate along diagnostic chords:

dcd = psitbxdcd('Demo');
psi = psitbxtcv(10000,.5,'01');
gdc = psitbxgrid('Diagnostic-Chords','Points',{intersect(dcd,psi.grid)},dcd);
rho = sqrt(psitbxf2f(psi,gdc));
plot(rho)

%%%%% Example 2

% How to compute the radial derivative of the poloidal flux which is related to the vertical magnetic field:

psi = psitbxtcv(10000);
dpsidz = psitbxf2f(psi,[1,0,0]);
plot(dpsidz)

%%%% Algorithm

% In a first step the points of "g2" where the function "f" must be interpolated are transformed to compute their coordinates in the coordinate system on which the function "f" is defined, using the method "psitbxg2g". Then the function "f" is interpolated in each spatial direction with cubic splines. The edge conditions are <CITE>Not-A-Knot</CITE> conditions except along angular coordinates which span 2&times;pi, in which case <CITE>periodic</CITE> edge conditions are used. Finaly this cubic spline interpolant is evaluated at the new positions, to the appropriate derivative order.

% Note that interpolation requires the function to be smooth. For example interpolating the square-root of the poloidal flux defined in cylindrical coordinates yields problem on the magnetic axis where its first derivative is not continuous.

%<HR>
%%% "psitbxfun/sum"
%<HR>

%%%% Purpose

% Method to integrate a "psitbxfun" object over spatial coordinates

%%%% Syntax

% <PRE>
% f = sum(f,dim)
% f = sum(f,dim,'metric')
% </PRE>

%%%% Description

% The method "sum" applies on a "psitbxfun" object defined on grid with "'Grid'" storage and returns a "psitbxfun" containing the function integrated along the specified dimensions "dim". For example if "dim = 1", the integral is performed along the first coordinate, if "dim = [2,3]", the surface integral along coordinates 2 and 3 is returned, if "dim = [1,2,3]", the volume integral is returned.
 
% Note that the grid on which the returned function is defined will become degenerated in the dimension along which integration is carried out.

% By specifying the option "'metric'", the appropriate metric will be taken into account; simply the function "f.*metric(f.grid,dim)" is integrated.

%%%% Example

% A non trivial way to compute the plasma current

jphi = psitbxtcv(10000,'JPHI');
ip = sum(jphi,[1,2]);
plot(ip)

%%%% Algorithm

% Uses a trapeze integration.

%<HR>
%%% "psitbxfun" Mathematical operators and functions
%<HR>

%%%% Purpose

% Methods for mathematical operations on "psitbxfun" objects

%%%% Syntax

% <PRE>
% f = f + a
% f = f + g
% f = f .* a
% f = sqrt(f)
% f = f > 0
% </PRE>

%%%% Description

% Most of the mathematical operators and functions have been overloaded for "psitbxfun" objects. For multi operand functions, the following rules apply:
 
% <UL>
% <LI>Some of them can be of class double, either scalar or with the appropriate size
% <LI>If all operands of class "psitbxfun" are not defined on the same grid, they will be interpolated on the definition grid of the first one using the "psitbxf2f" method.
% </UL>

%%%% Example

% Express the current density in MA/m<SUP>2</SUP>

jphi = psitbxtcv(10000,'JPHI')/1e6;

%<HR>
%%% "psitbxgrid"
%<HR>

%%%% Purpose

% Creates a "psitbxgrid" object containing the definition of a spatial grid

%%%% Syntax

% <PRE>
% g = psitbxgrid(type,storage,x[,p1,...])
% </PRE>

%%%% Description

% The function "psitbxgrid" is used to create a "psitbxgrid" object that contains the description of a spatial grid. The usual call is "psitbxgrid(type,storage,x)".

%%%%% "type"

% The character argument "type" specifies which coordinate system is used. Possible choices are:

% <DL>
% <DT>"'Orhogonal'"<DD>This is the usual orthogonal cartesian coordinate system:<DL><DT>"X"<DD>This axis lies on the midplane at toroidal position zero between sector #16 and #1, zero on machine axis<DT>"Y"<DD>This axis lies on the midplane at toroidal position pi/2, zero on machine axis<DT>"Z"<DD>This axis coincides with machine axis pointing upward, zero on the midplane<DD></DL>

% <DT>"'Cylindrical'"<DD>The cylindrical coordinate system whose axis coincides with the main axis of the machine:<DL><DT>"R"<DD>The radial coordinate<DT>"Z"<DD>The vertical coordinate, posivite upward, zero on the midplane<DT>"Phi"<DD>The azimutal coordinate, positive ccw seen from the top, zero between sector #16 and #1</DL>

% <DT>"'Toroidal'"<DD>The toroidal coordinate system whose main axis coincides with the main axis of the machine:<DL><DT>"A"<DD>The minor radius measured from the axis whose cylindrical coordinates are (R<SUB>0</SUB>,Z<SUB>0</SUB>)<DT>"Theta"<DD>The poloidal angle, zero on the HFS midplane, positive cw (This strange definition to get a right handed system; this also corresponds to the magnetic probe ordering)<DT>"Phi"<DD>The toroidal angle, positive ccw seen from the top, zero between sector #16 and #1</DL>

% <DT>"'Flux-coordinate'"<DL><DT>"Rho"<DD>The flux coordinate, usually the square-root of the poloidal flux normalised such that it is zero on magnetic axis and one on the last close flux surface. But can also be for example the toroidal flux<DT>"Theta"<DD>The poloidal angle, zero on the HFS midplane, positive cw (strange definition to get a right handed system; this also corresponds to the magnetic probe ordering)<DT>"Phi"<DD>The toroidal angle, positive ccw seen from the top, zero between sector #16 and #1</DL>

% <DT>"'Diagnostic-Chords'"<DD>This one dimensional coordinate system allows to scan a chord of a diagnostic:<DL><DT>"S"<DD>The linear distance along the chord measured from the chord reference point (see "psitbxdcd")</DL>
% </DL>

%%%%% "storage"

% The character argument "storage" specifies the way the grid points are stored in the "x" field of the "psitbxgrid" object. This field and argument is a cell array containing three double arrays. Possible values for "storage" are:

% <DL>
% <DT>"'Grid'"<DD>In this case the grid is made of the points of a mesh specified by the three vectors of "x". Note that this storage method can not be used with "Diagnostic-Chords"<DT>"'Points'"<DD>This storage method is used for points which does not align on a grid. In that case all the arrays of "x" must have the same size and the grid is made of the points whose coordinates by each triplet. Note that the shape of these arrays is free and can reflect the organisation of the problem<DT>"'Time-Points'"<DD>This form is used for points which move in time. All the arrays of "x" must also be of the same size but in addition their last dimension represents time samples
% </DL>

% For 1-D or 2-D grids, the corresponding array in "x" can be replaced either by a constant, for example for things that lie entirely on a plane, or by a single "NaN" which means that the problem does not depend on this particular coordinate. A grid can also be extended to other dimension that the basic space and time coordinates. For example a fourth dimension can be added that would represents a frequency or a wave length if necessary.
 
% There are some combinations of "type" and "storage" for which the "x" argument can be given as the string "'Default'", which will returns a predefined grid, in particular "psitbxgrid('Flux-Coordinate','Grid','Default')"
 
%%%%% Optional arguments

% For its points to be absolutely positioned in space, some coordinate systems requires additional information. Such an information is given in a list of optional arguments to the function "psitbxgrid". Here is a list of these pieces of information:

% <DL>
% <DT>"t"<DD>The time base for "'Time-Points'" storage. It can be found in: <UL><LI>A double array<LI>The "t" field of a "psitbxgrid" object<LI>The "t" field of a "psitbxfun" object<LI>The "t" field of a "psitbxpsi" object<LI>The "t" field of a "psitbxdcd" object</UL>

% <DT>"r0,z0"<DD>The position of the axis for a "'Toroidal'" coordinate system or of the magnetic axis for a "'Flux-Coordinate'" system. Note that this axis can move in time, in which case its size must be coherent with the time base. It can be found in: <UL><LI>A cell array made of two double array<LI>The "r0" and "z0" fields of a structure<LI>The "r0" and "z0" fields of the "par" field of a "psitbxgrid" object<LI>The "rmag" and "zmag" fields of a "psitbxpsi" object<LI>The "r0" and "z0" fields of the "grid" field of a "psitbxfun" object<LI>The "r0" and "z0" fields of the "grid" field of a "psitbxpsi" object</UL>

% <DT>"psi"<DD>A "psitbxpsi" object containing the distribution in cylindrical coordinate of the normalised poloidal flux. This information can be optionaly given for a "'Flux-Coordinate'" grid. It can be found in: <UL><LI>A "psitbxpsi" object of format "'01'"<LI>In the "psi" field of a structure</UL>

% <DT>"fsd"<DD>A "psitbxpsi" object containing the description of the constant poloidal flux contour in toroidal coordinates. This information can be optionaly given for a "'Flux-Coordinate'" grid. It can be found in: <UL><LI>A "psitbxpsi" object of format "'FS'"<LI>In the "fsd" field of a structure</UL>

% <DT>"dcd"<DD>The description of diagnostic chords. This information is required for grid of type "'Diagnostic-Chords'" and can be found in: <UL><LI>A "psitbxdcd" object<LI>In the "dcd" field of a structure</UL>
% </DL>

%%%% Example

%%%%% Example 1

% A "'Toroidal'" grid with "'Grid'" storage

x = {linspace(0,.2,6),linspace(-pi,pi,7),linspace(0,2*pi,25)};
r0z0 = {1,0};
g = psitbxgrid('Toroidal','Grid',x,r0z0)
plot(g)
waitforbuttonpress

%%%%% Example 2

% Two probes moving along the Z direction

t = linspace(0,.1,11);
x = {repmat([.6;.7],1,11),[t;t]-.75,pi/3};
g = psitbxgrid('Cylindrical','Time-Points',x,t)
plot(g)
waitforbuttonpress

%%%% Listing

type psitbxgrid

%<HR>
%%% "psitbxgrid/metric"
%<HR>

%%%% Purpose

% Method to compute the metric of a "psitbxgrid" object

%%%% Syntax

% <PRE>
% f = metric(g,dim)
% </PRE>

%%%% Description

% The method "metric" applies on a "psitbxgrid" object with "'Grid'" storage and returns a "psitbxfun" containing the metric along the specified dimensions "dim". For example if "dim = 1", the metric along the first coordinate is returned, if "dim = [2,3]", the metric of the surface element along coordinates 2 and 3 is returned, if "dim = [1,2,3]", the volume element is returned. In addition if "dim" is negative, the norm of the gradient of the specified coordinate is returned.

%%%% Example

% The very important gradient of rho entering the expression of the SEF (Shape Enhancement Factor)

fsd = psitbxtcv(10000,'FS');
grho = metric(fsd,-1);
plot(grho)

%%%% Algorithm

% For the flux-coordinate system, the metric matrix is not diagonal; for the others it is trivial. Here are the formulae:

% <IMG SRC=PSITBXG2G-2.GIF>

%<HR>
%%% "psitbxgrid/psitbxg2g"
%<HR>

%%%% Purpose

% Method to transform a "psitbxgrid" object on a different coordinate system

%%%% Syntax

% <PRE>
% g2 = psitbxg2g(g,type)
% g2 = psitbxg2g(g,type,p1,...)
% </PRE>

%%%% Description

% The method "psitbxg2g" applies on a "psitbxgrid" object to transform the coordinates of its points in a new coordinate system. In the call "psitbxg2g(g,type)", "g" is a "psitbxgrid" object and "type" a character string that specifies the new coordinate system; it can take the value of the "type" argument of the "psitbxgrid" function. "psitbxg2g" returns a "psitbxgrid" object with the same spatial points but with coordinates expressed in the new system. 

% Note that the "storage" may be altered; in particular "'Grid'" storage may be transformed in "'Points'" or "'Time-Points'" storage because points in the new system do not align on a grid anymore.

% It is impossible to map points on a "'Diagnostic-Chords'" system, since this system is discrete in space.

% Some transformations require additional information that will be sought in the "par" field of the "g" argument or in additional optional arguments, following the same rules as for the "psitbxgrid" function. These transformations are:
 
% <DL>
% <DT>Cylindrical to toroidal and toroidal to cylindrical<DD>Requires "r0" and "z0"

% <DT>Cylindrical to Flux coordinate<DD>Requires "r0", "z0" and a "psitbxpsi" of form "'01'"

% <DT>Flux coordinate to toroidal<DD>Requires a "psitbxpsi" of form "'FS'"</DL>

%%%% Example

% Position of the chord points in toroidal coordinates:

dcd = psitbxdcd('Demo');
s = intersect(dcd,psitbxgrid('Cylindrical','Grid','Default'));
gdc = psitbxgrid('Diagnostic-Chords','Points',{s},dcd);
g = psitbxg2g(gdc,'Toroidal',{.89,0});
plot(g)

%%%% Algorithm

% There are a set of basic transformations defined. All the others are composed from the basic set, which is:
% <IMG SRC=PSITBXG2G-1.GIF>

%<HR>
%%% "psitbxpsi"
%<HR>

%%%% Purpose

% Creates a "psitbxpsi" object containing the poloidal flux description or the flux surface description.

%%%% Syntax

% <PRE>psi = psitbxpsi(x,grid,t,form)</PRE>
% <PRE>psi = psitbxpsi(x,grid,t,form,[rmag,zmag,psimag])</PRE>
% <PRE>psi = psitbxpsi(x,grid,t,form,[rmag,zmag,psimag],iter,tol)</PRE>

%%%% Description

% The function "psitbxpsi" is used to create a "psitbxpsi" object, which is a child of the "psitbxfun" class that contains either the poloidal flux in cylindrical coordinates or the flux surface contours in toroidal coordinates. In the call "psitbxpsi(x,grid,t,form)" depending on the "form" argument:

% <DL>
% <DT>"'FS'"<DD>"x" is the minor radius in toroidal geometry of the surfaces of constant poloidal flux. "x" must be defined on the "psitbxgrid" object "grid" of type "'Flux-Coordinate'" with storage "'Grid'". "t" is the time base.

% <DT>"'10','+0','-0','01'"<DD>"x" is the poloidal flux defined on the "psitbxgrid" object "grid" of type "'Cylindrical'" with storage "'Grid'". The first character of "form" represents the value or the sign of the poloidal flux on the magnetic axis, the second on the last close flux surface. "'+0'" is the LIUQE representation, "'01'" the normalised flux used for flux coordinate.
% </DL>

% Optional arguments allows to specify the position of the magnetic axis "rmag" and "zmag" and the axis flux "psimag" and convergence parameters for the method "psitbxp2p", "iter" and "tol".
 
%<HR>
%%% "psitbxpsi/psitbxfsd" (private method)
%<HR>

%%%% Purpose

% Find the magnetic surface radius from "psitbxpsi" object

%%%% Syntax

% <PRE>
% a = psitbxfsd(psi)
% </PRE>

%%%% Description

% The function "psitbxfsd" is used to find the magnetic surface radius and return a "psitbxfun" on a "'Flux-coordinate'" grid which describes a(&rho;,&theta;).

%%%% Algorithm

% The coefficients for the cubic spline base function combination are computed using the matrices stored with the grid of the "psitbxpsi" object as c<SUB>kl</SUB> = &sum;<SUB>ij</SUB> M<SUP>R</SUP><SUB>ki</SUB> &psi;<SUB>ij</SUB> M<SUP>Z</SUP><SUB>jl</SUB><SUP>t</SUP>. These coefficients are then used to estimate &psi; at any point with &psi;(R,Z) = &sum;<SUB>kl</SUB> c<SUB>kl</SUB>b<SUB>k</SUB>(R)b<SUB>l</SUB>(Z), the gradient of the &psi; distribution &nabla;&psi;(R,Z) = &sum;<SUB>kl</SUB> c<SUB>kl</SUB> [b<SUB>k</SUB>'(R)b<SUB>l</SUB>(Z),b<SUB>k</SUB>(R)b<SUB>l</SUB>'(Z)] and its Laplacian &Delta;&psi;(R,Z) = &sum;<SUB>kl</SUB> c<SUB>kl</SUB> [b<SUB>k</SUB>''(R)b<SUB>l</SUB>(Z),b<SUB>k</SUB>'(R)b<SUB>l</SUB>'(Z);b<SUB>k</SUB>'(R)b<SUB>l</SUB>'(Z),b<SUB>k</SUB>(R)b<SUB>l</SUB>''(Z)].

% The derivatives of &psi; along the radius are computed using d&psi;(R(a,&theta;),Z(a,&theta;))/da = &psi;'(a) = -&part;&psi;/&part;R cos&theta; + &part;&psi;/&part;Z sin&theta; and d<SUP>2</SUP>&psi;(R(a,&theta;),Z(a,&theta;))/da<SUP>2</SUP>/2 = &psi;''(a) = &part;<SUP>2</SUP>&psi;/&part;R<SUP>2</SUP> cos<SUP>2</SUP>&theta; - 2&part;<SUP>2</SUP>&psi;/&part;R&part;Z sin&theta; cos&theta; + &part;<SUP>2</SUP>&psi;/&part;Z<SUP>2</SUP> sin<SUP>2</SUP>&theta;

% The value of &psi; along a radius is approximated by &psi;(a+da) = &psi;(a) + &psi;'(a) da + &psi;''(a)/2 da<SUP>2</SUP>. Then the radius of the flux surface with label &psi;<SUB>0</SUB> is found by solving the equation &psi;(a+da) = &psi;<SUB>0</SUB> whose solution is da = (&plusmn;(&psi;'<SUP>2</SUP>+2&psi;''(&psi;<SUB>0</SUB>-&psi;))<SUP>1/2</SUP> - &psi;')/&psi;''. Then a is iteratively replaces by a &rarr; a + da until the update drops below a threshold given by "psi.tol". Normaly the function &psi; in the format "'01'" increases along a radius, except in the private region of the divertor. To treat this case one of the two solutions for da is selected according to the sign of &psi;': da = (sign(&psi;')(&psi;'<SUP>2</SUP>+2&psi;''(&psi;<SUB>0</SUB>-&psi;))<SUP>1/2</SUP> - &psi;')/&psi;''.

%<HR>
%%% "psitbxpsi/psitbxfsg"
%<HR>

%%%% Purpose

% Method to compute geometrical quantities from a "psitbxpsi" object

%%%% Syntax

% <PRE>
% fsg = psitbxfsg(fsd)
% </PRE>

%%%% Description

% The method "psitbxfsg" applies on a "psitbxpsi" object with format "'FS'" and returns a structure with fields of class "psitbxfun":

% <DL><DT>"grho"<DD>flux surface average of the gradient of rho<DT>"surf"<DD>flux surface area<DT>"area"<DD>flux surface cross-section<DT>"vol"<DD>volume inside a flux surface<DT>"darea"<DD>cross-section between two flux surfaces divided by the flux coordinate increment<DT>"dvol"<DD>volume between two flux surfaces divided by the flux coordinate increment
% </DL>

%%%% Example

fsd = psitbxtcv(10000,'FS');
fsg = psitbxfsg(fsd);
plot(fsg.grho)

%%%% Algorithm

% Uses the following definitions:

% <DL><DT>"grho"<DD>"mean(metric(g,-1),[2,3])"<DT>"surf"<DD> "sum(metric(g,[2,3]),[2,3])"<DT>"area"<DD>"cumsum(sum(metric(g,[1,2]),2),1)" <DT>"vol"<DD>"cumsum(sum(metric(g,[1,2,3]),[2,3]),1)"<DT>"darea"<DD>"sum(metric(g,[1,2]),2)"<DT>"dvol"<DD>"sum(metric(g,[1,2,3]),[2,3])"
% </DL>

%<HR>
%%% "psitbxpsi/psitbxmag" (private method)
%<HR>

%%%% Purpose

% Find the magnetic axis and normalise a "psitbxpsi" object

%%%% Syntax

% <PRE>
% p = psitbxmag(p)
% </PRE>

%%%% Description

% The function "psitbxmag" is used to find the magnetic axis for a given flux distribution and then to normalise the flux and return a "psitbxpsi" object with format "'01'".

%%%% Algorithm

% The coefficients for the cubic spline base function combination are computed using the matrices stored with the grid of the "psitbxpsi" object as c<SUB>kl</SUB> = &sum;<SUB>ij</SUB> M<SUP>R</SUP><SUB>ki</SUB> &psi;<SUB>ij</SUB> M<SUP>Z</SUP><SUB>jl</SUB><SUP>t</SUP>. These coefficients are then used to estimate &psi; at any point with &psi;(R,Z) = &sum;<SUB>kl</SUB> c<SUB>kl</SUB>b<SUB>k</SUB>(R)b<SUB>l</SUB>(Z), the gradient of the &psi; distribution &nabla;&psi;(R,Z) = &sum;<SUB>kl</SUB> c<SUB>kl</SUB> [b<SUB>k</SUB>'(R)b<SUB>l</SUB>(Z),b<SUB>k</SUB>(R)b<SUB>l</SUB>'(Z)] and its Laplacian &Delta;&psi;(R,Z) = &sum;<SUB>kl</SUB> c<SUB>kl</SUB> [b<SUB>k</SUB>''(R)b<SUB>l</SUB>(Z),b<SUB>k</SUB>'(R)b<SUB>l</SUB>'(Z);b<SUB>k</SUB>'(R)b<SUB>l</SUB>'(Z),b<SUB>k</SUB>(R)b<SUB>l</SUB>''(Z)].

% If the gradient at (R+&Delta;R,Z+&Delta;Z) is approximated with &nabla;&psi;(R+&Delta;R,Z+&Delta;Z) = &nabla;&psi;(R,Z) + &Delta;&psi;(R,Z)(&Delta;R,&Delta;Z), the extremum is found where the condition &nabla;&psi = 0 is satisfied, that is (&Delta;R,&Delta;Z) = -&Delta;&psi;(R,Z)<SUP>-1</SUP>&nabla;&psi;(R,Z). Thus (R,Z) is iteratively replaced by (R,Z) &rarr; (R,Z) - &Delta;&psi;(R,Z)<SUP>-1</SUP>&nabla;&psi;(R,Z) until the update falls bellow a critical threshold specified by "p.tol".

%<HR>
%%% "psitbxpsi/psitbxp2p"
%<HR>

%%%% Purpose

% Method to change the format of a "psitbxpsi" object

%%%% Syntax

% <PRE>
% psi = psitbxp2p(psi,'01')
% fsd = psitbxp2p(psi,'FS')
% fsd = psitbxp2p(psi,'FS',g)
% </PRE>

%%%% Description

% The method "psitbxp2p" applies on a "psitbxpsi" object and returns a "psitbxpsi" object with the specified form:

% <DL><DT>"'01'"<DD>The poloidal flux contained in the "psi" argument wil be searched for its maximum, normalised and the location of the maximum stored in the "rmag" and "zmag" fields
 
% <DT>"'FS'"<DD>The contours of constant normalised poloidal flux will be searched. A third argument allows to specify the flux-coordinate grid on which these contours are to be sought.</DL>

%%%% Example

psi = psitbxtcv(10000,.5);
fsd = psitbxp2p(psi,'FS');
plot(fsd)

%%%% Algorithm

% Uses Gauss-Newton iterations to locate either the maximum or the contour level. Iteration convergence is controlled by the "iter" and "tol" fields of the "psitbxpsi" object.

%<HR>
%%% "psitbxtcv"
%<HR>

%%%% Purpose

% Easy access to quantities relevant for the Psi-Toolbox directly from the TCV shot files

%%%% Syntax

% <PRE>psitbxtcv(shot,form)</PRE>
% <PRE>psitbxtcv(shot,time,form)</PRE>

%%%% Description

% The function "psitbxtcv" builds object compatible with the Psi-Toolbox from the TCV shot files. Depending on the "form" argument it returns:

% <DL>
% <DT>"'+0'"<DD>A "psitbxpsi" object of form "'+0'" containing the poloidal flux from the equilibrium reconstruction by LIUQE
% <DT>"'01'"<DD>A "psitbxpsi" object of form "'01'" containing the poloidal flux from the equilibrium reconstruction by LIUQE normalised using the "psitbxp2p" method
% <DT>"'FS'"<DD>A "psitbxpsi" object of form "'FS'" containing the flux surface description as precalculated by the Psi-Toolbox and stored in the shot file
% <DT>"'JPHI'"<DD>A "psitbxfun" object containing the toroidal current density from the equilibrium reconstruction by LIUQE
% </DL>

% An optional time argument can be specified. Only the times for which a sample in the shot file at most 1e-6 second away exists are returned.
 
%%%% Example

psi = psitbxtcv(10000,'+0')
plot(psi)
waitforbuttonpress

%<HR>
%%% "psitbxtcv"
%<HR>

%%% Cubic spline interpolation

% In some operations, specificaly the interpolation of a "psitbxfun" on a new set of points or grid as well as the calculation of the flux surfaces or of the intersection of diagnostic chords with the flux surfaces, a function specified by values on a given grid must be interpolated on new points. This is done by fitting N-dimensional cubic splines on the data points.

%%%% Cubic spline base functions

% For a given set of knots {t<SUB>1</SUB>,...,t<SUB>N</SUB>} a set of base functions {b<SUB>1</SUB>,...b<SUB>N+2</SUB>} can be derived which satisfy the following conditions:

% <DL>
% <DD>b<SUB>i</SUB>(x) are piecewise cubic polynomials
% <DD>b<SUB>i</SUB>(x) = 0 for x &le; t<SUB>i-3</SUB> or x &ge; t<SUB>i+1</SUB>
% <DD>b<SUB>i</SUB>(x), b<SUB>i</SUB>'(x) and b<SUB>i</SUB>''(x) are continuous
% </DL>

t = (0:5)';
x = linspace(0,5,51);
b = bspbase(-3:8,4,x);
plot(x,b), set(gca,'xtick',t)
waitforbuttonpress

%%%% Cubic spline combination

% A function of x can be approximated by a linear combination of these base functions with appropriate coefficients {c<SUB>1</SUB>,...,c<SUB>N+2</SUB>} as f(x) = &sum;<SUB>i</SUB> b<SUB>i</SUB>(x)c<SUB>i</SUB>. For a 2-D function, the coefficients form a matrix: f(x,y) = &sum;<SUB>ij</SUB> c<SUB>ij</SUB>b<SUB>i</SUB>(x)b<SUB>j</SUB>(y).

%%%% Edge conditions

% Since there are 2 more base functions than knots, the exact interpolation of a function given by its values on the knots by a linear combination of cubic spline base functions requires two additional constraints. These are specified as two edge conditions with the following choice:

% <DL>
% <DT>not-a-knot<DD>f''' is continuous at t<SUB>2</SUB> and t<SUB>N-1</SUB>
% <DT>natural<DD>f''(t<SUB>1</SUB>) = 0 and f''(t<SUB>N</SUB>) = 0
% <DT>symetric<DD>f'(t<SUB>1</SUB>) = 0 and f'''(t<SUB>1</SUB>) = 0
% <DT>periodic<DD>f'(t<SUB>1</SUB>) = f'(t<SUB>N</SUB>) and f''(t<SUB>1</SUB>) = f''(t<SUB>N</SUB>)
% </DL>

%%%% Fast calculation of the coefficients

% So for a function given by its values on the knots f(t<SUB>j</SUB>) = f<SUB>j</SUB> and a selected edge condition pair, a set of coefficients {c<SUB>1</SUB>,...,c<SUB>N+2</SUB>} can be calcultated such that f<SUB>j</SUB> = &sum;<SUB>i</SUB> c<SUB>i</SUB>b<SUB>i</SUB>(t<SUB>j</SUB>). Inversely, the calculation of the coefficients can be casted in a matrix such that c<SUB>i</SUB> = &sum;<SUB>j</SUB> M<SUB>ij</SUB>f<SUB>j</SUB>. Then the interpolation of the function f(x) is easily performed by f(x) = &sum;<SUB>ij</SUB> b<SUB>i</SUB>(x)M<SUB>ij</SUB>f<SUB>j</SUB>. In addition the derivatives of b<SUB>i</SUB>(x) are very easy to compute, so are the derivatives of f, for example f'(x) = &sum;<SUB>ij</SUB> b<SUB>i</SUB>'(x)M<SUB>ij</SUB>f<SUB>j</SUB>. For a 2-D function defined on a grid (x<SUB>k</SUB>,y<SUB>l</SUB>) f(x<SUB>k</SUB>,y<SUB>l</SUB>) = f<SUB>kl</SUB>, the calculation of the coefficients requires two matrices c<SUB>ij</SUB> = &sum;<SUB>kl</SUB> M<SUB>ik</SUB>f<SUB>kl</SUB>M<SUB>lj</SUB><SUP>t</SUP>

M = csdec(t,t);
f = t.^3;
plot(t,f,'or',x,bspsum(-3:8,M*f,x),'-r',x,bspsum(-3:8,M*f,x,2),'--r'), set(gca,'xtick',t)
waitforbuttonpress
