%% The TCV Psi-Toolbox tutorial

% Insert_Contents_Here

%%% Introduction

%%%% History and future

% <STRONG>The Psi-Toolbox provides tools for manipulating space and time dependant quantities in the variable geometry of TCV.</STRONG> Targeted tasks are the analysis of the diagnostic data, analysis and simulation for transport and heating,... In summer 1995, extensive discussions with the potential users allowed to define the various wishes. Then in collaboration with Mathias Anton, we proposed a set of basic functionalities. In 1996, Christian Deschenaux coded a first version of these basic routines, which I later published as version 1. The result was pretty well described by Richard's words: <CITE>I have spent the whole morning battling against the Psi-Toolbox.</CITE> The fact is that the functionality was certainly adequat but the lack of detailed documentation made it unusable.
 
% The advent of Matlab 5 encouraged me to completely revisite the Psi-Toolbox. Having written the first version, I realised that the various functionalities could be unified in a couple of simple concepts and the new capabilities of Matlab 5 would make the use of these concepts very clear. So users, if any, of the first Psi-Toolbox should forget their old coding habits and only retain what they want to achieve before starting reading this tutorial.
 
% Although as you will see the concepts are very simple, the number of possibilities and combinations they generate is hughe. So please let me know if you feel at some point that this tutorial is not clear enough. I will try to find better or more explicite examples. Also be tolerant with the numerous bugs you certainly will discover. 

% <P ALIGN=CENTER><IMG SRC=richard.jpg></P>

% <P ALIGN=CENTER>Staring Richard A. Pitts<BR>Thanks to Basil for giffing the figures</P>

%%%% Classes, objects and methods

% Version 2.0 of the TCV Psi-Toolbox intensively uses the capability of Matlab 5 for object oriented programming. Thus to better understand how to use this toolbox and how to get profit from it, you need to acquire a few concepts.

% A <DFN>class</DFN> is a new data type reserved for a specific purpose. An <DFN>ojbect</DFN> is an instance of a class. The fisrt main step in the definition of a class is the choice of its data structure. In Matlab, objects are stored in structures whose fields contain the various information specific to the class. However these fields are neither visible nor accessible for the user. Instead, together with a class, <DFN>methods</DFN> are created which in principle safely manipulate these fields to perform operations specific to the class. One of these methods, the most fundamental, is the <DFN>class constructor</DFN> which creates the corresponding objects. In Matlab, the methods for a given class are gathered in a directory whose name is the class name preceeded by the string "$8A" and the constructor is the method with the same name as the class.

% Some basic methods, such as "size" or "ndims", may apply to different classes. In such cases Matlab calls the appropriate method on the basis of the class of its arguments. Such methods are called <DFN>overloaded</DFN>.

% A class can also be the <DFN>child</DFN> or a <DFN>parent</DFN> class. In that case, the child class inherits the fields and methods of its parent, to which specific fields and/or methods are added. It means that if a method is not defined for the child but only for the parent, the parent method will be used.

% This is almost all what you need to know to use the Psi-Toolbox but if necessary further details can be found in chapter 14 of <CITE>Using Matlab</CITE>.

%%% Coordinate systems

%%%% Coordinate systems: the psitbxgrid class

% The first basic tools to take out from the Psi-Toolbox are the different coordinate systems. The information defining the systems you will use are stored in objects of class "psitbxgrid". There are different coordinate systems and we will start with the "'Cylindrical'" coordinates (r,z,phi) whose axis corresponds to the main vertical axis of the tokamak. Suppose Richard, an edge physicist, wants to study in details what happens around an inner wall tile. It is certainly a good idea for him to define a spatial grid in cylindrical coordinates around his tile, let say from 0.62 to 0.7 m in radial direction, from -0.05 to 0.05 m in vertical direction and from -3 to 3 degrees in toroidal angle. Thus he needs to construct an object of class "psitbxgrid".

rzphi = {0.62:.02:0.7,-.05:.02:.05,pi/180*(-3:2:3)};
GridCyl = psitbxgrid('Cylindrical','Grid',rzphi)
plot(GridCyl)
drawnow

%%%% How to manipulate psitbxgrid objects: some methods for that class

% You see that Matlab displays for Richard the fields of the "psitbxgrid" object and the corresponding contents. To use these fields, access with the usual dot syntax is defined. However to protect the integrity of the object, Richard is not allowed to directly modify them. The different methods applying on a "psitbxgrid" object will do that safely.

% For example the field "x" is a cell array, each cell being the vector of the mesh points in each spatial direction. I will discuss the other fields later.

% The function "size" has been redefined to give the effective size of the grid defined by the "psitbxgrid" object. It is a method of the "psitbxgrid" class.

GridCyl.label
GridCyl.x
GridCyl.x{1},GridCyl.x{2},GridCyl.x{3}
size(GridCyl)
which size -all
help psitbxgrid/size

%%%% When things do not align on a grid: the storage field

% Richard has also installed on his tile two sets of ten probes of some kind. Of course these do not lie on a regular grid as the previously defined abstract points. So their position is specified by a triplet of coordinates, named "'Points'". It is the meaning of the "storage" field of a "psitbxgrid" object. In such a case the organisation of each cell can be arbitrary, for example here 2 by 10 matrices to represent the two sets of two probes. In addition, the numbers come from the drawing workshop and are specified on a cartesian "'Orthogonal'" system, whose z-axis coincides with the machine main axis and whose x-z plane is at phi = 0.

xyz = {repmat(.625,2,10),[zeros(1,10);linspace(-.01,.01,10)],[linspace(-.01,.01,10);zeros(1,10)]};
ProbePointsOrt = psitbxgrid('Orthogonal','Points',xyz)
plot(GridCyl), hold on, plot(ProbePointsOrt,'ob'), hold off
drawnow

%%%% When things move: the Time-Points storage

% Richard will tell you that things are never so simple. In fact his probes move in time. So the position of Richard's probe should be stored in multi-dimensional arrays whose last index is the time index. In this example each three spatial coordinate should be a 2 by 10 by number of times array. Then Richard must build a "psitbxgrid" object with "'Time-Points'" storage:

vx = .1;
t = [.1,.3,.5];
dim = [ones(1,ndims(xyz{1})),length(t)];
xyzt = {repmat(xyz{1},dim)+vx*repmat(reshape(t,dim),[size(xyz{1}),1]),repmat(xyz{2},dim),repmat(xyz{3},dim)};
MovingProbePointsOrt = psitbxgrid('Orthogonal','Time-Points',xyzt,t)
plot(GridCyl), hold on, plot(MovingProbePointsOrt,'ob'), hold off
drawnow

%%%% Coordinate transformation: the psitbxg2g method

% Of course the first thing Richard wants to know is: where are my probes in cylindrical coordinates. The arithmetic is to be found in text book. The receipe to be used is the "psitbxg2g" method, which will hopefully do the job for him.

ProbePointsCyl = psitbxg2g(ProbePointsOrt,'Cylindrical')
ProbePointsCyl.x{1},ProbePointsCyl.x{2},ProbePointsCyl.x{3}

%%%% One more step toward the tokamak: the toroidal coordinate system

% Toroidal coordinates is quite natural for a toroidal device. To transform orthogonal cartesian or cylindrical coordinates to toroidal coordinates, Richard must specify the position of the todoidal axis. This is achieved with an additional argument to the "psitbxgrid" constructor. This argument can take several form. The one shown here is a cell array with the (r,z) position of the axis. This information is then stored in the "par" field of the "psitbxgrid" object.

athetaphi = {linspace(0,.2,6),linspace(-pi,pi,7),linspace(0,2*pi,25)};
r0z0 = {1,0};
GridTor = psitbxgrid('Toroidal','Grid',athetaphi,r0z0)
GridTor.par
plot(GridTor)
drawnow

%%%% Parameter in coordinate transformation

% The same kind of additional parameter may be necessary for transformation to toroidal coordinate. In this example a "psitbxgrid" object with toroidal coordinates is used, which contains in its "par" field the necessary information:

MovingProbePointsTor = psitbxg2g(MovingProbePointsOrt,'Toroidal',GridTor)

%%%% One more step toward variable configuration: moving magnetic axis

% In TCV configuration varies in time, in particular the magnetic axis moves. In the Psi-Toolbox, this motion can be taken into account with a time vector for the radial and vertical position of the axis. In this case the coordinate system changes at each time step. Note for ECH people that relativistic effect are not included.
 
r0z0 = {ones(1,length(t)),t*10-1};
MovingGridTor = psitbxgrid('Toroidal','Grid',athetaphi,r0z0,t)
plot(MovingGridTor)
drawnow

%%%% Points in a time dependant coordinate system

% However relativity tells us that even a static point can have varying coordinates if the coordinate system moves. So points stored as "'Points'" must be transformed in "'Time-Points'". Note that in this example Richard has used directly the "psitbxgrid" object "GridTor" instead of the axis radial and vertical positions.

ProbePointsTor = psitbxg2g(ProbePointsOrt,'Toroidal',GridTor)

%%%% Time-Points in a time dependant coordinate system

% If your points are already stored as "'Time-Points'", transformation to a time dependant coordinate system can be achieved only if the number of points in time are the same for both. The rule is that the Psi-Toolbox only handles a unique time base.

MovingProbePointsTor = psitbxg2g(MovingProbePointsOrt,'Toroidal',GridTor)

%%%% Now something simpler: 0-D, 1-D and 2-D coordinate systems

% Now that the difficult situation have been treated, we can come back to simpler things. Very often, the quantity you are dealing with does not depend on all spatial coordinates. For this purpose you can place a "NaN" on the implied spatial direction:
 
rzphi = {0.62:.02:0.7,-.05:.02:.05,NaN};
Grid2DCyl = psitbxgrid('Cylindrical','Grid',rzphi)
plot(Grid2DCyl)
drawnow

%%%% Finaly the ultimate complexity: N-D coordinate systems

% Having gone from 2D to 3D and then to time dependant coordinates, it was nothing for me to generalise. For example if Richard wants to study edge parameter fluctuation at different frequencies, he can define an additional dimension to his already complicated problem:
 
rzphifreq = {0.62:.02:0.7,-.05:.02:.05,pi/180*(-3:2:3),[1e4,2e4,5e4,1e5]};
Grid4DCyl = psitbxgrid('Cylindrical','Grid',rzphifreq)
Grid4DOrt = psitbxg2g(Grid4DCyl,'Orthogonal')
Grid4DOrt.x{4}(1,1,1,:)

%%% Functions of the spatial coordinates and others

%%%% Functions of the spatial position: the psitbxfun class

% It is highly improbable that Richard will spend his life with coordinate systems. So it is time for him to open the Psi-Toolbox and take the "psitbxfun" class. This class will allow Richard to represent quantity that depend on the spatial position. Let us go back to the first grid Richard defined in cylindrical coordinates and suppose for simplicity that Richard's moving probes measure a given simple quantity. A "psitbxfun" object is then constructed by specifying the value of this function at the grid points and the grid where it is evaluated:

rzphi = GridCyl.x
[rzphi{:}] = ndgrid(rzphi{:})
x = exp(-300*(rzphi{1} - GridCyl.x{1}(1)).^2) .* exp(-3000*rzphi{2}.^2);
FunCyl = psitbxfun(x,GridCyl)
plot(FunCyl)
drawnow

%%%% Overloaded math operators for psitbxfun class
 
% Almost all applicable mathematical operators "+ - * / .* ./ \ .\ '", relational operators "== ~= > < >= <=", logical operators "~ & |", mathematical function "sin cos ..." are overloaded for "psitbxfun" objects. Experience:
 
plot(log(10*FunCyl))
drawnow

%%%% Functions of the spatial position and of the time

% Quantities that varies also in time can be represented by 4-D arrays whos first three dimensions correspond to the spatial coordinates and the last one to time samples. The time array is then also passed to the "psitbxfun" constructor and stored in the appropriate field. Of course the last dimension and the length of the time vector must agree.

rzphi = GridCyl.x
rzphit = {rzphi{:},t};
[rzphit{:}] = ndgrid(rzphit{:})
x = exp(-300*(rzphit{1} - GridCyl.x{1}(1)).^2) .* exp(-3000*rzphit{2}.^2) .* exp(3*rzphit{4});
FunTimeCyl = psitbxfun(x,GridCyl,t)

%%%% General forms from 0-D to N-D functions

% As for the grid, generalisation is straitforward. On an N-D grid "GridND" a function can have the size "[size(GridND),np1,np2,...[,nt]]". However I have no real example available yet. In addition if one of the grid coordinate is degenerated, being "NaN" or scalar, the corresponding singleton in the function value can be ommited. One common example are functions on a 2-D grid:
 
rz = Grid2DCyl.x
[rz{:}] = ndgrid(rz{:})
x = ((rz{1} - mean(Grid2DCyl.x{1})).^2 + (rz{2} - mean(Grid2DCyl.x{2})).^2);
Fun2DCyl = psitbxfun(x,Grid2DCyl)

%%%% Function mapping: the psitbxf2f method

% Now this obvious question is: what does Richard measures with his probes.  This can easily be answered using the "psitbxf2f" method: the quantity defined in cylindrical coordinates must be estimated at the position of the probes. Let us start with fixed probes and a quantity constant in time:

FunProbePointsOrt = psitbxf2f(FunCyl,ProbePointsOrt)
FunProbePointsOrt.x

%%%% Mapping of a function on Time-Points

% Now a quantity constant in time measured by moving probes:

FunMovingProbeOrt = psitbxf2f(FunCyl,MovingProbePointsOrt)
plot(FunMovingProbeOrt,'o','markersize',12)
drawnow

%%%% Mapping of a time dependant function on Time-Points

% Finaly a time dependant quantity measured by moving probes:

FunTimeMovingProbeOrt = psitbxf2f(FunTimeCyl,MovingProbePointsOrt)
plot(FunTimeMovingProbeOrt,'o','markersize',12)
drawnow

%%%% A last extra: spatial derivative of functions

% The "psitbxf2f" method can also be used to compute the spatial derivative of functions, by specifying the derivation order along each spatial coordinate. Note that since functions are interpolated by cubic spline, the third derivative is the last non zero derivatives.

DerDRFunTimeMovingProbeOrt = psitbxf2f(FunTimeCyl,MovingProbePointsOrt,[1,0,0])
plot(DerDRFunTimeMovingProbeOrt,'o','markersize',12)
drawnow

%%%% Function integrals: the sum method

% Integration along the definition variables of functions defined on a "'Grid'" is the task of the "sum" method. Summation along a dimension corresponding to one of the grid coordinates or the time results in an integration. Along the other dimensions, it is a simple summation. The dimension argument can be "'t'" with obvious effect. As an example take the "psitbxtcv" function with mode "'JPHI'", which will return a "psitbxfun" object containing the toroidal current distribution.

jphi = psitbxtcv(10000,t,'JPHI')
plot(jphi)
drawnow

% Certainly the most non trivial way to obtain the plasma current is:

Ip = sum(jphi,[1,2])
Ip.grid.x
plot(Ip) 
drawnow

%%%% Function average: the mean method

% Average is the task of the "mean" method, which works more or less like "sum". However it accepts an optional argument, a "psitbxfun" object which will play the role of a weigthing function. In fact in that case "mean" is simply "sum(x.*w)./sum(w)":

javg = mean(jphi,[1,2],jphi>0)
plot(javg)
drawnow

%%% Equilibrium geometry description

% Tokamak equilibrium is very often described with the spatial distribution of the poloidal flux in cylindrical coordinate. Of course since Tokamak is axisymetric, this flux does not depend on the aximuthal coordinate.

%%%% First description of the equilibrium: the psitbxpsi class

% It is pretty obvious that the poloidal flux used to describe the Tokamak equilibrium could be stored in a "psitbxfun" object with cylindrical coordinates. However to add some functionality I defined a new class, the "psitbxpsi" class, child of the "psitbxfun" class. This means that all methods defined for the "psitbxfun" class apply to the "psitbxpsi" class and that a "psitbxpsi" object contains all the fields of a "psitbxfun" object plus others. The first additional field is the "format" field which specifies the way the poloidal flux is stored. It is a two character string each with possible values of "'-'", "'0'" or "'+'". The first character gives the sign of the poloidal flux in the plasma center, the second on the plasma edge. A useless example of a "psitbxpsi" object would be:
 
psi = psitbxpsi(ndgrid(Grid2DCyl.x{1},Grid2DCyl.x{2}),Grid2DCyl,[],'-0')

%%%% Equilibrium description from TCV shot file

% Much more usefull is the description of a real TCV shot. This can be easily built with the "psitbxtcv" function. Here Richard selected times in his favorite standard shot #10000:

psi = psitbxtcv(10000,t)
plot(psi)
drawnow

%%%% Oh how easy it is now

% A long time ago Richard wanted to know the poloidal magnetic field near his tiles, basicaly given by the derivative of the poloidal flux along the radial and vertical directions. At this point you realise how easy it is:

drpsi = psitbxf2f(psi,Grid2DCyl,[1 0])
dzpsi = psitbxf2f(psi,Grid2DCyl,[0 1])
plot(drpsi)
drawnow

%%%% Normalised poloidal flux and magnetic axis: the psitbxp2p method

% Most of the Psi-Toolbox functions requires normalised poloidal flux which will later serve as radial coordinate: A normalised poloidal flux is one on the last close flux surface and zero on magnetic axis. Normalisation is achieved with the method "psitbxp2p" on a "psitbxpsi" object. It locates the magnetic axis position stored in "rmag" and "zmag" fields and modifies the value of psi:

psi = psitbxp2p(psi,'01')
plot(psi)
drawnow

% Note that specifying the format in the "psitbxtcv" call would yield the same result:

psi = psitbxtcv(10000,t,'01')

%%%% The Flux coordinate system

% Let me go back to the coordinate systems and introduce the very important flux coordinate system. This is basicaly a toroidal coordinate system in which the distance from the axis is replaced by the squaroot of the normalised poloidal flux. Although not mandatory, it may be very useful to associate to such a coordinate system the value of the normalised poloidal flux stored in a "psitbxpsi" object:

GridFluxPsi = psitbxgrid('Flux','Grid','Default',psi)
size(GridFluxPsi)

%%%% Quantities depending only on the flux: profiles

% A very common concept in Tokamak physics is the profile, that is a quantity constant on a flux surface and therefor depending only on the flux coordinate. Such a profile is trivialy represented by a function on a 1-D grid, here the canonical parabolic profile:

Grid1DRho = psitbxgrid('Flux','Grid',{linspace(0,1,41),NaN,NaN},psi)
ParabolicRho = psitbxfun(1-linspace(0,1,41)'.^2,Grid1DRho)
plot(ParabolicRho)
drawnow

%%%% Coordinate transformation to flux coordinate system

% Of course discovering this new coordinate system, Richard wants to know where his moving probes lies in term of flux coordinate. The "psitbxg2g" method applies. Transformation from any coordinate system to flux coordinate requires a normalised poloidal flux description:

MovingProbePointsFlux = psitbxg2g(MovingProbePointsOrt,'Flux',psi)
MovingProbePointsFlux.x{1}

%%%% Flux surface description

% An alternative but very important representation of the equilibrium geometry in the Psi-Toolbox is the description of the flux surfaces themselves. These are contours of constant poloidal flux; this is formalised as the distance from the magnetic axis as a function of the squaroot of the normalised flux and the poloidal angle. It can be stored in a "psitbxfun" object with a flux coordinate grid. But here once more to get addional functionality, it is stored in a "psitbxpsi" object with format "'FS'". Conversion from a "psitbxpsi" object containing the distribution of the poloidal flux is achieved with the "psitbxp2p" method:

fsd = psitbxp2p(psi,'FS')
plot(fsd)
drawnow

% Now try to catch this:

plot(fsd.psitbxfun)
drawnow

%%%% Flux surface description from the TCV shot file

% Richard will very soon realise that the "psitbxp2p" method takes a lot of time. To spare Richard's time, the necessary information is stored in the TCV shot file once the equilibrium has been reconstructed. Thus a fast way to recover the flux surface contour for an old shot is to use the "psitbxtcv" function with format "'FS'":
 
fsd = psitbxtcv(10000,t,'FS')

%%%% Coordinate transformation from flux coordinate system

% The first usage of a "psitbxpsi" object of format "FS" is for coordinate transformation from flux coordinates:
 
rhothetaphi = {[0,.5,1],linspace(-pi,pi,9),linspace(0,pi/2,4)};
plot(psitbxg2g(psitbxgrid('Flux','Grid',rhothetaphi),'Orthogonal',fsd),'-or','markersize',12)
drawnow

%%%% Flux coordinate metric

% The main difficulty with flux coordinates is that it is not an orthogonal coordinate system, which renders its metric non trivial. The Psi-Toolbox provides a method for taking the appropriate metric into account, simply called "metric". It applies to a "psitbxgrid" object with "'Grid'" storage or a "psitbxpsi" object with format "'FS'" and takes as second argument the combination of dimension for which you need the metric. For example if you want the metric on the flux surface, use "[2,3]" to get the unit area on (theta,phi) planes; for unit volume, use "[1,2,3]". You can also get the norm of the coordinate gradient with a corresponding negative dimension. In the example the very important gradient of rho entering the expression of the SEF (Shape Enhancement Factor) is plotted. This method returns a "psitbxfun" object and also works for coordinate systems with diagonal metric.
 
GradRho = metric(fsd,-1)
plot(GradRho)
drawnow

%%%% Function integrals with metric

% Using the "sum" method with the option "'metric'" allows to perform spatial integrals with the correct metric. In the example the toroidal current distribution in cylindrical coordinates is first mapped in flux coordinates and then surface averaged:

GridFluxFSD = psitbxgrid('Flux','Grid','Default',fsd)
jphiFlux = psitbxf2f(jphi,GridFluxFSD)
IpFlux = sum(jphiFlux,[1,2],'metric')
plot(Ip), hold on, plot(IpFlux,'-og'), hold off
drawnow

%%%% Function average with metric

% The "mean" method can also be used with a "'metric'" option. The simplified definitions are: <DL><DT> "mean(f,dim,w,'metric')" <DD> "sum(f.*w.*metric(f.grid,dim),dim)./sum(w.*metric(f.grid,dim),dim)"</DL>

mean(jphiFlux,[1,2],'metric')
mean(jphi,[1,2],jphi > 0)

%%%% Flux surface geometrical parameters

% Quantities such as the flux surface area, cross section, inside volume or the volume or cross-section between two adjacent flux surfaces are widely used. Also very important since it enters in the expression of the SEF is the flux surface average of the gradient of rho. They are now obvious to compute using the following definitions: <DL><DT> "grho" flux surface average of the gradient of rho <DD> mean(metric(g,-1),[2,3]) <DT> "surf" flux surface area <DD> "sum(metric(g,[2,3]),[2,3])" <DT> "area" flux surface cross-section <DD> "cumsum(sum(metric(g,[1,2]),2),1)" <DT> "vol" flux surface inside volume <DD> "cumsum(sum(metric(g,[1,2,3]),[2,3]),1)" <DT> "darea" cross-section between two flux surfaces divided by the flux coordinate increment <DD> "sum(metric(g,[1,2]),2)" <DT> "dvol" volume between two flux surfaces divided by the flux coordinate increment <DD> "sum(metric(g,[1,2,3]),[2,3])" </DL> All these quantities are returned in a structure by the method "psitbxfsg":

fsg = psitbxfsg(fsd)
plot(fsg.grho)
drawnow

%%% Diagnostic chords

%%%% Diagnostic Chord Description: the psitbxdcd class

% Describing the geometry of a multi-chord diagnostic require to gather some information. This will be arranged in a "psitbxdcd" object with the following properties:

% <DL COMPATCT>
% <DT>"rd"<DD>The radial position in cylindrical coordinate of a reference point along the chord, very often the pin-hole
% <DT>"zd"<DD>The vertical position of the same point <DT>"phid"<DD>The toroidal angle of the same point
% <DT>"pvd"<DD>The poloidal viewing angle, measured between the horizontal plane and the chord, positive if the chord is pointing above the plane 
% <DT>"tvd"<DD>The toroidal viewing angle, measured between the vertical plane containing the machine axis and the chord reference point and the vertical plane containing the chord itself, positive if the chord is looking in the direction of increasing toroidal angle
% <DT>"nd"<DD>The number of points used to discretise the chord
% <DT>"t"<DD>The time base for moving chords
% </DL>
 
% Here is an example with two "cameras" with 20 detectors each, the first one with a poloidal coverage, the second one with a toroidal coverage. Note that the fields of a "psitbxdcd" object can retain an appropriate organisation, here 2 by 20, or 1028 by 1028 for CCD fans.

rd = 1.2;
zd = 0;
phid = pi/8;
pvd = [linspace(-pi/3,pi/3,20);zeros(1,20)];
tvd = [zeros(1,20);linspace(-pi/3,pi/3,20)];
dcd = psitbxdcd(rd,zd,phid,pvd,tvd);
plot(dcd)
drawnow

%%%% Moving chords

% The possibility to have real-time scaning diagnostic or beam is included. In that case the last dimension of the fields of the "psitbxdcd" object represents the time base. Given as an example <CITE>Le cauchemar du gyrotroniste</CITE>:

td = linspace(0,1,73);
rd = 1.2;
zd = 0;
phid = pi/2;
pvd = pi/6*cos(pi*td);
tvd = pi/6*sin(2*pi*td);
LeCauchemarDuGyrotroniste = psitbxdcd(rd,zd,phid,pvd,tvd,[],td);
plot(LeCauchemarDuGyrotroniste)
drawnow

%%%% Diagnostic-Chords coordinate

% In essence a diagnostic chord is an object with only one dimension, defined as the distance from the reference point. Although certainly possible, I have not left the possibility to have a constant grid along all the chords. Only "'Points'" or "'Time-Points'" storage is possible. The size of the array representing the chord coordinate must be "[dcd.nd,size(dcd.rd)]". Since it is somehow fastidious to chose a good grid for each chord, I have made the "intersect" method which returns the chord coordinate so that it covers the intersection of the chords and a cylindrical coordinate grid.

GridDiagChords = psitbxgrid('Diagnostic-Chords','Points',{intersect(dcd,psi.grid)},dcd)
plot(psitbxgrid('Cylindrical','Grid',{[.614,1.147],[-.76,.76],linspace(-pi/2,pi,10)}))
hold on, plot(GridDiagChords,'-b'), hold off
drawnow

%%%% Coordinate transformation from Diagnostic-Chords coordinates

% Any coordinate transformation from Diagnostic-Chords coordinate is possible. The inverse is not possible since chords are discrete and do not cover the whole space.

GridDiagChordsFlux = psitbxg2g(GridDiagChords,'Flux',psi)
subplot(211)
plot(squeeze(GridDiagChords.x{1}(:,1,:)),squeeze(GridDiagChordsFlux.x{1}(:,1,:,1)))
subplot(212)
plot(squeeze(GridDiagChords.x{1}(:,2,:)),squeeze(GridDiagChordsFlux.x{1}(:,2,:,1)))
xlabel('S'), ylabel('Rho')
drawnow
subplot

%%%% Function mapping on a Diagnostic-Chords coordinates

% Mapping a function of the position on the chords of a diagnostic is a very essential operation in the Psi-Toolbox. The example shows the value of the canonical parabolic profile along the various chords of the diagnostic. Note that the value at the points along the chord and outside the definition domain of the function are set to "NaN".

ParabolicRhoDiagChords = psitbxf2f(ParabolicRho,GridDiagChords)
plot(ParabolicRhoDiagChords,'-')
drawnow

% or in another representation

plot(ParabolicRhoDiagChords)
drawnow

%%%% Line integrals and averages

% A function defined along Diagnostic-Chords can be integrated or averaged along the line with the "sum" and "mean" methods. In the example the line integral and average of the parabolic profile seen by each camera at different times is ploted as a function of the detector number.

ParabolicRhodInt = sum(ParabolicRhoDiagChords,1)
ParabolicRhodAvg = mean(ParabolicRhoDiagChords,1)
subplot(121), plot(reshape(permute(ParabolicRhodAvg.x,[2,1,3]),[20,6]))
subplot(122), plot(reshape(permute(ParabolicRhodInt.x,[2,1,3]),[20,6]))
drawnow
subplot

%%%% Transfer matrix and image reconstruction

% The transfer matrix from a local emissivity defined on a cylindrical grid and chords can also be computed by the "tmat" method. In the following example an attractive plasma emitting the &psi; letter is seen by a 64 by 64 pixel CCD camera:

t = linspace(0,60,64)/180*pi;
p = linspace(-30,30,64)'/180*pi;
tvd = repmat(t,length(p),1);
pvd = repmat(p,1,length(t));
dcd = psitbxdcd(1,0,0,pvd,tvd);
grid = psitbxgrid('c','g',{linspace(.3,.9,21),linspace(-.6,.6,41),NaN});
T = tmat(dcd,grid);
x = repmat(' ',20,40);
x(:,15:32) = [
' ++     ++     ++   '
'  +++   ++   +++    '
'   ++   ++   ++     '
'   +++  ++  +++     '
'   +++  ++  +++     '
'   +++  ++  +++     '
'   +++  ++  +++     '
'   +++  ++  +++     '
'   +++  ++  +++     '
'    ++  ++  ++      '
'    ++  ++  ++      '
'     ++ ++ ++       '
'      ++++++        '
'        ++          '
'        ++          '
'        ++          '
'        ++          '
'        ++          '
]';
x = x == '+';
x = reshape(T*x(:),size(dcd));
set(image(x),'cdatamapping','scaled'), axis equal off
drawnow

% "intersect" may also precalculate matrices that allow fast reconstruction of the emmissivity profile.

% END
