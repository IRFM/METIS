__________________________________________________________________________
"StructBrowser_gui_g" v1.1 for MatLab 6.5, 
GUI for browsing any structure and plotting its fields


Syntax: StructBrowser_gui_g(cell_struct)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H. Lahdili,(hassan.lahdili@crc.ca)
% Communications Research Centre (CRC) | Advanced Audio Systems (AAS)
% www.crc.ca | www.crc.ca/aas
% Ottawa. Canada
%
% CRC Advanced Audio Systems - Ottawa © 2002-2003
% 16/05/2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__________________________________________________________________________


0. Introduction.



The Graphics User Interface (GUI) is designed to fill a gap in
the MatLab 6.5 Array Editor which pops up when you double click
on a struct array in the MatLab workspace (or using openvar command). 
The Array editor exposes the surface structure (top level structure) 
of the struct array but it does not allow you to delve any deeper.
In order to expose the lower layers of the array struct, it
is necessary to abandon this array editor and type in commands
in the MatLab command window. This can be tedious when you
are doing this many times while running the MatLab debugger.



The CRC StructBrowser, was designed to expose the contents
of all variables of class struct to any depth and plot any 
of the components.



1. Launching the GUI


Ensure that the CRC StructBrowser sources are in your path.
The easiest way of starting it up is to type the command
"StructBrowser" (without the double quotes) while the workspace
contains some information (the same command can be used in the 
debug mode).
Alternatively, you can also start it less directly by typing
a sequence of commands described below. It is unfortunately,
necessary to retype StructBrowser whenever the workspace
has been altered since the MatLab script does not automatically
update. However, you do not need to terminate the StructBrowser
each time you restart it.



The GUI can also be launched by typing "StructBrowser_gui_g(cell_struct)" 
in the MatLab command. cell_struct is a cell of size (2 X N), where N
is the number of structures to browse. The first row of cell_struct
contains all the structures' names, and the second row contains all
the corresponding values. Assuming your workspace contains the 3
structures  struct_1, struct_2 and struct_3, (you can verify this
by typing the MatLab command, "workspace" or clicking the
view menu -> workspace item) then you set cell_stuct as follows:



cell_struct = {'struct_1', 'struct_2', 'struct_3'; ...
                   struct_1, struct_2, struct_3};


and then start the StructBrowser with the command


StructBrowser_gui_g(cell_struct).




In the case of the base workspace, the function "copy_all_struct"
is provided to copy all variables of class struct from workspace 
and store them in a cell array. The syntax of this function is:


cell_struct = copy_all_struct;


and then launch StructBrowser using the command


StructBrowser_gui_g(cell_struct).


Note this mode does not work in the debugger.





2. Browsing the structures


The GUI loads all the structures in the cell_struct array (or all
the structures in the workspace) and lists their names into the 
listbox which appears on the left side. By selecting any structure 
from the listbox, all its fields are displayed in the contents list.
If one of the fields is of class struct, a '+' is appended to the 
fields in the left listbox.
Double clicking on the name of this field in the left listbox,
will expand this structure ('+' turn to '-') and all its subfields are 
listed in a tree format.  All the fields and subfields down to the
lowest level can be displayed. At the lowest level, the fields'
contents are displayed in the contents listbox.  This includes
the size in bytes and the class of the item selected.


In the case where the variable selected is an array of structures (1xN),
double-clicking this selection will expand it in a tree format, showing
the name of the structure's array with subscript (1) to (N).



3. Plot control


3.0  Controls



The radio buttons and buttons appearing in the box below the listbox
on the left provide controls for plotting the contents of the right
listbox.



3.1. Plot selected pushbutton



By pressing this control, the item selected (assuming it has numeric
values) is plotted.  Two options are available: Plot in the 
same figure (Clear figure is checked) or Plot in a new figure
(New figure is checked), so one can keep old plots intact.  If more
than one item in the left listbox are selected, the radio buttons
(SubPlot and Hold on) are activated. The user, can then choose 
between subplotting all the items selected or plotting them in 
a superposed form.



3.2. Close all control



If the"New figure" option was chosen, many figures might be 
generated leading to a very cluttered screen.  The user can close
all the figures by pressing close all pushbutton.




4. Context menu



The plot control is also available by right-clicking any selection in 
the left listbox. A context menu appears: "plot selected" 



5. Demo file



A demo file (sbrowse_demo.m) is provided to make the user familiar with 
the GUI. By typing sbrowse_demo in the MatLab command, two variables of
class struct (struct1 & struct2) are created in the workspace and the
StructBrowser GUI is launched.



6. Copying agreement



We are providing this package as freeware under the
GNU General Public License with the added conditions
that the CRC logo is not removed from the GUI and
credit remains with our laboratory.


Thus you may copy, modify and freely distribute
the program for noncommercial applications.


We assume no liabilities for any damages arising from
the use of this MatLab script.


