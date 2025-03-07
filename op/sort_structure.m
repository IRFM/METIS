function structure_out = sort_structure(structure)
%structure_out = sort_structure(structure)
%This function sorts the entries in a
%structure in an ascending order.
%Structures within the structure are sorted
%seperately. They are put before the other
%fields. They themselves are also sorted.
%Input:     -Structure, with or without substructures
%Output:    -Structure, sorted
%Written on 12th july 2002 by Mark Jobse


%Declarations
structure_out = [];
struct_cell = [];
non_struct_cell = [];

%Get all fields out of the structure
cell = fieldnames(structure);

%Now split those fields into structures and other fields
counter = 1;
counter2 = 1;
for i = 1:length(cell)
    boolean = isstruct(getfield(structure, cell{i}));
    if (boolean)
        struct_cell{counter} = cell{i};
        counter = counter + 1;
    else
        non_struct_cell{counter2} = cell{i};
        counter2 = counter2 + 1;
    end
end

%Sort both cells
struct_cell = sort(struct_cell);
non_struct_cell = sort(non_struct_cell);

%Now sort all entries in the found structures (with this function)
for i = 1:length(struct_cell)
    value = sort_structure(getfield(structure, struct_cell{i}));
    structure_out = setfield(structure_out, struct_cell{i}, value);
end

%Now do the same for the non-structure fields.
for i = 1:length(non_struct_cell)
    value = getfield(structure, non_struct_cell{i});
    structure_out = setfield(structure_out,non_struct_cell{i}, value);
end