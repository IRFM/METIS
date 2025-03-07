function out = copy_composition_to_compositionstype(in)

Nspec=length(in.amn);

for j=1:Nspec
  comps.nuclei(j).amn=-1;
  comps.nuclei(j).zn=-1;
end

% Fill in comps().nuclei
Nnuclei = 0;
for j=1:Nspec
       tmp_nuc.amn = in.amn(j);
       tmp_nuc.zn  = in.zn(j);
       if j==1 
          Nnuclei=Nnuclei+1;
          comps.nuclei(Nnuclei).amn = in.amn(j);
          comps.nuclei(Nnuclei).zn  = in.zn(j);
       else
          if find_same_nuclei(tmp_nuc, comps.nuclei ) == -1
             Nnuclei=Nnuclei+1;
             comps.nuclei(Nnuclei).amn = in.amn(j);
             comps.nuclei(Nnuclei).zn  = in.zn(j);
          end
       end
end

for j=1:Nnuclei
       out.nuclei(j).amn = comps.nuclei(j).amn;
       out.nuclei(j).zn  = comps.nuclei(j).zn;
end

for j=1:Nspec
       % Dummy nuclei; used only to find the index in out.nuclei
       tmp_nuc.amn = in.amn(j);
       tmp_nuc.zn  = in.zn(j);

       out.ions(j).zion   = in.zion(j);
       out.ions(j).nucindex  = find_same_nuclei( tmp_nuc , out.nuclei );
end


function rep = find_same_nuclei(arg1,arg2)

if ~isempty(arg2)
       for j=1:length(arg2)
          if same_nuclei(arg1,arg2(j))
             rep = j;
             return
          end
       end
end
rep = -1;  % same nuclei not found


%> Are the arg1 and arg2 the same nuclei?
function rep = same_nuclei(arg1,arg2)

charge_tol  = 0.01;
mass_tol    = 0.01;
rep = (abs(arg1.amn - arg2.amn) < mass_tol) & ...
              (abs(arg1.zn  - arg2.zn ) < charge_tol);

 
 
