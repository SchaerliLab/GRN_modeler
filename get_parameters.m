function p = get_parameters(Mobj,name)
%GET_PARAMETERS get the constant parameters from the model and put them
%into a structure
% if name is provided, rename them like name1, name2, etc
% otherwise try to use the SimBiology names (it is not always possible)

n_ext = 0;
for i = 1:numel(Mobj.Parameters)
    if Mobj.Parameters(i).Constant == true
        if nargin == 1
            p.(Mobj.Parameters(i).Name) = Mobj.Parameters(i).Value;
        else
            n_ext = n_ext+1;
            p.([name, int2str(n_ext)]) = Mobj.Parameters(i).Value;
        end
    end
end

