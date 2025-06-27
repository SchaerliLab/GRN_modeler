function [Ecoli] = build_Tomazou_repressilator(model_name,N)
%BUILD_REPRESSILATOR build the repressilator model
% model_name: the name of the selected model
% N: the number of nodes

% build repressilator
[Ecoli] = build_repressilator(model_name,N);

% add a common protease to every node
for i = 1:N
    Ecoli = Ecoli.add_protease(['N',int2str(i)],'PROT1','type1');
end

end

