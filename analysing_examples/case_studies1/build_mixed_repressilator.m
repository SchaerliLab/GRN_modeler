function [Ecoli] = build_mixed_repressilator(model_name,N,A)
%BUILD_MIXED_REPRESSILATOR build the mixed repressilator model,
% when we have activations in given positions
% model_name: the name of the selected model
% N: the number of nodes
% A: vector for the number of the nodes which will activate the next one

% clean start
clean_up_GRN

% create the model
Ecoli = Cell(model_name);
% add the nodes
for i = 1:N
    Ecoli = Ecoli.add_node(['N' int2str(i)],'type1');
end
% add regulators
for i = 1:N-1
    if any(A==i)
        Ecoli = Ecoli.add_regulator('Activation_in',['N',int2str(i+1)],'HILL',['P_N',int2str(i)]);
    else
        Ecoli = Ecoli.add_regulator('Repression_in',['N',int2str(i+1)],'HILL',['P_N',int2str(i)]);
    end
end
% close the loop
if any(A==N)
    Ecoli = Ecoli.add_regulator('Activation_in','N1','HILL',['P_N',int2str(N)]);
else
    Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL',['P_N',int2str(N)]);
end

% follow every protein
for i = 1:N
    Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {['P_N',int2str(i)]}];
end

end

