function [Ecoli] = build_repressilator(model_name,N)
%BUILD_REPRESSILATOR build the repressilator model
% model_name: the name of the selected model
% N: the number of nodes

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
    Ecoli = Ecoli.add_regulator('Repression_in',['N',int2str(i+1)],'HILL',['P_N',int2str(i)]);
end
% close the loop
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL',['P_N',int2str(N)]);

% % follow every protein
% for i = 1:N
%     Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {['P_N',int2str(i)]}];
% end

end

