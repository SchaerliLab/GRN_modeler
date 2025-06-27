function Ecoli = build_double_repressilator(model_name,N1,N2)
%BUILD_DOUBLE_REPRESSILATOR build the repressilator model
% model_name: the name of the selected model
% N1: the number of nodes in the first ring
% N2: the number of nodes in the second ring

% build a repressilator
Ecoli = build_repressilator(model_name,N1);

% if the second loop does not contain new nodes,
% there is nothing to do
if N2<3
    return
end

% create the second ring
% N1 and N2 is part of this ring, we need less extra nodes
% add the nodes
for i = N1+1:N1+N2-2
    Ecoli = Ecoli.add_node(['N' int2str(i)],'type1');
end

% add regulators
for i = N1+1:N1+N2-3
    Ecoli = Ecoli.add_regulator('Repression_in',['N',int2str(i+1)],'HILL',['P_N',int2str(i)]);
end
% close the loop
Ecoli = Ecoli.add_regulator('Repression_in',['N',int2str(N1+1)],'HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL',['P_N',int2str(N1+N2-2)]);