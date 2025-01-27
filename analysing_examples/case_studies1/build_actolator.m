function [Ecoli] = build_actolator(model_name,N)
%BUILD_REPRESSILATOR build the repressilator model
% model_name: the name of the selected model
% N: the number of nodes

if rem(N,2)~=0
    error('The number of nodes must be even in the actolator!')
end

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
    Ecoli = Ecoli.add_regulator('Activation_in',['N',int2str(i+1)],'HILL',['P_N',int2str(i)]);
end
% close the loop
Ecoli = Ecoli.add_regulator('Activation_in','N1','HILL',['P_N',int2str(N)]);

% modify it with toggle swithes
% add regulators
for i = 1:N/2
    Ecoli = Ecoli.add_regulator('Repression_in',['N',int2str(i+N/2)],'HILL',['P_N',int2str(i)]);
    Ecoli = Ecoli.add_regulator('Repression_in',['N',int2str(i)],'HILL',['P_N',int2str(i+N/2)]);
end

end

