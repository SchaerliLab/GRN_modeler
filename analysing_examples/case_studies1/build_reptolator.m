function [Ecoli] = build_reptolator(model_name,N)
%BUILD_REPTOLATOR build the reptolator model
% model_name: the name of the selected model
% N: the number of nodes

if rem(N,2)~=0
    error('The number of nodes must be even in the reptolator!')
end

% fist build the repressilator
[Ecoli] = build_repressilator(model_name,N);

% modify it with toggle swithes
% add regulators
for i = 1:N/2
    Ecoli = Ecoli.add_regulator('Repression_in',['N',int2str(i+N/2)],'HILL',['P_N',int2str(i)]);
    Ecoli = Ecoli.add_regulator('Repression_in',['N',int2str(i)],'HILL',['P_N',int2str(i+N/2)]);
end

end

