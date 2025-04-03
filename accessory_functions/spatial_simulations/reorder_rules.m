function rule_order = reorder_rules(model)
% model: SimBiology model object

% find a permutation for the order of the rules when we can caluclate it
% e.g. if we need the 2nd rule to calculate the 1st one, then we should
% calculate the 2nd first
% number of the rules
nrule = numel(model.Rules);
% variables set by the rules
rule_names = cell(nrule,1);
for i = 1:nrule
    % find "=" sign
    pos = find(model.Rules(i).Rule=='=',1,'first')-2;
    rule_names{i} = model.Rules(i).Rule(1:pos);
end
% build a matrix reflecting the "interactions" between the rules
Mrule = zeros(nrule,nrule);
for i = 1:nrule
    for j = 1:nrule
        if i~=j
            Mrule(i,j) = contains(model.Rules(i).Rule,rule_names{j});
        end
    end
end

% reorder ithe matrix to get an under triangular matrix
% generate increasing weights when the new weight is bigger then the sum of
% every weight before it
weights = zeros(1,nrule);
weights(1) = 1;
weight_sum = 1;
for i = 2:nrule
    weight_sum = weight_sum+weights(i-1);
    weights(i) = weight_sum+1;
end
% give weights to the matrix and sum it to the rows
B = sum(Mrule.*weights,2);
% sort it according to the weights
[~,rule_order] = sort(B);

end

