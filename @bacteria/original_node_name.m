function original_name = original_node_name(obj,input_name)
% calculate the original name of a node without a prefix

name = strsplit(input_name,'_');
for i = numel(name):-1:1
    original_name = strcat(name{i:end},'_');
    % original_name = [original_name{:}];
    original_name(end) = [];

    if obj.isNode(original_name)
        % if the name is ['node_prefix','original_name'] then we are ready
        proposed_name = [obj.nodes.(original_name).Mobj.UserData.output_prefix,original_name];
        proposed_name = [proposed_name{:}];
        if strcmp(proposed_name,input_name)
            return;
        end
    end
end

% if we did not found another name
original_name = input_name;