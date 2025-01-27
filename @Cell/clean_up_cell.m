function cell = clean_up_cell(cell)
% clean up every deleted model element in every object

% nodes
node_names = fieldnames(cell.nodes);
for i = 1:numel(node_names)
    cell.nodes.(node_names{i}) = cell.nodes.(node_names{i}).clean_up_obj();
end
% proteases
protease_names = fieldnames(cell.proteases);
for i = 1:numel(protease_names)
    cell.proteases.(protease_names{i}) = cell.proteases.(protease_names{i}).clean_up_obj();
end
% regulators
regulator_names = fieldnames(cell.data.regulator_routes);
for i = 1:numel(regulator_names)
    regulation = cell.data.regulator_routes.(regulator_names{i});
    for j = 1:numel(regulation)
        % find the object
        regulated_obj = eval(cell.route2fieldnames(regulation{j}));
        regulated_obj = regulated_obj.clean_up_obj();
    end
end

end

