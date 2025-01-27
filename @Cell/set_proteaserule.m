function [] = set_proteaserule(cell,prot_name)
%SET_PROTEASERULE set the rule for the substrates for a given protease

rule = ['Substrates_',cell.proteases.(prot_name).name,' = '];
if isempty(cell.proteases.(prot_name).species_list)
    rule = cat(2,rule,'0');
else
    rule_substrates = [];
    for i = 1:numel(cell.proteases.(prot_name).species_list)
        rule_substrates = cat(2,rule_substrates,['+' cell.proteases.(prot_name).species_list{i}]);
    end
    % it should not start with a "+"
    rule_substrates(1) = [];
    rule = cat(2,rule,rule_substrates);
end
set(sbioselect(cell.proteases.(prot_name).Mobj.Rules,'Name',...
    ['substrate_rule_' cell.proteases.(prot_name).name]),'Rule',rule);

end

