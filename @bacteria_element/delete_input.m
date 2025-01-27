function [] = delete_input(obj,Mobj)
% delete input ('HILL' parameters) from the reactions and rules in the Mobj model

% we will go through on every possible 'HILL' input of the object
for n_hill = 1:numel(obj.Mobj.UserData.input)

    % % the name of the actual HILL parameter
    % Hill_target = obj.Mobj.UserData.input{n_hill}; % e.g. 'HILL1'

    % the name of the actual HILL parameter
    % rename the hill target
    if isa(obj,'node') || isa(obj,'protease')
        Hill_target = [obj.Mobj.UserData.input{n_hill},'_',obj.name,'-input'];
    else
        Hill_target = [obj.Mobj.UserData.input{n_hill},'_',obj.suffix_simple{:},'-input'];
    end


    % delete *hillname if the next character is not '_'.
    % e.g.: *HILL_A_R1_R2*HILL_A => *HILL_A_R1_R2
    pattern = ['*[' Hill_target ']'];
    % reactions
    for i = 1:numel(Mobj.Reactions)
        set(Mobj.Reactions(i),'Reactionrate',strrep(get(Mobj.Reactions(i),'Reactionrate'),pattern,''))
    end
    % rules
    for i = 1:numel(Mobj.Rules)
        set(Mobj.Rules(i),'Rule',strrep(get(Mobj.Rules(i),'Rule'),pattern,''))
    end

    delete(sbioselect(Mobj.Parameters,'Name',Hill_target));

    % % delete *hillname if the next character is not '_'.
    % % e.g.: *HILL_A_R1_R2*HILL_A => *HILL_A_R1_R2
    % pattern = [hill_name '(?!_)'];
    % % reactions
    % for i = 1:numel(Mobj.Reactions)
    %     set(Mobj.Reactions(i),'Reactionrate',regexprep(get(Mobj.Reactions(i),'Reactionrate'),pattern,''))
    % end
    % % rules
    % for i = 1:numel(Mobj.Rules)
    %     set(Mobj.Rules(i),'Rule',regexprep(get(Mobj.Rules(i),'Rule'),pattern,''))
    % end
    % 
    % hill_name = hill_name(2:end); % hill_name without '*'
    % delete(sbioselect(Mobj.Parameters,'Name',hill_name));
end

end
