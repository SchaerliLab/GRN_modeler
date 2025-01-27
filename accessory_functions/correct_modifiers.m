function Mobj = correct_modifiers(Mobj_in)
%CORRECT_MODIFIERS is SBML files species in the rate law should set as
%modifiers. Here we will add them to both sides
% e.g. null -> A (r = k*b) => B -> A + B

% copy the simbiology object
Mobj = copyobj(Mobj_in);

% add an extra "-" to the end of the species, so in the ratelaw they will
% have a bracket around them [specname-], and we can find them
% even though if the names contains each other (e.g. A, AD, AAD etc)
for i = 1:numel(Mobj.Species)
    rename(Mobj.Species(i),[Mobj.Species(i).Name,'-']);
end

% every species in the model
allSpecies = {Mobj.Species.Name};

% go through on every reactions
for i = 1:numel(Mobj_in.Reactions)

    % get species in the reaction
    reactionSpecies = getSpeciesInReaction(Mobj.Reactions(i));

    % species included in the reaction rate but not in the reaction
    speciesNames = extractSpeciesNames(Mobj.Reactions(i), allSpecies, reactionSpecies);

    % if we do not have this type of species
    if isempty(speciesNames)
        continue;
    end

    % the new species in the reaction will not change its amount
    % it will not be included on the differential equations
    new_spec = strcat({' + '},speciesNames(:));
    new_spec = [new_spec{:}];

    if Mobj.Reactions(1).Reversible == true
        sign = '<->';
    else
        sign = '->';
    end

    % get the reaction complexes
    complexes = strsplit(Mobj.Reactions(i).Reaction,sign);
    for j = 1:2
        % get rid of null if we have it
        complexes{j} = strrep(complexes{j},'null','');
        % add the new species
        complexes{j} = [complexes{j},new_spec];
        % get rid of white spaces at the beginning and at the end
        complexes{j} = strtrim(complexes{j});
        % it should not start with a + sign
        if complexes{j}(1) == '+'
            complexes{j}(1) = [];
        end
    end

    %rebuild the ratelaw
    Mobj.Reactions(i).Reaction = [complexes{1},' ',sign,' ',complexes{2}];

end

% delete the extra "-" character from the end of the species names
% what we have added at the beginning
for i = 1:numel(Mobj.Species)
    rename(Mobj.Species(i),Mobj.Species(i).Name(1:end-1));
end

end

% collect the species which are present in the reactionrate but not in the
% reaction
function speciesNames = extractSpeciesNames(reactionRate, allSpecies, reactionSpecies)
    % Initialize an empty cell array to store species names
    speciesNames = {};
    
    % Loop through all species and check if they are in the reaction rate but not in the reaction
    for i = 1:length(allSpecies)
        speciesName = allSpecies{i};
        % Check if the species name is present in the reaction rate and not in the reaction species
        if contains(reactionRate.ReactionRate,['[',speciesName,']']) && ~any(strcmp(reactionSpecies,speciesName))
            speciesNames{end+1} = speciesName; %#ok<AGROW>
        end
    end
end

% get every species in the reaction
function reactionSpecies = getSpeciesInReaction(reaction)
    % Get reactants and products
    reactants = get(reaction.Reactants, 'Name');
    products = get(reaction.Products, 'Name');

    if ~iscell(reactants) && ~isempty(reactants)
        reactants = {reactants};
    end
    if ~iscell(products) && ~isempty(products)
        products = {products};
    end
    
    if isempty(reactants)
        reactionSpecies = products;
    elseif isempty(products)
        reactionSpecies = reactants;
    else
        reactionSpecies = [reshape(reactants,[1,numel(reactants)]),...
            reshape(products,[1,numel(products)])];
    end
    
end