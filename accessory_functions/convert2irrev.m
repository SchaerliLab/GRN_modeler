function Mobj = convert2irrev(Mobj_in)
%CONVERT2IRREV convert the reactions to irreversible reactions
% we are just separating naively the positive and negative terms

% copy the simbiology object
Mobj = copyobj(Mobj_in);

% go through on every reactions
for i = 1:numel(Mobj_in.Reactions)
    % if it is a reversible reaction
    if Mobj_in.Reactions(i).Reversible == true

        % the rate equation
        rate = Mobj_in.Reactions(i).ReactionRate;

        % they should not be in a species name surrounded with []
        % e.g. [HILL_G|-R3] or in a bracket e.g. the plus sign in
        % -(k_mRNA_degr+dilution)*mRNA_G
        positive = true(size(rate));
        pos_sign = true;

        % how many curvy () and rectengular [] brackets are open
        % if there are no open ones, we can change the sign
        open_curvy = 0;
        position = 1;
        while position <= numel(rate)

            % check the brackets
            if rate(position) == '('
                open_curvy = open_curvy+1;
            elseif rate(position) == ')'
                open_curvy = open_curvy-1;
            elseif rate(position) == '['
                % it is around spetial species names, we will go to the end
                while rate(position) ~= ']'
                    % we store the actual sign
                    positive(position) = pos_sign;
                    position = position+1;
                end
            else
                % if every bracket is closed
                if open_curvy==0
                    % if we have a sign
                    if rate(position) == '+'
                        pos_sign = true;
                    elseif rate(position) == '-'
                        pos_sign = false;
                    end
                end
            end
            % we store the actual sign
            positive(position) = pos_sign;
            % next position in the reaction rate
            position = position+1;
        end

        % forward reaction rate
        rate_forward = rate(positive);
        if ~isempty(rate_forward) && rate_forward(1) == '+'
            rate_forward(1) = [];
        end
        % backward reaction rate
        rate_backward = rate(~positive);

        % create two irreversible reactions
        complexes = strsplit(Mobj_in.Reactions(i).Reaction,'<->');
        % forward reaction
        if ~isempty(rate_forward)
            addreaction(Mobj,[complexes{1} '->' complexes{2}],'ReactionRate',rate_forward,'Name',[Mobj_in.Reactions(i).Name,'_forward']);
        end
        % backward reaction
        if ~isempty(rate_backward)
            addreaction(Mobj,[complexes{2} ' -> ' complexes{1}],'ReactionRate',['-(' rate_backward ')'],'Name',[Mobj_in.Reactions(i).Name,'_backward']);
        end

        % delete the original reaction
        delete(sbioselect(Mobj.Reactions,'Name',Mobj_in.Reactions(i).Name));

    end
end

