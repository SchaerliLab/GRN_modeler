function [Mobj,types] = Tomazou(Mobj,obj_type,n_type)
% Mobj: the simbiology model object which will contain our model
% obj_type: type of the object ('node','protease','regulator')
% n_type: determines the type of the given object

% types: cell array, define the types of the model types, n_type is the 
% number selects the actual type from this cell array

% 'Notes': 'Individual' or 'Common'
% Every species, parameter, reaction, rule must have a 'Note'
% 'Individual': the parameter will be renamed, and different objects (e.g.
% nodes) will have their own parameters
% 'Common': common parameters between every object (e.g. dilution)

% 'Tag': 'INPUT'
% the parameter will be changed by the regulators (e.g. 'HILL' parameters)
% (it should have an 'Individual' 'Note' and it should not be constant)

% nodes:
% 'Tag': 'REGULATOR'
% The species which can act as a regulator

% protease:
% 'Tag': 'PROTEASE'
% The sepecies which will be used as a protease
% 'Tag': 'SUBSTRATE'
% The parameter which will count the sum of the substrates of the protease

% regulator:
% 'REGULATOR':this word will be replaced with the corresponding regulator
% name (e.g.: regulator: 'sgRNA', 'dCas:REGULATOR' => 'dCas:sgRNA')
% 'REGULATED': this word will be replaced with the regulated object name
% (e.g.: regulated object name: 'node1', 'DNA_REGULATED' => 'DNA_node1')
% OUTPUT will be renamed to the approprialte 'INPUT' ('HILL') parameter 
% (according to the regulated object, e.g.: 'OUTPUT' => 'HILL_node1')

switch obj_type
    case 'node' % node model

        % node type names
        types = {'type1'};
        type = types{n_type};

        switch type
            case 'type1'

                % "position", where the node can be regulated
                addparameter(Mobj,'HILL',1,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

                % === Species ===
                % the last species will be used as regulator for other nodes
                addspecies(Mobj,'mRNA','Units','molecule','Notes','Individual');
                addspecies(Mobj,'uP','Units','molecule','Notes','Individual');
                addspecies(Mobj,'P','Units','molecule','Notes','Individual','Tag','REGULATOR');

                % === Reactions ===
                % mRNA
                % transcription and mRNA degradation
                addreaction(Mobj,'null <-> mRNA','Reactionrate',...
                    'n_copy*(a0+a1*HILL) - (k_mRNA_degr+dilution)*mRNA','Notes','Individual');

                % Protein
                % translation and protein degradation
                % unfolded protein translation, degradation
                addreaction(Mobj,'null <-> uP','Reactionrate',...
                    'k_translation*mRNA - (dilution)*uP','Notes','Individual');

                % maturation
                addreaction(Mobj,'uP -> P','Reactionrate',...
                    'k_mat*uP','Notes','Individual');

                % folded protein degradation
                addreaction(Mobj,'null <-> P','Reactionrate',...
                    '-(dilution)*P','Notes','Individual');

                % === Parameters ===
                addparameter(Mobj,'n_copy',25,'Units','molecule','Notes','Individual');
                addparameter(Mobj,'a0',1e-3,'Units','1/minute','Notes','Individual');
                addparameter(Mobj,'a1',100,'Units','1/minute','Notes','Individual');
                addparameter(Mobj,'k_mRNA_degr',5e-1,'Units','1/minute','Notes','Common');
                addparameter(Mobj,'k_translation',50,'Units','1/minute','Notes','Individual');
                addparameter(Mobj,'k_mat',.4,'Units','1/minute','Notes','Common');

                addparameter(Mobj,'dilution','Value',0.01,'Unit','1/minute','Notes','Common');
        end

    case 'protease' % protease model

        % protease type names
        types = {'type1'};
        type = types{n_type};

        switch type
            case 'type1'

                % the position when the intearaction might be regulated
                addparameter(Mobj,'HILL',1,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

                % === Species ===
                addspecies(Mobj,'uP','Units','molecule','Notes','Individual');
                addspecies(Mobj,'P','Units','molecule','Notes','Individual');
                addspecies(Mobj,'PROT','Units','molecule','InitialAmount',100,'Constant',true,'Notes','Individual','Tag','PROTEASE');

                % === Reactions ===
                % Protein degradation
                % degradation with protease
                addreaction(Mobj,'uP -> null','Reactionrate',...
                    'protease_rate*uP','Notes','Individual');

                % folded protein degradation
                addreaction(Mobj,'P -> null','Reactionrate',...
                    'protease_rate*P','Notes','Individual');

                % === Rules ===
                % rule for the protease
                addrule(Mobj,'protease_rate = k_protease_max*PROT*HILL/(K_protease+Substrates)',...
                    'RuleType','repeatedAssignment','Name','protease_rule','Notes','Individual');
                % calculating the substrate concentraitons (the substrates will be added
                % later
                addrule(Mobj,'Substrates = 0','RuleType','repeatedAssignment',...
                    'Name','substrate_rule','Notes','Individual');

                % === Parameters ===
                % michaelis menten parameters (constant parameters)
                addparameter(Mobj,'K_protease',30,'Units','molecule','Notes','Common'); % MM parameter, half concentration
                addparameter(Mobj,'k_protease_max',50,'Units','1/minute','Notes','Common'); % max speed

                % for the rules (nonconstant parameters)
                addparameter(Mobj,'protease_rate',0,'Units','1/minute','Constant',false,'Notes','Individual');
                addparameter(Mobj,'Substrates',0,'Units','molecule','Constant',false,'Notes','Individual');
        end

    case 'regulator' % regulator model

        % the parameter determined by the Hill function
        addparameter(Mobj,'OUTPUT','Units','dimensionless','Constant',false,'Notes','Individual');

        % regulator types: contains 'A': activation, contains 'R': repression
        types = {'Activation_in','Repression_in','Activation_out','Repression_out'};

        type = types{n_type};

        switch type
            case {'Activation_in','Repression_in'}
                unit = 'molecule';
            case {'Activation_out','Repression_out'}
                unit = 'micromolarity';
        end

        switch unit
            case 'molecule'

                % add 'Hill' parameter
                addparameter(Mobj,'HILL','Value',1,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

                % === Parameters ===
                addparameter(Mobj,['K_' unit],'Value',5,'Units',unit,'Notes','Individual');
                addparameter(Mobj,['n_' unit],'Value',2,'Units','dimensionless','Notes','Individual');

                % === Species ===
                % 'REGULATOR' species (the regulator)
                addspecies(Mobj,'REGULATOR','InitialAmount',0,'Units',unit,'Notes','Individual','Constant',false);

            case 'micromolarity'

                % add 'Hill' parameter
                addparameter(Mobj,'HILL','Value',1,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

                % === Parameters ===
                addparameter(Mobj,['K_' unit],'Value',50,'Units',unit,'Notes','Individual');
                addparameter(Mobj,['n_' unit],'Value',2,'Units','dimensionless','Notes','Individual');

                % === Species ===
                % 'REGULATOR' species (the regulator)
                addspecies(Mobj,'REGULATOR','InitialAmount',0,'Units',unit,'Notes','Individual','Constant',true);
        end

        % rule for the regulation
        % HILL can be replaced later if the regulator has a regulator
        % OUTPUT will be renamed to the approprialte 'HILL' parameter
        % 'REGULATOR' will be renamed with the name of the regulator
        if contains(type,'A') % activation
            rule = ['(REGULATOR*HILL/K_' unit ')^n_' unit '/(1+(REGULATOR*HILL/K_' unit ')^n_' unit ')'];
        elseif contains(type,'R') % repression
            rule = ['1/(1+(REGULATOR*HILL/K_' unit ')^n_' unit ')'];
        end
        addrule(Mobj,['OUTPUT = ' rule],'RepeatedAssignment','Name','Rule_HILL','Notes','Individual');

end

end

