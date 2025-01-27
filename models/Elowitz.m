function [Mobj,types] = Elowitz(Mobj,obj_type,n_type)
% Mobj: the simbiology model object which will contain our model
% obj_type: type of the object ('node','protease','regulator')
% n_type: determines the type of the given object.
% types: define the types of the model types, n_type is the number selects
% the actual type from this cell array

% Elowitz: https://doi.org/10.1038/35002125

switch obj_type
    case 'node' % node model

        % node type names
        types = {'type1'};
        type = types{n_type};

        switch type
            case 'type1'

                % "position", where the node can be regulated
                addparameter(Mobj,'HILL',1,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');
                % addparameter(Mobj,'HILL2',1,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

                % === Species ===
                % the last species will be used as regulator for other nodes
                addspecies(Mobj,'mRNA','Units','molecule','Notes','Individual');
                addspecies(Mobj,'P','Units','molecule','Notes','Individual','Tag','REGULATOR');

                % === Reactions ===
                % mRNA
                addreaction(Mobj,'null <-> mRNA','ReactionRate',...
                    'k0 + k1*HILL - k2*[mRNA]','Notes','Individual');

                % Protein
                addreaction(Mobj,'null <-> P','ReactionRate',...
                    'k3*[mRNA] - k4*[P]','Notes','Individual');

                % === Parameters ===
                addparameter(Mobj,'k0',5e-4*60,'Units','molecule/minute','Notes','Individual');
                addparameter(Mobj,'k1',0.5*60,'Units','molecule/minute','Notes','Individual');
                addparameter(Mobj,'k2',log(2)/120*60,'Units','1/minute','Notes','Individual');
                addparameter(Mobj,'k3',20*log(2)/120*60,'Units','1/minute','Notes','Individual');
                addparameter(Mobj,'k4',log(2)/600*60,'Units','1/minute','Notes','Individual');
        end

    case 'protease' % protease model

        % it is not completely empty just to be able to switch between this
        % and the Tomazou model
        % protease type names
        types = {'type1'};
        type = types{n_type};

        switch type
            case 'type'
                % empty model
        end

    case 'regulator' % regulator model

        % the parameter determined by the Hill function
        addparameter(Mobj,'OUTPUT','Units','dimensionless','Constant',false,'Notes','Individual','Tag','REGULATOR');

        % regulator types: contains 'A': activation, contains 'R': repression
        types = {'Activation_in','Repression_in','Mixed_in','Activation_out','Repression_out'};

        type = types{n_type};

        % add 'Hill' parameter
        addparameter(Mobj,'HILL','Value',0,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

        switch type
            case {'Activation_in','Repression_in','Mixed_in'}

                % === Parameters ===
                addparameter(Mobj,'K','Value',40,'Units','molecule','Notes','Individual');
                addparameter(Mobj,'n','Value',2,'Units','dimensionless','Notes','Common');

                % === Species ===
                % 'REGULATOR' species (the regulator)
                addspecies(Mobj,'REGULATOR','InitialAmount',0,'Units','molecule','Notes','Individual','Constant',false);

            case {'Activation_out','Repression_out'}

                % === Parameters ===
                addparameter(Mobj,'K_out','Value',40,'Units','micromolarity','Notes','Individual');
                addparameter(Mobj,'n_out','Value',2,'Units','dimensionless','Notes','Individual');

                % === Species ===
                % 'REGULATOR' species (the regulator)
                addspecies(Mobj,'REGULATOR','InitialAmount',0,'Units','micromolarity','Notes','Individual','Constant',false);
        end

        % rule for the regulation
        % HILL can be replaced later if the regulator has a regulator
        % OUTPUT will be renamed to the approprialte 'HILL' parameter
        % 'REGULATOR' will be renamed with the name of the regulator
        switch type
            case 'Activation_in' % activation
                rule = '(REGULATOR*HILL/K)^n/(1+(REGULATOR*HILL/K)^n)';
            case 'Repression_in' % repression
                rule = '1/(1+(REGULATOR*HILL/K)^n)';
            case 'Mixed_in' % activation*repression
                rule = '(REGULATOR*HILL/K)^n/(1+(REGULATOR*HILL/K)^n).^2';
            case 'Activation_out' % activation
                rule = '(REGULATOR*HILL/K_out)^n_out/(1+(REGULATOR*HILL/K_out)^n_out)';
            case 'Repression_out' % repression
                rule = '1/(1+(REGULATOR*HILL/K_out)^n_out)';
        end
        addrule(Mobj,['OUTPUT = ' rule],'RepeatedAssignment','Name','Rule_HILL','Notes','Individual');

end

end



