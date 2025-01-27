function [Mobj,types] = CRISPR(Mobj,obj_type,n_type)
% Mobj: the simbiology model object which will contain our model
% obj_type: type of the object ('node','protease','regulator')
% n_type: determines the type of the given object.
% types: define the types of the model types, n_type is the number selects
% the actual type from this cell array

% Javier: https://doi.org/10.1038/s41467-023-38033-3
% "assuming that one copy is 1nM"
% Elowitz: https://doi.org/10.1038/35002125

dilution = 0.018; % dilution rate
ratio = 0.5*60/30/0.044; % Elowitz k1, Javier a1: k1/30/a1

switch obj_type
    case 'node' % node model

        % node type names
        types = {'type1','mixed','Elowitz'};
        type = types{n_type};

        % "position", where the node can be regulated
        if ~strcmp(type,'Elowitz')
            addparameter(Mobj,'NOHILL',1,'Units','dimensionless','Constant',false,'Notes','Common','Tag','INPUT');
        end
        addparameter(Mobj,'HILL',1,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

        % === Species ===
        if ~strcmp(type,'Elowitz')
            addspecies(Mobj,'dCas','Units','molecule','Notes','Common','InitialAmount',1434); % Javier
            addspecies(Mobj,'sgRNA','Units','molecule','Notes','Individual','Tag','REGULATOR');
            addspecies(Mobj,'dCas:sgRNA','Units','molecule','Notes','Individual');
        end
        addspecies(Mobj,'DNA','Units','molecule','Notes','Individual','InitialAmount',30);
        addspecies(Mobj,'mRNA','Units','molecule','Notes','Individual');
        addspecies(Mobj,'P','Units','molecule','Notes','Individual');

        switch type
            case {'mixed','Elowitz'} % 'P' is considered as an output
                set(sbioselect(Mobj.Species,'Name','P'),'Tag','REGULATOR')
        end

        % === Reactions ===
        % mRNA production, degradation
        addreaction(Mobj,'null <-> mRNA','ReactionRate',...
            'a0 + a1*(leak+(1-leak)*HILL)*[DNA] - (dilution+d_RNA)*[mRNA]','Notes','Individual');

        % inactive protein production, degradation
        addreaction(Mobj,'null <-> P','ReactionRate',...
            'k_P*[mRNA] - (dilution+d_P)*[P]','Notes','Individual');

        if ~strcmp(type,'Elowitz')
            % sgRNA production, degradation
            addreaction(Mobj,'null <-> sgRNA','ReactionRate',...
                'a0 + a1*(leak+(1-leak)*HILL)*[DNA] - (dilution+d_RNA)*[sgRNA]','Notes','Individual');

            % dCas:sgRNA production, dilution
            addreaction(Mobj,'dCas + sgRNA <-> dCas:sgRNA','ReactionRate',...
                'kfds*[dCas]*[sgRNA] - krds*[dCas:sgRNA]','Notes','Individual');
            % dCas:sgRNA dilution (sgRNA elimination!)
            addreaction(Mobj,'dCas:sgRNA -> dCas','ReactionRate',...
                'dilution*[dCas:sgRNA]','Notes','Individual');
        end

        % === Parameters ===
        addparameter(Mobj,'dilution','Value',dilution,'Unit','1/minute','Notes','Common');
        addparameter(Mobj,'k_P',20*log(2)/120*60,'Units','1/minute','Notes','Common'); % Elowitz k3
        addparameter(Mobj,'d_P',log(2)/600*60-dilution,'Units','1/minute','Notes','Common'); % Elowitz k4 - dilution
        addparameter(Mobj,'a0',5e-4*60,'Units','molecule/minute','Notes','Common'); % Elowitz k0
        addparameter(Mobj,'a1',0.5*60/30,'Units','1/minute','Notes','Individual'); % Elowitz k1/30
        addparameter(Mobj,'d_RNA',log(2)/120*60-dilution,'Units','1/minute','Notes','Common'); % Elowitz k2 - dilution
        if ~strcmp(type,'Elowitz')
            addparameter(Mobj,'kfds',10^-1.19*ratio,'Units','1/(molecule*minute)','Notes','Common'); % Javier kfSGCAS*ratio
            addparameter(Mobj,'krds',10^-1.11,'Units','1/minute','Notes','Common'); % Javier kbSGCAS
        end

        addparameter(Mobj,'leak',0,'Units','dimensionless','Notes','Common'); % ratio of the leakage

    case 'protease' % protease model

        % empty model
        types = [];

    case 'regulator' % regulator model

        % regulator types: contains 'A': activation, contains 'R': repression
        types = {'Repression_in','Activation_in_TF','Repression_in_TF','Activation_out','Repression_out','Light_Ara','Act_direct'};

        type = types{n_type};

        switch type
            case {'Repression_in','Activation_in_TF','Repression_in_TF'}
                unit = 'molecule';
            case {'Activation_out','Repression_out'}
                unit = 'micromolarity';
            case {'Light_Ara'}
                unit = 'micromolarity';
            case {'Act_direct'}
                unit = 'micromolarity';
        end

        switch type % regulation with Hill function
            case {'Activation_in_TF','Repression_in_TF','Activation_out','Repression_out'}

                % === Parameter ===
                % the parameter determined by the Hill function
                addparameter(Mobj,'OUTPUT','Units','dimensionless','Constant',false,'Notes','Individual','Tag','REGULATOR');

                % add 'Hill' parameter
                addparameter(Mobj,'HILL','Value',0,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

                % === Species ===
                % 'REGULATOR' species (the regulator)
                addspecies(Mobj,'REGULATOR','InitialAmount',0,'Units',unit,'Notes','Individual','Constant',true);

                % === Rule ===
                % rule for the regulation
                % HILL can be replaced later if the regulator has a regulator
                % OUTPUT will be renamed to the appropriate 'HILL' parameter
                % 'REGULATOR' will be replaced with the name of the regulator
                if contains(types{n_type},'A') % activation
                    rule = ['(REGULATOR*HILL/K_' unit ')^n_' unit '/(1+(REGULATOR*HILL/K_' unit ')^n_' unit ')'];
                elseif contains(types{n_type},'R') % repression
                    rule = ['1/(1+(REGULATOR*HILL/K_' unit ')^n_' unit ')'];
                end

                addrule(Mobj,['OUTPUT = ' rule],'RepeatedAssignment','Name','Rule_HILL','Notes','Individual');

            case 'Light_Ara'

                % === Parameter ===
                % the parameter determined by the Hill function
                addparameter(Mobj,'OUTPUT','Units','dimensionless','Constant',false,'Notes','Individual','Tag','REGULATOR');

                % add 'Hill' parameter
                addparameter(Mobj,'HILL_light','Value',0,'Units',unit,'Constant',false,'Notes','Individual','Tag','INPUT');
                addparameter(Mobj,'HILL_ara','Value',0,'Units',unit,'Constant',false,'Notes','Individual','Tag','INPUT');

                % === Species ===
                % we will not have a 'REGULATOR' species
                % just the Hill functions from the light and ara system

                % 'REGULATOR' species (the regulator)
                addspecies(Mobj,'REGULATOR','InitialAmount',0,'Units',unit,'Notes','Individual','Constant',true);

                % === Rule ===
                % rule for the regulation
                % HILL_light and Hill_ara must be replaced by the
                % corresponding regulators
                % OUTPUT will be renamed to the appropriate 'HILL' parameter
                rule = '((HILL_ara/K_ara)^n_ara+k_light_ratio*(HILL_light/K_light)^n_light)/((1+(HILL_ara/K_ara)^n_ara)*(1+(HILL_light/K_light)^n_light))';

                addrule(Mobj,['OUTPUT = ' rule],'RepeatedAssignment','Name','Rule_HILL','Notes','Individual');

            case 'Act_direct' % just replace HILL with the regulator

                % === Parameter ===
                % the parameter determined by the Hill function
                addparameter(Mobj,'OUTPUT','Units',unit,'Constant',false,'Notes','Individual','Tag','REGULATOR');

                % add 'Hill' parameter
                addparameter(Mobj,'HILL','Value',0,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

                % === Species ===
                % 'REGULATOR' species (the regulator)
                addspecies(Mobj,'REGULATOR','InitialAmount',0,'Units',unit,'Notes','Individual','Constant',true);

                % === Rule ===
                % rule for the regulation
                % HILL_light and Hill_ara must be replaced by the
                % corresponding regulators
                % OUTPUT will be renamed to the appropriate 'HILL' parameter
                rule = 'REGULATOR*HILL';

                addrule(Mobj,['OUTPUT = ' rule],'RepeatedAssignment','Name','Rule_HILL','Notes','Individual');

        end

        switch type
            case 'Repression_in' % regulation through CRISPR system

                % === Species ===
                % 'REGULATOR' species (the regulator: sgRNA)
                addspecies(Mobj,'DNA_REGULATED','Units','molecule','Notes','Individual','InitialAmount',0);
                addspecies(Mobj,'dCas','Units','molecule','Notes','Common','InitialAmount',0);
                addspecies(Mobj,'dCas:REGULATOR','InitialAmount',0,'Units',unit,'Notes','Individual','Constant',false);
                addspecies(Mobj,'dCas:REGULATOR:DNA_REGULATED','InitialAmount',0,'Units',unit,'Notes','Individual','Constant',false);

                % dCas:sgRNA:DNA production, degradation
                addreaction(Mobj,'dCas:REGULATOR + DNA_REGULATED <-> dCas:REGULATOR:DNA_REGULATED','ReactionRate',...
                    'kfdsd*[dCas:REGULATOR]*[DNA_REGULATED] - krdsd*[dCas:REGULATOR:DNA_REGULATED]','Notes','Individual');
                % dCas:sgRNA:DNA dilution (sgRNA elimination!)
                addreaction(Mobj,'dCas:REGULATOR:DNA_REGULATED -> dCas + DNA_REGULATED','ReactionRate',...
                    'dilution*[dCas:REGULATOR:DNA_REGULATED]','Notes','Individual');
                
                % === Parameters ===
                addparameter(Mobj,'kfdsd',10^-1.93*ratio,'Units','1/(molecule*minute)','Notes','Common'); % Javier kfGCSGsg2*ratio
                addparameter(Mobj,'krdsd',0,'Units','1/minute','Notes','Common'); % Javier oscillator: 0, kbGCSGsg1: 10^(-1.1e-1)

            case {'Activation_in_TF','Repression_in_TF'} % regulation through transcription factors

                % === Parameters ===
                addparameter(Mobj,['K_' unit],'Value',40,'Units','molecule','Notes','Individual');
                addparameter(Mobj,['n_' unit],'Value',2,'Units','dimensionless','Notes','Individual');

            case {'Activation_out','Repression_out'} % regulation through outer species (arabinose)

                % === Parameters ===
                addparameter(Mobj,['K_' unit],'Value',38.9045,'Units',unit,'Notes','Individual'); % Javier Km
                addparameter(Mobj,['n_' unit],'Value',1.0233,'Units','dimensionless','Notes','Individual'); % Javier ARAn

            case 'Light_Ara'

                % === Parameters ===
                addparameter(Mobj,'k_light_ratio','Value',0.31854,'Units','dimensionless','Notes','Individual');
                addparameter(Mobj,'K_light','Value',156.8606,'Units',unit,'Notes','Individual');
                addparameter(Mobj,'n_light','Value',0.98973,'Units','dimensionless','Notes','Individual');
                addparameter(Mobj,'K_ara','Value',0.00027567,'Units',unit,'Notes','Individual');
                addparameter(Mobj,'n_ara','Value',1.1177,'Units','dimensionless','Notes','Individual');


        end

end

end


