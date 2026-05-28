function [Mobj,types] = CRISPR_Murray(Mobj,obj_type,n_type)
% Mobj: the simbiology model object which will contain our model
% obj_type: type of the object ('node','protease','regulator')
% n_type: determines the type of the given object.
% types: define the types of the model types, n_type is the number selects
% the actual type from this cell array

% Murray: https://www.biorxiv.org/content/10.1101/225318v2
% dummy: https://www.biorxiv.org/content/10.1101/2021.09.26.461880v1
% Murray code: https://github.com/sclamons/plasmid_replication_modeling/blob/main/code/CRISPRessilator.ipynb
% "assuming that one copy is 1nM"

% dilution = 0.018; % dilution rate
% ratio = 0.5*60/30/0.044; % Elowitz k1, Javier a1: k1/30/a1

% for stochastic simulations we need dummy particles
STOCH_DUMMY = false;

% All parameters in units of 1/minute or 1/molecule/minute (1/nM/minute)
k_gamma       = log(2) / 30; % dilution
k_gRNA_prod   = 5.0;
k_gRNA_deg    = log(2) / 100 * 60;
k_dCas_prod   = 4.5 * 1.0;
k_dCas_deg    = 0 * k_gamma;
k_gRNA_bind   = log(2) / 375 *60;
k_gRNA_unbind = 0;
k_complex_bind     = log(2);
k_complex_unbind   = 0; % k_complex_bind / 1000.0 <-- Just a guess!
k_gRNA_leak   = 0; % gRNA_prod_rate / 100.0

% Murry introduced a dummy species, R, to describe stochastic DNA
% production, dilution properly:
% https://www.biorxiv.org/content/10.1101/2021.09.26.461880v1
DNA_eq  = 30; % equilibrium DNA concentration / molecule
if STOCH_DUMMY==true
    k_alpha = DNA_eq * k_gamma; % dummy molecule production
    k_dummy = 10 * k_gamma; % dummy molecule consumption
end


switch obj_type
    case 'node' % node model

        % node type names
        types = {'type1'};
        type = types{n_type};

        % "position", where the node can be regulated
        if ~strcmp(type,'Elowitz')
            addparameter(Mobj,'NOHILL',1,'Units','dimensionless','Constant',false,'Notes','Common','Tag','INPUT');
        end
        addparameter(Mobj,'HILL',1,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

        % === Species ===
        if ~strcmp(type,'Elowitz')
            addspecies(Mobj,'dCas','Units','molecule','Notes','Common','InitialAmount',k_dCas_prod/(k_dCas_deg+k_gamma));
            addspecies(Mobj,'sgRNA','Units','molecule','Notes','Individual','Tag','REGULATOR');
            addspecies(Mobj,'dCas:sgRNA','Units','molecule','Notes','Individual');
            if STOCH_DUMMY==true
                addspecies(Mobj,'Dummy','Units','molecule','Notes','Individual'); % dummy particle for DNA replication
            end
        end
        addspecies(Mobj,'DNA','Units','molecule','Notes','Individual','InitialAmount',DNA_eq);
%         addspecies(Mobj,'mRNA','Units','molecule','Notes','Individual');
%         addspecies(Mobj,'P','Units','molecule','Notes','Individual');

%         switch type
%             case {'mixed','Elowitz'} % 'P' is considered as an output
%                 set(sbioselect(Mobj.Species,'Name','P'),'Tag','REGULATOR')
%         end

        % === Reactions ===
%         % mRNA production, degradation
%         addreaction(Mobj,'null <-> mRNA','ReactionRate',...
%             'a0 + a1*(leak+(1-leak)*HILL)*[DNA] - (dilution+d_RNA)*[mRNA]','Notes','Individual');

%         % inactive protein production, degradation
%         addreaction(Mobj,'null <-> P','ReactionRate',...
%             'k_P*[mRNA] - (dilution+d_P)*[P]','Notes','Individual');

        if ~strcmp(type,'Elowitz')
            % sgRNA production, degradation
            addreaction(Mobj,'null <-> sgRNA','ReactionRate',...
                'k_gRNA_prod*HILL*[DNA] - (k_gamma+k_gRNA_deg)*[sgRNA]','Notes','Individual');

            % dCas:sgRNA production, dilution
            addreaction(Mobj,'dCas + sgRNA <-> dCas:sgRNA','ReactionRate',...
                'k_gRNA_bind*[dCas]*[sgRNA] - k_gRNA_unbind*[dCas:sgRNA]','Notes','Individual');
            % dCas:sgRNA dilution (sgRNA elimination!)
            addreaction(Mobj,'dCas:sgRNA -> null','ReactionRate',...
                '(k_dCas_deg+k_gamma)*[dCas:sgRNA]','Notes','Individual');

            if STOCH_DUMMY==true
                % dummy species production
                addreaction(Mobj,'null -> Dummy','ReactionRate',...
                    'k_alpha','Notes','Individual');
                % DNA replication with dummy species
                addreaction(Mobj,'DNA + Dummy -> 2 DNA','ReactionRate',...
                    'k_dummy*[DNA]*[Dummy]','Notes','Individual');
                % DNA dilution
                addreaction(Mobj,'DNA -> null','ReactionRate',...
                    'k_gamma*[DNA]','Notes','Individual');
            end

            % dCas production, degradation
            addreaction(Mobj,'null <-> dCas','ReactionRate',...
                'k_dCas_prod - (k_dCas_deg+k_gamma)*[dCas]','Notes','Common');

        end

        % === Parameters ===
        addparameter(Mobj,'k_gamma','Value',k_gamma,'Unit','1/minute','Notes','Common');
%         addparameter(Mobj,'k_P',20*log(2)/120*60,'Units','1/minute','Notes','Common'); % Elowitz k3
%         addparameter(Mobj,'d_P',log(2)/600*60-dilution,'Units','1/minute','Notes','Common'); % Elowitz k4 - dilution
%         addparameter(Mobj,'a0',5e-4*60,'Units','molecule/minute','Notes','Common'); % Elowitz k0
%         addparameter(Mobj,'a1',0.5*60/30,'Units','1/minute','Notes','Individual'); % Elowitz k1/30
%         addparameter(Mobj,'d_RNA',log(2)/120*60-dilution,'Units','1/minute','Notes','Common'); % Elowitz k2 - dilution
        if ~strcmp(type,'Elowitz')
            addparameter(Mobj,'k_gRNA_prod',k_gRNA_prod,'Units','1/minute','Notes','Common');
            addparameter(Mobj,'k_gRNA_deg',k_gRNA_deg,'Units','1/minute','Notes','Common');
            addparameter(Mobj,'k_gRNA_bind',k_gRNA_bind,'Units','1/(molecule*minute)','Notes','Common');
            addparameter(Mobj,'k_gRNA_unbind',k_gRNA_unbind,'Units','1/minute','Notes','Common');
            if STOCH_DUMMY==true
                addparameter(Mobj,'k_alpha',k_alpha,'Units','molecule/minute','Notes','Common');
                addparameter(Mobj,'k_dummy',k_dummy,'Units','1/(molecule*minute)','Notes','Common');
            end
            addparameter(Mobj,'k_dCas_prod',k_dCas_prod,'Units','molecule/minute','Notes','Common');
            addparameter(Mobj,'k_dCas_deg',k_dCas_deg,'Units','1/minute','Notes','Common');
        end

    case 'protease' % protease model

        % empty model
        types = [];

    case 'regulator' % regulator model

        % regulator types: contains 'A': activation, contains 'R': repression
%         types = {'Repression_in','Activation_in_TF','Repression_in_TF','Activation_out','Repression_out','Light_Ara','Act_direct'};
        types = {'Repression_in','Activation_out','Repression_out','Light_Ara','Act_direct'};

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
                addspecies(Mobj,'sgRNA_REGULATED','Units','molecule','Notes','Individual','InitialAmount',0);
                if STOCH_DUMMY==true
                    addspecies(Mobj,'Dummy_REGULATED','Units','molecule','Notes','Individual','InitialAmount',0);
                end

                % dCas:sgRNA:DNA production, degradation
                addreaction(Mobj,'dCas:REGULATOR + DNA_REGULATED <-> dCas:REGULATOR:DNA_REGULATED','ReactionRate',...
                    'k_complex_bind*[dCas:REGULATOR]*[DNA_REGULATED] - k_complex_unbind*[dCas:REGULATOR:DNA_REGULATED]','Notes','Individual');
                
                % leakage for sgRNA production
                addreaction(Mobj,'dCas:REGULATOR:DNA_REGULATED -> dCas:REGULATOR:DNA_REGULATED + sgRNA_REGULATED','ReactionRate',...
                    'k_gRNA_leak*[dCas:REGULATOR:DNA_REGULATED]','Notes','Individual');

                if STOCH_DUMMY==true
                    % DNA replication with dummy species
                    addreaction(Mobj,'dCas:REGULATOR:DNA_REGULATED + Dummy_REGULATED -> dCas:REGULATOR + 2 DNA_REGULATED','ReactionRate',...
                        'k_dummy*[dCas:REGULATOR:DNA_REGULATED]*[Dummy_REGULATED]','Notes','Individual');
                    % dCas:sgRNA:DNA dilution (sgRNA elimination!)
                    addreaction(Mobj,'dCas:REGULATOR:DNA_REGULATED -> null','ReactionRate',...
                        'k_gamma*[dCas:REGULATOR:DNA_REGULATED]','Notes','Individual');
                else
                    % dCas:sgRNA:DNA dilution (sgRNA elimination!) DNA_REGULATED
                    addreaction(Mobj,'dCas:REGULATOR:DNA_REGULATED -> DNA_REGULATED','ReactionRate',...
                        'k_gamma*[dCas:REGULATOR:DNA_REGULATED]','Notes','Individual');
                end

                
                % === Parameters ===
                addparameter(Mobj,'k_complex_bind',k_complex_bind,'Units','1/(molecule*minute)','Notes','Common');
                addparameter(Mobj,'k_complex_unbind',k_complex_unbind,'Units','1/minute','Notes','Common');
                addparameter(Mobj,'k_gRNA_leak',k_gRNA_leak,'Units','1/minute','Notes','Common');
                if STOCH_DUMMY==true
                    addparameter(Mobj,'k_dummy',k_dummy,'Units','1/(molecule*minute)','Notes','Common');
                end

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


