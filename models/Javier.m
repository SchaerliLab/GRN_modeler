function [Mobj,types] = Javier(Mobj,obj_type,n_type)
% Mobj: the simbiology model object which will contain our model
% obj_type: type of the object ('node','protease','regulator')
% n_type: determines the type of the given object.
% types: define the types of the model types, n_type is the number selects
% the actual type from this cell array

% https://doi.org/10.1038/s41467-023-38033-3:
% "assuming that one copy is 1nM"

switch obj_type
    case 'node' % node model

        % node type names
        types = {'type1'};%{'A','B','C'};
        % type = types{n_type};


                % "position", where the node can be regulated
                addparameter(Mobj,'NOHILL',1,'Units','dimensionless','Constant',false,'Notes','Common','Tag','INPUT');
                addparameter(Mobj,'HILL',1,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

                % === Species ===
                addspecies(Mobj,'dCas','Units','molecule','Notes','Common','InitialAmount',1434.0390711255027); % 2.1677e3

                addspecies(Mobj,'DNA','Units','molecule','Notes','Individual','InitialAmount',30);
                addspecies(Mobj,'mRNA','Units','molecule','Notes','Individual');
                addspecies(Mobj,'PI','Units','molecule','Notes','Individual');
                addspecies(Mobj,'P','Units','molecule','Notes','Individual');
                addspecies(Mobj,'sgRNA','Units','molecule','Notes','Individual','Tag','REGULATOR');
                addspecies(Mobj,'dCas:sgRNA','Units','molecule','Notes','Individual');
                % addspecies(Mobj,'dCas:sgRNA:DNA','Units','molecule','Notes','Individual');

                % === Reactions ===
                % mRNA production, degradation
                addreaction(Mobj,'null <-> mRNA','ReactionRate',...
                    'a0_mRNA + a1*HILL*[DNA] - d_RNA*[mRNA]','Notes','Individual');

                % inactive protein production, degradation
                addreaction(Mobj,'null <-> PI','ReactionRate',...
                    'k_PI*[mRNA] - (dilution+d_PI)*[PI]','Notes','Individual');

                % inactive-active protein transformation
                addreaction(Mobj,'PI -> P','ReactionRate',...
                    'm_PI*[PI]','Notes','Individual');

                % active protein degradation
                addreaction(Mobj,'P -> null','ReactionRate',...
                    '(dilution+d_PI)*[P]','Notes','Individual');

                % sgRNA production, degradation
                addreaction(Mobj,'null <-> sgRNA','ReactionRate',...
                    'a0_sgRNA + a1*HILL*[DNA] - d_RNA*[sgRNA]','Notes','Individual');

                % dCas:sgRNA production, dilution
                addreaction(Mobj,'dCas + sgRNA <-> dCas:sgRNA','ReactionRate',...
                    'kfds*[dCas]*[sgRNA] - krds*[dCas:sgRNA]','Notes','Individual');
                % dCas:sgRNA dilution (sgRNA elimination!)
                addreaction(Mobj,'dCas:sgRNA -> dCas','ReactionRate',...
                    'dilution*[dCas:sgRNA]','Notes','Individual');


% 
%                         % rename sgRNA to sgRNA1
%                         rename(sbioselect(Mobj.Species,'Name','sgRNA'),'sgRNA1');
%                         rename(sbioselect(Mobj.Species,'Name','dCas:sgRNA'),'dCas:sgRNA1');
%                         rename(sbioselect(Mobj.Species,'Name','dCas:sgRNA:DNA'),'dCas:sgRNA1:DNA');
% 
%                         % sgRNA production, degradation
%                         addreaction(Mobj,'null <-> sgRNA2','ReactionRate',...
%                             'a0_sgRNA2 + a1*HILL*[DNA] - d_RNA*[sgRNA2]','Notes','Individual');
% 
%                         % dCas:sgRNA production, dilution
%                         addreaction(Mobj,'dCas + sgRNA2 <-> dCas:sgRNA2','ReactionRate',...
%                             'kfds*[dCas]*[sgRNA2] - krds*[dCas:sgRNA2]','Notes','Individual');
%                         % dCas:sgRNA dilution (sgRNA elimination!)
%                         addreaction(Mobj,'dCas:sgRNA2 -> dCas','ReactionRate',...
%                             'dilution*[dCas:sgRNA2]','Notes','Individual');


                        % === Parameters ===
                        addparameter(Mobj,'k_PI',10^-1.46,'Units','1/minute','Notes','Common'); %  kMaPIa
                        addparameter(Mobj,'d_PI',10^-3.58176608841291,'Units','1/minute','Notes','Common'); %  dPa
                        addparameter(Mobj,'m_PI',10^-1.9421824463322,'Units','1/minute','Notes','Common'); % kPIaPa

                        % switch type
                        %     case 'A'
                        %         addparameter(Mobj,'k_PI',10^-1.46,'Units','1/minute','Notes','Common'); %  kMaPIa
                        %         addparameter(Mobj,'d_PI',10^-3.58176608841291,'Units','1/minute','Notes','Common'); %  dPa
                        %         addparameter(Mobj,'m_PI',10^-1.9421824463322,'Units','1/minute','Notes','Common'); % kPIaPa
                        %     case 'B'
                        %         addparameter(Mobj,'k_PI',10^-1.24631340961862,'Units','1/minute','Notes','Common'); % kMbPIb
                        %         addparameter(Mobj,'d_PI',10^-4.17145764268855,'Units','1/minute','Notes','Common'); % dPb
                        %         addparameter(Mobj,'m_PI',10^-1.00007653476072,'Units','1/minute','Notes','Common'); % kPIbPb
                        %     case 'C'
                        %         addparameter(Mobj,'k_PI',10^-0.135660911341568,'Units','1/minute','Notes','Common'); % kMcPIc
                        %         addparameter(Mobj,'d_PI',10^-2.15301105071788,'Units','1/minute','Notes','Common'); % dPc
                        %         addparameter(Mobj,'m_PI',10^1.97211569359306,'Units','1/minute','Notes','Common'); % kPIcPc
                        % end

                        addparameter(Mobj,'a0_mRNA',9.7004e-04,'Units','molecule/minute','Notes','Common'); % bA*VMAX*GA, bA:10^-3.9302785791727 VMAX: 10^-0.56 
                        % modified !!!
                        addparameter(Mobj,'a1',0.044,'Units','1/minute','Notes','Common'); % KMC: 10^0.17
                        addparameter(Mobj,'d_RNA',10^0.860098805220079,'Units','1/minute','Notes','Common'); %dmRNAcaseA, dmRNAcaseC: 6.2946
                        addparameter(Mobj,'a0_sgRNA',1.1405e-04,'Units','molecule/minute','Notes','Common');% bSG1*VMAX*GA, bSG1:1.3804e-05
                        addparameter(Mobj,'kfds',10^-1.19,'Units','1/(molecule*minute)','Notes','Common'); % kfSGCAS
                        addparameter(Mobj,'krds',10^-1.11,'Units','1/minute','Notes','Common'); % kbSGCAS


                        % modified !!!
                        addparameter(Mobj,'dilution','Value',0.006,'Unit','1/minute','Notes','Common'); % m: 1e-3, oscillator parameter


    case 'protease' % protease model

        % empty model
        types = [];

    case 'regulator' % regulator model

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
                addparameter(Mobj,'kfdsd',10^-1.93,'Units','1/(molecule*minute)','Notes','Common'); %kfGCSGsg2, kfGCSGsg4: 0.0224, kfGCSGsg1: 0.014454397707459272
                addparameter(Mobj,'krdsd',0*10^-0.41,'Units','1/minute','Notes','Common'); % kbGCSGsg2, kbGCSGsg4:0.7762 (0), kbGCSGsg1:0; oscillator: 0


            case 'micromolarity'

                % the parameter determined by the Hill function
                addparameter(Mobj,'OUTPUT','Units','dimensionless','Constant',false,'Notes','Individual');

                % add 'Hill' parameter
                addparameter(Mobj,'HILL','Value',0,'Units','dimensionless','Constant',false,'Notes','Individual','Tag','INPUT');

                % === Parameters ===
                addparameter(Mobj,['K_' unit],'Value',38.9045,'Units',unit,'Notes','Individual'); % Km
                addparameter(Mobj,['n_' unit],'Value',1.0233,'Units','dimensionless','Notes','Individual'); % ARAn

                % === Species ===
                % 'REGULATOR' species (the regulator)
                addspecies(Mobj,'REGULATOR','InitialAmount',0,'Units',unit,'Notes','Individual','Constant',true);

                % === Rule ===
                % rule for the regulation
                % HILL can be replaced later if the regulator has a regulator
                % OUTPUT will be renamed to the approprialte 'HILL' parameter
                % 'REGULATOR' will be renamed with the name of the regulator
                if contains(types{n_type},'A') % activation
                    rule = ['(REGULATOR*HILL/K_' unit ')^n_' unit '/(1+(REGULATOR*HILL/K_' unit ')^n_' unit ')'];
                elseif contains(types{n_type},'R') % repression
                    rule = ['1/(1+(REGULATOR*HILL/K_' unit ')^n_' unit ')'];
                end
                addrule(Mobj,['OUTPUT = ' rule],'RepeatedAssignment','Name','Rule_HILL','Notes','Individual');
        end

end

end


