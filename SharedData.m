classdef SharedData < matlab.mixin.Copyable % copyable handle class
    % common value of the nodes
    % Mobj: the simbiology model obect
    properties 

      % name/type of the reacdtion kinetics model
      model_name (1,:) char;

      % model of the whole cell
      Mobj SimBiology.Model;

      % additional model what we can load separately and we can add it to
      % our model when we export the whole model
      Mobj_extern SimBiology.Model;

      % struct for the different node models
      node_models struct;

      % struct for the different protease models
      protease_models struct;

      % struct for the different regulator models
      regulator_models struct;

      % common  properties
      name (1,:) char = 'GRN';
      compartment_name (1,:) char = 'Ecoli';
      volume SimBiology.Parameter;

      % non redundant struct of the regulator objects (handle object)
      regulators struct;
      % route fot the regulators (useful in deleting regulators)
      regulator_routes struct;

      %% Simulation settings
      StatesToLog (1,:) cell = {};
      Accelerate (1,1) logical = false;

      %% Graph properties

      G; % the graph of the system

      original_graph (1,1) logical = false;      % to plot only the original graph

      node_size (1,1) double = 15;               % MarkerSize for the node
      inner_regulator_size (1,1) double = 3;     % MarkerSize for the regulator produced by a node
      outer_regulator_size (1,1) double = 10;    % MarkerSize for the not node produced regulator
      protease_size (1,1) double = 10;           % MarkerSize for the protease
      linewidth_protease (1,1) double = 1.5;     % protease edge width
      linewidth_regulator(1,1) double  = 2;      % regulator edge width
      node_color (1,:) char = 'red';            % node color
      protease_color (1,:) char = 'blue';       % protease color
      regulation_color (1,:) char = 'yellow';   % regulation color
      regulator_color (1,:) char = 'green';     % regulator color
      node_fontsize (1,1) double = 8;           % FontSize of the nodes
      node_fontweight (1,:) char = 'bold';      % FontWeight of the nodes

      act_edge_color (1,:) char = '#77AC30';    % activation edge color
      repr_edge_color (1,:) char = '#A2142F';   % repression edge color
      other_edge_color (1,:) char = '#7E2F8E';  % other edge type color
      protease_edge_color (1,:) char = '#B2B2B2'; % protease color

      layout_type (1,:) char = 'force';%'force';%'layered'; % layout type of the graph

    end

    methods

        % constructor
        function obj = SharedData(loadgraph)

            obj.clean_start();
            obj.set_Mobj();

            % load graph settings 
            if nargin==0 || loadgraph==true
                obj.set_graph();
            end
        end
    
        % after sbioreset the simbiology objects will be deleted
        % set up them again
        function set_Mobj(obj)
            obj.volume = addparameter(sbiomodel('model'),'volume','Value',0.7,'Unit','micrometer^3','Notes','Common');
        end

        % set changable properties to initial state
        function clean_start(obj)
            % delete the regulator lists
            obj.regulators = struct;
            obj.regulator_routes = struct;
            % delete the object model lists
            obj.node_models = struct;
            obj.protease_models = struct;
            obj.regulator_models = struct;

            obj.StatesToLog = {};
            obj.Accelerate = false;
        end

        % set every properties according to another ShardData object
        function set_properties(obj,other)
            properties = fieldnames(obj);
            for i = 1:numel(properties)
                obj.(properties{i}) = other.(properties{i});
            end
        end

        % set graph properties acccording to the saved settings
        function set_graph(obj)
            % find the correct folder
            [folderPath, ~, ~] = fileparts(which('SharedData'));
            if ispc % Win
                folderPath  = [folderPath '\app\'];
            else % linux, mac
                folderPath  = [folderPath '/app/'];
            end
            % check if we have the setting file
            if exist([folderPath 'graph_settings.mat'],'file')==2
                % get the settings saved in the file
                settings = load([folderPath 'graph_settings.mat']);
                names = fieldnames(settings);
                for i = 1:numel(names)
                    obj.(names{i}) = settings.(names{i});
                end
            end
        end

        % save and apply default graph settings
        function default_graph_settings(obj)

            % create a new instance with the default graph settings
            default_settings = SharedData(false);
            % copy the properties into a structure
            properties = fieldnames(default_settings);
            for i = 1:numel(properties)
                s.(properties{i}) = default_settings.(properties{i});
            end
            % find the correct folder
            [folderPath, ~, ~] = fileparts(which('SharedData'));
            if ispc % Win
                folderPath  = [folderPath '\app\'];
            else
                folderPath  = [folderPath '/app/'];
            end % linux, mac
            % save the settings
            save([folderPath 'graph_settings.mat'],'-struct','s','node_size','inner_regulator_size', ...
                'outer_regulator_size','protease_size','linewidth_protease', ...
                'linewidth_regulator','node_color','protease_color', ...
                'regulation_color','regulator_color','node_fontsize', ...
                'node_fontweight','act_edge_color','repr_edge_color', ...
                'other_edge_color','protease_edge_color','layout_type','original_graph');
            % apply the sLettings
            obj.set_graph();

        end

    end
    
end

