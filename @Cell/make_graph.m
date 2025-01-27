function [h] = make_graph(obj,ha)
% createand plot graph of the system
% ha: axes handle

weight.protease = 1;  % protease weight
weight.regulator = 2; % regulator weight

obj.data.G = [];
obj.data.G = digraph();

%% Add nodes and edges to the graph
% add nodes as graph nodes
node_names = fieldnames(obj.nodes);
for i = 1:size(node_names,1)
    obj.data.G = addnode(obj.data.G,node_names{i});
end

% add protease (graph nodes)
protease_names = fieldnames(obj.proteases);
for i = 1:size(protease_names,1)
    obj.data.G = addnode(obj.data.G,protease_names{i});
end

% if there are no nodes and proteases then we can plot nothing
if isempty(node_names) && isempty(protease_names)
    if nargin > 1
        delete(get(ha,'Children'))
    end
    return;
end

% protease reaction (edges)
for i = 1:size(node_names,1)
    % there is a protease reaction, add the edge
    if ~isempty(obj.nodes.(node_names{i}).protease)
        obj.data.G = addedge(obj.data.G,node_names{i},obj.nodes.(node_names{i}).protease,weight.protease);
    end
end

% add regulators (graph nodes) and their interactions (edges)
% positive edge weight: activation
% negative edge weight: inhibition
% regulators for the nodes and the regulators of the regs iteratively
for i = 1:numel(node_names)
    obj.reg2graph(obj.nodes.(node_names{i}),weight);
end
% regulators for the protease and the regulators of the regs iteratively
for i = 1:numel(protease_names)
    obj.reg2graph(obj.proteases.(protease_names{i}),weight);
end

%% Delete unused Graph Nodes
% e.g. a regulator which is not present in rules

% find the names of the inner regulators (produced by a node) 
% and the outer regulators (not produced by nodes)
regulator_names = obj.data.G.Nodes.Name(numel(node_names)+numel(protease_names)+1:end);
inner_regulator_names = false(length(regulator_names),1);
for i = 1:length(regulator_names)
    inner_regulator_names(i) = any(strcmp(regulator_names(i),obj.data.G.Edges.EndNodes(:,2)));
end
outer_regulator_names = ~inner_regulator_names;
inner_regulator_names = regulator_names(inner_regulator_names);
outer_regulator_names = regulator_names(outer_regulator_names);

% find unused species
unused_spec = findUnusedComponents(obj.data.Mobj);
% just the species, parameters are not interesting
unused_spec = unused_spec(strcmp({unused_spec.Type},'species'));
% the name of these species
unused_spec_name = {unused_spec.Name};
% the species which correspond to a node and we want to remove it
remove_spec = false(numel(unused_spec),1);

% if we have a Node in the Graph with this name, delete it
for i = 1:numel(unused_spec)
    if any(strcmp(obj.data.G.Nodes.Name,unused_spec_name{i}))

        % select the node what we want to remove
        remove_spec(i) = true;

        % delete from nodes
        node_names(strcmp(node_names,unused_spec_name{i})) = [];
        % delete from inner regulators
        inner_regulator_names(strcmp(inner_regulator_names,unused_spec_name{i})) = [];
        % delete from outer regulators
        outer_regulator_names(strcmp(outer_regulator_names,unused_spec_name{i})) = [];
        % delete from proteases
        protease_names(strcmp(protease_names,unused_spec_name{i})) = [];
    end
end
% remove the corresponding nodes
obj.data.G = rmnode(obj.data.G,unused_spec_name(remove_spec));

%% plot the graph
% interpreter: none to see the special characters
if nargin > 1
    h = plot(ha,obj.data.G, 'Interpreter', 'none');
else
    h = plot(obj.data.G, 'Interpreter', 'none');
end
disableDefaultInteractivity(get(h,'Parent'))
% change layout
layout(h,obj.data.layout_type) % 'force'

%% set colors and labels
% nodes
highlight(h,node_names,'NodeColor',obj.data.node_color,'MarkerSize',obj.data.node_size)
% regulators
highlight(h,inner_regulator_names,'NodeColor',obj.data.regulation_color,'MarkerSize',obj.data.inner_regulator_size)
highlight(h,outer_regulator_names,'NodeColor',obj.data.regulator_color,'MarkerSize',obj.data.outer_regulator_size)
% proteases
highlight(h,protease_names,'NodeColor',obj.data.protease_color,'MarkerSize',obj.data.protease_size)

set(h,'NodeFontSize',obj.data.node_fontsize,'NodeFontWeight',obj.data.node_fontweight)

if ~isempty(obj.data.G.Edges)
    % find and chategorize the edges
    [sOut,tOut] = findedge(obj.data.G);
    % regulator edges
    activator_edges = obj.data.G.Edges.Weight==weight.regulator;
    repressor_edges = obj.data.G.Edges.Weight==-weight.regulator;
    other_edges = obj.data.G.Edges.Weight==-2*weight.regulator;
    % protease edges
    protease_edges = obj.data.G.Edges.Weight==weight.protease;

    % change the settings for the edges
    highlight(h,sOut(activator_edges),tOut(activator_edges),'EdgeColor',obj.data.act_edge_color,'LineWidth',obj.data.linewidth_regulator,'LineStyle','-')
    highlight(h,sOut(repressor_edges),tOut(repressor_edges),'EdgeColor',obj.data.repr_edge_color,'LineWidth',obj.data.linewidth_regulator,'LineStyle','-')
    highlight(h,sOut(other_edges),tOut(other_edges),'EdgeColor',obj.data.other_edge_color,'LineWidth',obj.data.linewidth_regulator,'LineStyle','-')
    highlight(h,sOut(protease_edges),tOut(protease_edges),'EdgeColor',obj.data.protease_edge_color,'LineWidth',obj.data.linewidth_protease,'LineStyle',':')

    % use our improved repression signs
    if obj.data.original_graph == false


        % we do not give arrow head or |- symbol to inner regulator nodes
        % only if it comes from another inner regulator
        % nohead_edges = false(length(activator_edges),1);
        for i = 1:length(inner_regulator_names)
            % edge coming into the inner regulator node
            incoming_edge = inedges(obj.data.G, inner_regulator_names(i));
            % the node where the edge comes from
            source_node = obj.data.G.Edges.EndNodes(incoming_edge, 1);
            % if the edge does not come from an other inner node, we will
            % not draw a head to the arrow in the graph
            for j = 1:length(source_node)
                if ~any(strcmp(inner_regulator_names,source_node{j}))
                    % nohead_edges(incoming_edge(j)) = 1;
                    activator_edges(incoming_edge(j)) = 0;
                    repressor_edges(incoming_edge(j)) = 0;
                end
            end
        end
        % every edge which is not activator or repressor will be no head:
        nohead_edges = ~(activator_edges|repressor_edges);


        % if we want to redraw the arrows with |- symbols
        % get rid of the original arrows
        set(h,'ArrowSize', 0)
        % activator edges with arrow head
        customizeSelectedArrowheads(h, obj.data.G, findedge(obj.data.G,sOut(activator_edges),tOut(activator_edges)), 'arrow')
        % activator edges with |- head
        customizeSelectedArrowheads(h, obj.data.G, findedge(obj.data.G,sOut(repressor_edges),tOut(repressor_edges)), 'perpendicular')
        % edges to inner regulators without head
        customizeSelectedArrowheads(h, obj.data.G, findedge(obj.data.G,sOut(nohead_edges),tOut(nohead_edges)), 'nohead')
        % redraw protease edges
        customizeSelectedArrowheads(h, obj.data.G, findedge(obj.data.G,sOut(protease_edges),tOut(protease_edges)), 'nohead')
        % bring graph to the top
        uistack(h,'top')
        % get rid of the original graph edges
        h.EdgeColor = 'none';
    end

end

% % regulator routes signed with '<-' and '|-'
% regrouts_signed = obj.calc_regreouts_signed();

% % change the label
% for i = 1:numel(regrouts_signed)
%     regrout = strrep(regrouts_signed{i},'<-','_');
%     regrout = strrep(regrout,'|-','_');
%     labelnode(h,[regrout{:}],[regrouts_signed{i}{:}])
% end



end