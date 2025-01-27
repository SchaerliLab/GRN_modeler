function [route_out,route_out_simplified] = addinput2route(obj,route)
% if the 'input' is not specified in the route, put the first
% possible 'input' in it (except for the last one)
% in the simplified we put the input only if we have multiple inputs,
% otherwise its an empty character {''}.

% if the 'input' names are not specified in the route
% lets put the first input names in it
first_object = eval(strrep(obj.route2fieldnames(route(1)),'cell','obj'));
if numel(route)==1 || ~any(contains(first_object.Mobj.UserData.input,route{2}))
    route_out = cell(2*numel(route)-1,1);
    route_out(1:2:end) = route;
    route_out_simplified = route;
    % add the inputs step by step
    for i = 1:numel(route)-1
        act_object = eval(strrep(obj.route2fieldnames(route_out(1:2*i-1)),'cell','obj'));
        route_out{2*i} = act_object.Mobj.UserData.input{1};

    end
else
    route_out = route;
    route_out_simplified = route;
end
% now every second element in route is the input except for the last
% one, for the actual new regulator

end

