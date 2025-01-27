function success = checkstruct(obj,str,start_pos)
% checks wheter the fields in a nested structure or object exist or not
% it gives a warning at the missing fields
% str: a string to the field in question ('l.g.k')
% startpos: a position where we start checking (without it we start at the
% beginning)

% str can have a format like l.('g').('k') => l.g.k
str = strrep(str,'(''','');
str = strrep(str,''')','');

% success = true/false
if nargin == 2
    start_pos = 1;
end
% rename struct
eval([str(1:find(str=='.',1,'first')-1) ' = obj;'])
% position of the next dots in str
% first field(s) in str
pos1 = find(str(start_pos:end)=='.',1,'first');
pos1 = start_pos+pos1-1;
% next field in str
pos2 = find(str(pos1+1:end)=='.',1,'first');
if isempty(pos2)
    % if that was the last field to check
    pos2 = numel(str)+1;
else
    pos2 = pos1+pos2;
end
% check wheter we have the next field or not
if ~isfield(eval(str(1:pos1-1)),str(pos1+1:pos2-1)) && all(~isprop(eval(str(1:pos1-1)),str(pos1+1:pos2-1)))
    warning('We did not find the field ''struct%s''',str(find(str=='.',1,'first'):pos2-1));
    success = false;
    return
elseif pos2 == numel(str)+1
    % if that was the last one, we are happy
    success = true;
    return
else
    % call the next chunk recursively
    start_pos = pos2;
    success = checkstruct(obj,str,start_pos);
    return;
end