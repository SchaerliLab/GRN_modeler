function str_out = strrep_smart(str_in,str_old,str_new)
% replace str_old to str_new in str_in if the previous or next character is
% not a letter or a '_' or the previous is not a '.' character or the next
% is not a number
% if str_old is surrounded with [] brackets, they will be replaced as well

    str_out = str_in;
    
    if ~contains(str_in,str_old)
        return;
    end

    if contains(str_in,['[' str_old ']'])
        str_old = ['[' str_old ']'];
    end

    pos = strfind(str_in,str_old);

    % if there is a letter (or '_' or '.')  before or after str_new, 
    % or we have a number after that, then str_old will be not
    % replaced
    replace = false(length(pos),1);
    for i = 1:length(pos)
        if pos(i)~=1 && (isletter(str_in(pos(i)-1)) || str_in(pos(i)-1)=='_' || str_in(pos(i)-1)=='.')
            continue;
        end
        if pos(i)+length(str_old)-1<length(str_in) && (isletter(str_in(pos(i)+length(str_old))) || str_in(pos(i)+length(str_old))=='_' || (isstrprop(str_in(pos(i)+length(str_old)),'digit')))
            continue;
        end
        replace(i) = true;
    end

    if ~any(replace)
        return;
    end

    if all(replace)
        str_out = strrep(str_in,str_old,str_new);
        return;
    end

    pos(end+1) = length(str_in)+1;
    str_out = str_in(1:pos(1)-1);
    for i = 2:length(pos)
        if replace(i-1) == true
            str_out = [str_out,strrep(str_in(pos(i-1):pos(i)-1),str_old,str_new)];
        else
            str_out = [str_out,str_in(pos(i-1):pos(i)-1)];
        end
    end
end

