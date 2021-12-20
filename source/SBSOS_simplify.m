% simplifies a SBSOS statement
% param:
%   input: input statement
% return:
%   res: simplified statement

function [res] = SBSOS_simplify(input)
    [num_terms, num_variables] = size(input);
    res = sortrows(input, -(1:num_variables));
    
    idx = 1;
    for i = 2:num_terms
        if res(i,1:num_variables-1) == res(idx, 1:num_variables-1)
            res(idx,num_variables) = res(i, num_variables) + res(idx, num_variables);
        else
            idx = idx + 1;
            res(idx,:) = res(i,:);
        end
    end
    
    res = res(1:idx, :);
    res = res(res(:, num_variables) ~= 0, :);
end