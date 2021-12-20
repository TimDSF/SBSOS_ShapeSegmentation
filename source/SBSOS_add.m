% multiplis two SBSOS statements and return the product
% param:
%   statement_1: augend
%   statement_2: addend
% return:
%   sum: sum of augend and addend

function [sum] = SBSOS_add(statement_1, statement_2)
    sum = SBSOS_simplify([statement_1; statement_2]); % return the simplifies statement
end