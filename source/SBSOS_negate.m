% multiplis two SBSOS statements and return the product
% param:
%   statement: input
% return:
%   negation: negation of the input

function [negation] = SBSOS_negate(statement)
    num_variables = size(statement, 2); % number of variables
    
    negation = statement;
    negation(:, num_variables) = - negation(:, num_variables); % negate the coefficents of each monomial
end