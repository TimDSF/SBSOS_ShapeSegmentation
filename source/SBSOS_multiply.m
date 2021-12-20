% multiplis two SBSOS statements and return the product
% param:
%   statement_1: multiplicand
%   statement_2: multiplier
% return:
%   product: product

function [product] = SBSOS_multiply(statement_1, statement_2)
    num_variables = size(statement_1, 2); % number of variables, with last one being constant terms
    num_terms_1 = size(statement_1, 1); % number of terms in state_1
    num_terms_2 = size(statement_2, 1); % number of terms in state_2
    
    product = zeros(num_terms_1 * num_terms_2, num_variables); % initialize the product
    
    idx = 0; % index of the resulting terms in the product
    for i = 1:num_terms_1
        for j = 1:num_terms_2
            idx = idx + 1;
            
            product(idx, :) = statement_1(i, :) + statement_2(j, :); % sum of the powers
            product(idx, num_variables) = statement_1(i, num_variables) * statement_2(j, num_variables); % product for the coefficients
        end
    end
    
    product = SBSOS_simplify(product); % return the simplifies product
end