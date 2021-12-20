% multiplis two SBSOS statements and return the product
% param:
%   constant: multiplicand
%   statement: multiplier
% return:
%   product: product

function [product] = SBSOS_constant_multiply(constant, statement)
    product = statement;
    product(:, size(product, 2)) = product(:, size(product, 2)) * constant;
end