function [result] = split_line_linspace(A, B, n)
    result = {};
    step = (B - A) / n;
    for i=1:n-1
        result = {result{:} A + step * (i)};        
    end
end