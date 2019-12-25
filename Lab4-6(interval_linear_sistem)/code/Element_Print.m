function [] = Element_Print(elem)
    N = length(elem.points);
    for i = 1:N
        display(elem.points(i));  
    end
end