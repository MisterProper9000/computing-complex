function [result] = Point_get_interface()
    result = struct('create', @Point_Create, 'dist', @Point_Dist);
end
