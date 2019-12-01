function [result] = Element_Create()
POINT = Point_get_interface();
p = POINT.create(0, 0);
result = struct('points', [p], 'index', -1);    
end