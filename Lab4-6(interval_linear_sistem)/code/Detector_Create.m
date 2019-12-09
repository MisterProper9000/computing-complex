function [result] = Detector_Create()
POINT = Point_get_interface();
p = POINT.create(0, 0);
result = struct('start', p, 'end', p,'step', [], 'direction', p, 'aperture_offset', 0, 'center', p, 'aperture_pos', p, 'horizontal_step', p, 'z_start', 0, 'z_step', 0);   
end