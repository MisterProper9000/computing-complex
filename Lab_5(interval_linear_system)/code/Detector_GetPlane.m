function [d] = Detector_GetPlane(det, col)
POINT = Point_get_interface();
DETECTOR = Detector_get_interface();
p = POINT.create(0, 0);
A = DETECTOR.get_col_pos(det, col);

B = det.aperture_pos;
[k, b] = get_line(A, B);
d = abs(k * p.x + (-1) * p.y + b) / sqrt(k^2 + (-1)^2) ;
end