function [bottom] = Detector_Get_Bottom(det)
bottom = 0 - det.step.y * 8;
end