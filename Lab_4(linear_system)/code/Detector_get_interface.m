function [result] = Detector_get_interface()
    result = struct('create', @Detector_Create, 'get_ray', @Detector_GetRay, 'get_plane', @Detector_GetPlane, 'get_bottom', @Detector_Get_Bottom, 'get_col_pos', @Detector_Get_Col_Pos);   
end