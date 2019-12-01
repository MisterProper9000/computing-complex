function [R_segments_arr, Z_segments_arr, lines_start, lines_end] = get_web_grid(cur_cut_y, cur_cut_z, cur_magnet_axis, N, L)
    cur_rad = get_curv_radius_arr( cur_cut_y,  cur_cut_z, length(cur_cut_y));
    cur_rad_smoothed = medfilt1(cur_rad);

    
    init_sector = 1:length(cur_cut_y);
    [s1, s2, s3, s4] = split_sector(init_sector, cur_rad_smoothed);
    sectors = {s1, s2, s3, s4};
    result = {};
    for i = 1:length(sectors)
        tmp_sector = sectors{i};
        tmp = split_linspace(tmp_sector, N);
        result = {result{:} tmp{:}};
    end
    sectors = result;

    segments = [];
    middle = cur_magnet_axis;
    R_segments_arr = {};
    Z_segments_arr = {}; 
    R_segments = [];
    Z_segments = [];


    for i = 1:length(sectors)
        tmp_sector = sectors{i};
        cur_ind = tmp_sector(1);
        cur_point = [cur_cut_y(cur_ind),cur_cut_z(cur_ind)];
        split_points = split_line_linspace(cur_point, middle, L);
        for k = 1:length(split_points)
            cur_split_point = split_points{k};
            %mid_point = get_center(cur_point, middle);
            R_segments_arr{k}(i) = cur_split_point(1);
            Z_segments_arr{k}(i) = cur_split_point(2);
        end
    end
    
    lines_start = {};
    lines_end = {};
    for i = 1:length(sectors)
        tmp_sector = sectors{i};
        cur_ind = tmp_sector(1);
        cur_point = [cur_cut_y(cur_ind),cur_cut_z(cur_ind)];
        lines_start{i} = cur_point;
        lines_end{i} = middle;
    end
    
end