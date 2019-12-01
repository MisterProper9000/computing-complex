[flux,RBDRY,ZBDRY,NBDRY,R,Z,time,rdim,zdim] = gfile_extractor_1t(34363, 000162, 65);

[arr, ind_arr] = min(flux);
[flux_min, min_j] = min (arr);
min_i = ind_arr(min_j);
 
figure()
grid on
hold on
plot([RBDRY, RBDRY(1)], [ZBDRY, ZBDRY(1)]);
magnet_axis = [R(min_j), Z(min_i)];
plot(magnet_axis(1), magnet_axis(2), "o");
title("separatrix")
legend("separatrix", "magnetic axis")
 

a = [RBDRY(NBDRY),ZBDRY(NBDRY)];
b = [RBDRY(1),ZBDRY(1)];
c = [RBDRY(2),ZBDRY(2)];      
r = get_curv_radius(a, b, c);
for i=2:NBDRY - 1
    a = [RBDRY(i-1),ZBDRY(i-1)];
    b = [RBDRY(i),ZBDRY(i)];
    c = [RBDRY(i + 1),ZBDRY(i + 1)];
    tmp_r = get_curv_radius(a, b, c);
    r = [r, tmp_r];
end
 
a = [RBDRY(NBDRY-1),ZBDRY(NBDRY-1)];
b = [RBDRY(NBDRY),ZBDRY(NBDRY)];
c = [RBDRY(1),ZBDRY(1)];
tmp_r = get_curv_radius(a, b, c);
r = [r, tmp_r];
figure()
grid on
hold on
plot(r);
title("radius of curvature")
 
%%
 
r_smoothed = medfilt1(r);
figure()
grid on
hold on
plot(r_smoothed);
title("smoothered radius of curvature")
%%
% sector = 1:NBDRY;
% sectors = {sector};
% N = 2;
% while length(sectors) < 4 * N
%     result = {};
%     for i = 1:length(sectors)
%         tmp_sector = sectors{i};
%         if( length(tmp_sector) <= 4)
%             result = {result{:} tmp_sector};
%         else
%     
%         [s1, s2, s3, s4] = split_sector(tmp_sector, r);
%         result = {result{:} s1 s2 s3 s4};
%         end
%     end
%     sectors = result
% end
 
N = 5;
init_sector = 1:NBDRY;
[s1, s2, s3, s4] = split_sector_easy(init_sector, r_smoothed);
sectors = {s1, s2, s3, s4};
result = {};
for i = 1:length(sectors)
    tmp_sector = sectors{i};
    tmp = split_linspace(tmp_sector, N);
    result = {result{:} tmp{:}};
end
sectors = result;
 
segments = [];
middle = magnet_axis;
R_segments = [];
Z_segments = [];
 
% src_segments = []
% for i = 1:length(sectors)
%     tmp_sector = sectors{i};
%     cur_ind = tmp_sector(1);
%     tmp_segment = [R_segments, mid_point(1)]
%     [Z_segments, mid_point(2)]
%     
%         cur_point = [RBDRY(cur_ind),ZBDRY(cur_ind)];
%     plot([cur_point(1), middle(1)], [cur_point(2), middle(2)], 'b')
%     
% end
for i = 1:length(sectors)
    tmp_sector = sectors{i};
    cur_ind = tmp_sector(1);
    cur_point = [RBDRY(cur_ind),ZBDRY(cur_ind)];
    mid_point = get_center(cur_point, middle);
    R_segments = [R_segments, mid_point(1)];
    Z_segments = [Z_segments, mid_point(2)];
end
 
 
figure()
grid on
hold on
plot([RBDRY, RBDRY(1)], [ZBDRY, ZBDRY(1)], "r");
plot(magnet_axis(1), magnet_axis(2), "or");
title("grid")
 
 
for i = 1:length(sectors)
    tmp_sector = sectors{i};
    cur_ind = tmp_sector(1);
    cur_point = [RBDRY(cur_ind),ZBDRY(cur_ind)];
    plot([cur_point(1), middle(1)], [cur_point(2), middle(2)], 'b')
end
 
for i = 1:length(R_segments)-1
    plot([R_segments(i), R_segments(i + 1)], [Z_segments(i), Z_segments(i + 1)], 'b')
end
plot([R_segments(length(R_segments)), R_segments(1)], [Z_segments(length(R_segments)), Z_segments(1)], 'b')
plot(magnet_axis(1), magnet_axis(2), "or");
legend("separatrix", "magnetic axis", "grid")
 


function r = get_curv_radius(A, B, C)
    a = B - A;
    b = C - B;
    c = C - A;
    cos = sum(c.*b)/(norm(c)*norm(b));
    sin = sqrt(1 - cos^2);
    r = norm(a) / (2 *sin);
end


function [s1, s2, s3, s4] = split_sector(sector, r)
    separator_ind = fix(length(sector) / 2);  % делим нацело на 2 части
 
    first_sector = sector(1:separator_ind);
    second_sector = sector(separator_ind + 1:length(sector));
    [min_val , first_min_ind] = min(r(first_sector));
    [max_val , first_max_ind] = max(r(first_sector));
    
    [min_val , second_min_ind] = min(r(second_sector));
    [max_val , second_max_ind] = max(r(second_sector));
    second_min_ind = second_min_ind + length(first_sector); % т.к. индекс получаем относительно массива r(second_sector) (36 элементов), а хотим иметь абсолютный
    second_max_ind = second_max_ind + length(first_sector);
 
    sep_1 = min(first_min_ind, first_max_ind);
    sep_2 = max(first_min_ind, first_max_ind);
    sep_3 = min(second_min_ind, second_max_ind);
    sep_4 = max(second_min_ind, second_max_ind);
    if(sep_4 == length(sector))
        s1 = sector(1 : sep_1);
    else
        s1 = sector([sep_4 + 1: length(sector), 1 : sep_1]);
    end        
    s2 = sector(sep_1 + 1 : sep_2);
    s3 = sector(sep_2 + 1 : sep_3);
    s4 = sector(sep_3 + 1: sep_4);
end


function [s1, s2, s3, s4] = split_sector_easy(sector, r)
    separator_ind = fix(length(sector) / 2);  % делим нацело на 2 части
 
    first_sector = sector(1:separator_ind);
    second_sector = sector(separator_ind + 1:length(sector));
    [max_val , first_max_ind] = max(r(first_sector));
    [max_val , second_max_ind] = max(r(second_sector));
   
    second_max_ind = second_max_ind + length(first_sector);% т.к. индекс получаем относительно массива r(second_sector) (36 элементов), а хотим иметь абсолютный
 
    sep_1 = first_max_ind;
    sep_2 = length(first_sector);
    sep_3 = second_max_ind;
    sep_4 = length(sector);
    
    s1 = sector(1 : sep_1);    
    s2 = sector(sep_1 + 1 : sep_2);
    s3 = sector(sep_2 + 1 : sep_3);
    s4 = sector(sep_3 + 1: sep_4);
end


function [result] = split_linspace(sector, N)
    result = {};
    cur_step = ceil(length(sector) / N);
    cur_left = 1;
    cur_right = cur_step;
    reduse = ceil(length(sector));
    for i = 1:N
        if(reduse == 0)
            return
        end
        tmp = sector(cur_left:cur_right);
        result = {result{:} tmp};
        reduse = reduse - length(tmp);
        cur_left = cur_right + 1;
        cur_step = ceil(reduse / (N - i));
        cur_right = cur_right + cur_step;  
    end 
end


function [result] = split_linspace(A, B)
    result = A + (B - A) / 2;  
end
