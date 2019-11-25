%%
clear all
x_lim_left = -0.015;
x_lim_right = 0.825;  % эперически подобранные константы для выравнивания масштаба осей
%%
[flux,RBDRY,ZBDRY,NBDRY,R,Z,time,rdim,zdim] = gfile_extractor_1t(34363, 000162, 65);
%%
[arr, ind_arr] = min(flux);
[flux_min, min_j] = min (arr);
min_i = ind_arr(min_j);

figure()
grid on
hold on
plot([RBDRY, RBDRY(1)], [ZBDRY, ZBDRY(1)], "o");
magnet_axis = [R(min_j), Z(min_i)];
plot(magnet_axis(1), magnet_axis(2), "*");
xlim([x_lim_left x_lim_right]) 
title("separatrix")
legend("separatrix", "magnetic axis")

%%
%вычисляем отдельно для первой и последней точки, учитывая замкнутость кривой
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
plot([RBDRY, RBDRY(1)], [ZBDRY, ZBDRY(1)], "or");
xlim([x_lim_left x_lim_right]) 
plot(magnet_axis(1), magnet_axis(2), "*r");
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
plot(magnet_axis(1), magnet_axis(2), "*r");
legend("separatrix", "magnetic axis", "grid")

%%