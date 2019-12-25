%%
clear all
%close all
POINT = Point_get_interface();
ELEMENT = Element_get_interface();
DETECTOR = Detector_get_interface();

%--------------      lab 4      --------------
is_need_plot(1) = false;    %график сепаратрисы и магнитной оси
is_need_plot(2) = false;    %график сетки разбиения
is_need_plot(3) = false;    %график сетки разбиения сепаратриссы в сечении плоскость x = H
is_need_plot(4) = false;    %график плоскостей сечения токамака
is_need_plot(5) = false;    %график лучей для 16-го сечения
is_need_plot(6) = false;    %графики лучей и пересечений для всех плоскостей сечений
is_need_plot(7) = false;    %график плоскости детектора
is_need_plot(8) = true;    %построение матрицы длинн хорд по всем плоскостям
is_need_plot(9) = false;    %легенда для графиков лучей и пересечений для всех плоскостей сечений
%--------------      lab 5      --------------
% для 5 лабы нужен is_need_plot(8)= true
is_need_plot(10) = false;    %обзор интегральной светимости
is_need_plot(11) = false;    % x = b / A
is_need_plot(12) = false;    % x = (A' * A)^(-1) * A' * b
is_need_plot(13) = false;    % tolsolvty
is_need_plot(14) = false;    % графики Ax, sup_b, inf_b для каждой итерации
is_need_plot(15) = false;    % графики полученного решения x_i от i
is_need_plot(16) = false;    % оценка числа обусловленности матрицы А и вычисление IVE
is_need_plot(17) = false;    % вычисление IVE интервальной матрицы А
%--------------      lab 6      --------------
% для 6 лабы нужен is_need_plot(8)= true
is_need_plot(18) = false;    % гистограмма x решения ЗЛП
is_need_plot(20) = false;    % график x решения ЗЛП
is_need_plot(19) = false;    % график w решения ЗЛП
is_need_plot(21) = false;    % график значений Ax решения ЗЛП

input_file_index = 37000;
input_file_name = strcat("data\", num2str(37000), "_SPD16x16.mat");
input_time_period = 000162;
%%
[flux, RBDRY, ZBDRY, NBDRY, R, Z, time, rdim, zdim] = gfile_extractor_1t(input_file_index, input_time_period, 65);

%изменим начало обхода и измени обход сепаратрисы на противчасовой 
%т.к. последняя точка равна первой, учтём это
%новое начало - 32-ой элемент
[RBDRY,ZBDRY] = circle_spin_and_reverse(RBDRY,ZBDRY,NBDRY, 32);

% вычисления магнитной оси как минимум магнитного потока flux
[arr, ind_arr] = min(flux);
[flux_min, min_j] = min (arr);
min_i = ind_arr(min_j);
magnet_axis = [R(min_j), Z(min_i)];

%построение графика сепаратриссы и магнитной оси
if(is_need_plot(1))
    figure()
    grid on
    hold on
    plot(RBDRY, ZBDRY, "o");
    plot(magnet_axis(1), magnet_axis(2), "*");
    [x_limit, y_limit] = get_plot_lim(RBDRY, ZBDRY);
    xlim(x_limit) 
    ylim(y_limit) 
    axis equal
    title("separatrix")
    legend("separatrix", "magnetic axis")
% вывод индексы точек на графике
%     for i = 1:NBDRY
%         text(RBDRY(i), ZBDRY(i), num2str(i))
%     end
end
%%
%%
%num_sectors - число радиальных линий сетки внутри каждой четверти
%num_circle - число кольцевых линий сетки
num_sectors = 8;
num_circle = 6;
cur_cut_y = RBDRY;
cur_cut_z = ZBDRY;
cur_magnet_axis = magnet_axis;

%get_web_grid - разбиение сепаратрисы на сетку
%CreateElementsFromGrid - построение секторов разбиения по сетке
[R_segments_arr, Z_segments_arr, lines_start, lines_end] = get_web_grid(cur_cut_y, cur_cut_z, cur_magnet_axis, num_sectors, num_circle);
elements = CreateElementsFromGrid(R_segments_arr, Z_segments_arr, lines_start, lines_end);

%построение графика сетки разбиения сепаратриссы
if(is_need_plot(2))
    figure()
    grid on
    hold on
    axis equal
    N = length(elements);
    for i = 1:N
        ELEMENT.draw(elements(i), "b", i);
    end
    legend("grid");
end

%H - расстояние от цетра токомака до плоскости сечения
min_x_sep = min(RBDRY);
max_x_sep = max(RBDRY);
%H=(min_x_sep + max_x_sep)/2;
%H=min_x_sep - (max_x_sep - min_x_sep)/16;
%H=min_x_sep;
%H=min_x_sep;
%H=0.4;
H = 0;
%H = 0.36351
%Element.get_cut - возвращает результат пересечения фигуры вращения сектора
%и плоскости x = H (возвращает секторы и их количество)
h = length(elements);
cut_elements = [];
left_cut_elem = [];
right_cut_elem = [];
for i = 1:h
    [elem, count] = ELEMENT.get_cut(elements(i), H);
    if(count == 2)
        left_cut_elem = [left_cut_elem, elem(1)];
        right_cut_elem = [right_cut_elem, elem(2)];
    else
        left_cut_elem = [left_cut_elem, elem]; 
    end
    %cut_elements = [cut_elements, elem];
end

cut_elements = [left_cut_elem, right_cut_elem];

%построение графика сетки разбиения сепаратриссы в сечении плоскость x = H
if(is_need_plot(3))
    figure()
    grid on
    hold on
    axis equal
    N = length(cut_elements);
    %Element.draw - строит сектор на текущей figure, первый аргумент - сектор,
    %второй - цвет
    for i = 1:N
        ELEMENT.draw(cut_elements(i), "b", i);
    end
    title(strcat("grid Cut H = ", num2str(H)));
    legend("grid");
end
%%
%% магические константа от преподавателя
% угол между направлением камеры-обскуры и направлением на центр (между 8 и 9 лучами)
ang = acos((708^2 + 720^2 - 31^2) / (2 * 708 * 720));
% положение края детектора (1-го столбца)
spd_start = [0, -0.708];
% положение 16-го столбца
spd_end = [0.72 * sin(ang), 0.72 * -cos(ang)];
% вектор направления камеры-обскуры в экваториальной плоскости
spd_vect = (spd_end - spd_start) / norm(spd_end - spd_start);
% шаг между столбцами в плоскости детектора, 2 числа
spd_xy_step = [2.3375 - 0.88 , 3.81 - 2.3375 + 0.88 ] * 1e-03;
% центр детектора
pp = spd_start + spd_vect * ((spd_xy_step(1) + spd_xy_step(2)) * 8 + 0.52 * 1e-03) / 2;% + spd_vect * 0.35 / 2 * 1e-03;
% отступ апертуры от центра детектора
aperture_xy_offset = 0.0395;
% координата апертуры
aperture_xy = [pp(1) - spd_vect(2) * aperture_xy_offset, pp(2) + spd_vect(1) * aperture_xy_offset];
%spd xz – устройство детектора в меридиональной плоскости
spd_z_start = (27.52 - 0.49) / 2 * 1e-03;
spd_z_step = -1.72 * 1e-03;
spd_xy = spd_start + spd_vect * (spd_xy_step(2) / 2 + 0.26 * 1e-03);
%%

detector = DETECTOR.create();
detector.start = POINT.create(spd_start(1), spd_start(2));
detector.end = POINT.create(spd_end(1), spd_end(2));
detector.step = spd_xy_step;
detector.direction = POINT.create(spd_vect(1), spd_vect(2));
detector.center = POINT.create(pp(1), pp(2));
detector.aperture_offset = 0.0395;
detector.aperture_pos = POINT.create(aperture_xy(1), aperture_xy(2)); 
detector.z_start = spd_z_start;
detector.z_step = spd_z_step;
%тут установил равномерный шаг, чтобы в плоскости детектора поместили 16 столбцов
detector.horizontal_step = POINT.create((detector.end.x - detector.start.x) / 16 , (detector.end.y - detector.start.y) / 16 ); 

%построение графика сечений токомака
if(is_need_plot(4))
    
    figure()
    grid on
    hold on
    axis equal
    phi = -pi:pi/360:pi;
    R = max(RBDRY);
    r = min(RBDRY);
    plot(R * cos(phi), R * sin(phi));
    plot(r * cos(phi), r * sin(phi));
    plot(detector.aperture_pos.x, detector.aperture_pos.y, "ob")
    
    
    x = detector.start.x;
    y = detector.start.y;
    for i=1:16
        Pos = DETECTOR.get_col_pos(detector, i);
        plot(Pos.x, Pos.y, "or")  
    end
    plot(detector.center.x, detector.center.y, "*b");
    plot(detector.start.x, detector.start.y, "*b") 
    plot(detector.end.x, detector.end.y, "*b")
    
    for i=1:16
        
        A = DETECTOR.get_col_pos(detector, i);
        x = A.x:0.01:0.8;
        B = detector.aperture_pos;
        [k, b] = get_line(A, B);
        plot(x, k*x+b, "r")
        %просто красивый вывод надписи
        text_R = R*1.2;
        text_x = (-2*b*k + sqrt(4*b^2*k^2 - 4*(k^2+1)*(b^2 - text_R^2)))/(2*(k^2+1));
        text_y = k*text_x + b;
        text(text_x, text_y, num2str(i));
    end

    xlim([-0.8, 0.8])
    ylim([-0.8, 0.8])
    
    %xlim([-0.2, 0.2])
    %ylim([-0.8, -0.5])
end
%%
%посроение графика лучей для 16-го сечения
if(is_need_plot(5))
    cut_ind = 16;
    figure()
    grid on
    hold on
    axis equal
    H = DETECTOR.get_plane(detector, cut_ind);
    title(strcat("grid Cut H = ", num2str(H)))

    h = length(elements);
    cut_elements = [];
    draw_cut_elements = [];
    for i = 1:h
        [elem, count] = ELEMENT.get_cut(elements(i), H);
        draw_cut_elements = [draw_cut_elements, elem];
        if(H < min_x_sep)
        %область двусвязная, не вычисляем хорды для правого графика
            if(count == 2) 
                cut_elements = [cut_elements, elem(2)];
            else
                cut_elements = [cut_elements, elem];
            end
        else
        %область односвязная
            cut_elements = [cut_elements, elem];
        end
    end
    N = length(draw_cut_elements);
    for i = 1:N
        ELEMENT.draw(draw_cut_elements(i), "b");
    end
    N = length(cut_elements); 
    for ray_ind=1:16
        [k, b, det_pos, apper_pos] = DETECTOR.get_ray(detector, cut_ind, ray_ind, true);
        x = det_pos.x:0.01:0.7;
        plot(x, k*x+b, "r")
        plot(det_pos.x, det_pos.y, "*r")

        for t = 1:N
            ELEMENT.draw(cut_elements(t), "b")
            [hord, intersection] = ELEMENT.get_hord(cut_elements(t), k, b);          
            for g = 1:length(intersection)
                plot(intersection(g).x, intersection(g).y, "om");
            end
        end
        %вывод номеров лучей
        text_x = 0.6;
        text_y = text_x * k + b;
        text(text_x, text_y, num2str(ray_ind));
    end
    plot(apper_pos.x, apper_pos.y, "ob")
    %xlim([-0.8, -0.6])
    %ylim([-0.2, 0.2])
end
%%

%графики лучей и пересечений для всех плоскостей сечений
if(is_need_plot(6))
    for cut_ind=1:16
        figure()
        grid on
        hold on
        axis equal
        H = DETECTOR.get_plane(detector, cut_ind);
        title(strcat("grid Cut H = ", num2str(H)));
        element_num = length(elements);
        cut_elements = [];
        for i = 1:element_num
            [elem, count] = ELEMENT.get_cut(elements(i), H);
            if(H < min_x_sep)
            %область двусвязная, не вычисляем хорды для правого графика
                if(count == 2)
                    cut_elements = [cut_elements, elem(2)];
                else
                    cut_elements = [cut_elements, elem];
                end
            else
                %область односвязная
                cut_elements = [cut_elements, elem];
            end
        end
        N = length(cut_elements);
        
        for i = 1:N
            %ELEMENT.draw(cut_elements(i), "b", i);
        end
        
        for ray_ind=1:16
            [k, b, det_pos, apper_pos] = DETECTOR.get_ray(detector, cut_ind, ray_ind);
            x = det_pos.x:0.01:0.6;
            plot(x, k*x+b, "r")
            plot(det_pos.x, det_pos.y, "*r")
            for t = 1:N
                ELEMENT.draw(cut_elements(t), "b")
                [hord, intersection] = ELEMENT.get_hord(cut_elements(t), k, b);          
                for g = 1:length(intersection)
                    plot(intersection(g).x, intersection(g).y, "og");
                end
            end
        end    
        plot(apper_pos.x, apper_pos.y, "ob")
    end
end
%%
%график плоскости детектора
if(is_need_plot(7))
    figure()
	grid on
    hold on
	%axis equal
    spd_sizes = [0.88, 1.23] * 1e-3;

    y_square_up = [0, spd_sizes(1), spd_sizes(1), 0, 0];
    y_square_down = [0, 0, spd_sizes(1), spd_sizes(1), 0];
    z_square_up = [0, 0, spd_sizes(2), spd_sizes(2), 0];
    z_square_down = [0, -spd_sizes(2), -spd_sizes(2), 0, 0];
    color = "k";
    point_color = ".k";
    for j = 0:7
        for i = 0:7
            x = (spd_sizes(1) + j * sum(spd_xy_step) + y_square_up);
            y = (-spd_z_step - spd_sizes(2)) / 2 + i * (-spd_z_step) + z_square_up;
            plot(x, y, color);
            x = x(1:4);
            x = sum(x)/4;
            y = y(1:4);
            y = sum(y)/4;
            plot(x, y, point_color);
            x = spd_sizes(1) + spd_xy_step(1) + j * sum(spd_xy_step) + y_square_up;
            y = (-spd_z_step - spd_sizes(2)) / 2 + i * (-spd_z_step) + z_square_up;
            plot(x, y, color);
            x = x(1:4);
            x = sum(x)/4;
            y = y(1:4);
            y = sum(y)/4;
            plot(x, y, point_color);
            x = spd_sizes(1) + j * sum(spd_xy_step) + y_square_down;
            y = (spd_z_step + spd_sizes(2)) / 2 + i * spd_z_step + z_square_down;
            plot(x, y, color);
            x = x(1:4);
            x = sum(x)/4;
            y = y(1:4);
            y = sum(y)/4;
            plot(x, y, point_color);
            x = spd_sizes(1) + spd_xy_step(1) + j * sum(spd_xy_step) + y_square_down;
            y = (spd_z_step + spd_sizes(2)) / 2 + i * spd_z_step + z_square_down;
            plot(x, y, color);
            x = x(1:4);
            x = sum(x)/4;
            y = y(1:4);
            y = sum(y)/4;
            plot(x, y, point_color);
        end
    end
end
%%
%%

%построение матрицы длинн хорд по всем плоскостям
element_num = length(elements);
hord_matrix = zeros(256, element_num);
cur_row = 1;
if(is_need_plot(8))
    for cut_ind=1:16
        H = DETECTOR.get_plane(detector, cut_ind);
        cut_elements = [];
        for i = 1:element_num
            [elem, count] = ELEMENT.get_cut(elements(i), H);
            if(H < min_x_sep)
            %область двусвязная, не вычисляем хорды для правого графика
                if(count == 2)
                    cut_elements = [cut_elements, elem(2)];
                else
                    cut_elements = [cut_elements, elem];
                end
            else
                %область односвязная
                cut_elements = [cut_elements, elem];
            end
        end
        N = length(cut_elements);
        for ray_ind=1:16
            [k, b, det_pos, apper_pos] = DETECTOR.get_ray(detector, 16, ray_ind);
            for t = 1:N
                [hord, intersection] = ELEMENT.get_hord(cut_elements(t), k, b); 
                hord_matrix(cur_row, cut_elements(t).index) = hord_matrix(cur_row, cut_elements(t).index) + hord;
            end
            %переходим к следующей строке т.к. будем обрабатывать следующий луч
            cur_row = cur_row + 1;
        end
    end
end
%%
%легенда для графиков сечений
if(is_need_plot(9))
    figure()
    grid on
    hold on
    plot([1,2], [1,2], "b");
    plot([1,2], [1,2], "r");
    plot([1,2], [1,2], "*r");
    plot([1,2], [1,2], "ob");
    plot([1,2], [1,2], "og");

    legend("grid", "ray", "detector`s pixel", "aperture", "intersection point")
end
%%
%--------------------------------------------------------------------------
% lab 5
%--------------------------------------------------------------------------
%%
%матрица построена в is_need_plot(8)
%A = hord_matrix;
input_data = load(input_file_name);
sign_bb = input_data.sign_bb(:, :, :);
cnt_meas = size(sign_bb, 3);
tp = cell2mat(input_data.Data(1, 2)) * 1e-3;% - шаг по времени
tz = cell2mat(input_data.Data(2, 2));
t_start = tz;
t_end = t_start + (cnt_meas - 1) * tp;
t_i = t_start:tp:t_end;

%%
%обзор интегральной светимости
if(is_need_plot(10))
    t_cons_start = 125;
    t_cons_end = 200;

    dt_cons = 1;
    start_efit_time_i = t_cons_start:dt_cons:t_cons_end;

    B = [];
    for start_efit_time1 = t_cons_start:t_cons_end
        ind = find(abs(t_i - start_efit_time1) < tp/2);
        b = [];
        for i = 16:-1:1
            b = [b; sign_bb(16:-1:1, i, ind(1))];
        end
        b = double(b);
        Bnew = sum(b(:));
        B=[B, Bnew];
    end
    figure()
    plot([t_cons_start:t_cons_end],B);
    title_str_name=['36917SPD16x16.matt', ' Sum b'];
    title(title_str_name);
    xlabel('start efit time');
    figure_name_out = strcat(title_str_name,'.png');
    %print('-dpng', '-r300', figure_name_out), pwd;
end   

%%
% выбираем "окно" по котору мычисляем границы b
b_time_window = 1;

%извлекаем значения по нужным периода времени 
b_data = [];
for start_efit_time1 = input_time_period - b_time_window:input_time_period + b_time_window
    ind = find(abs(t_i - start_efit_time1) < tp/2);
    b = [];
    for i = 16:-1:1
        b = [b; sign_bb(16:-1:1, i, ind(1))];
    end
    b = double(b);
    b_data = [b_data, b];
end

% вычисляем верхнюю и нижнюю границе столбца свободных членов
N = length(b_data);
inf_b = zeros(N, 1); 
sup_b = zeros(N, 1); 
for i = 1:N
    inf_b(i) = min(b_data(i, :));
    sup_b(i)  = max(b_data(i, :));
end
%%
b = (sup_b + inf_b)/2;
A = hord_matrix;


%x = b / A
if(is_need_plot(11)) 
   
    b_tmp = b
    x1 = b'/A'
    figure()
    hold on;
    grid on  
    plot(x1, "o");
    ylabel("x_i")
    xlabel("i")
    figure()
    hold on;
    grid on;
    hist(x1);
    ylabel("numder")
    xlabel("x_i")
    title("x = b'/A'")
    
    b_tmp = b
    x11 = A\b;
    figure()
    hold on;
    grid on  
    plot(x11, "o");
    ylabel("x_i")
    xlabel("i")
    figure()
    hold on;
    grid on;
    hist(x11);
    ylabel("numder")
    xlabel("x_i")
    title("x = b\A")
            
end

%x = (A'A) * A' * b
if(is_need_plot(12))
    disp(strcat("cond(A) = ", num2str(cond(A))));

    disp(strcat("cond(A'A) = ", num2str(cond(A'*A))));
    x2 = inv((A'*A)) *A' *b; 
    lambda = eig(A'*A);
    figure()
    hist(lambda)  
    title("Собственные значения матрицы А'A")
    ylabel("numder")
    xlabel("lambda")
    title("x = inv((A'*A)) *A' *b")
    
    ind = find(lambda > 0.2);
    disp(strcat("numder of lambda > 0.2: ", num2str(length(ind))));
    
end


%tolsolvty
if(is_need_plot(13))
    b_sup_tmp = sup_b;
    b_inf_tmp = inf_b;
    A_inf = A;
    A_sup = A;
    [tolmax, argmax, envs, ccode] =  tolsolvty(A_inf, A_sup, b_inf_tmp,  b_sup_tmp);
    disp(strcat("tolmax = ", num2str(tolmax)));
    iter = 0;
    disp(strcat("iteration = ", num2str(iter)));
    % графики Ax, sup_b, inf_b для каждой итерации
    if(is_need_plot(14))
        figure()
        hold on;
        grid on;    
        % т.к. сейчас А - не интервальная матрица, то (A_inf == A_sup)
        % и для вычисления значения можем взять любую
        plot_A = A_sup * argmax;
        x_axis = 1:length(b_sup_tmp);
        plot(x_axis, plot_A');
        plot(x_axis, b_inf_tmp')
        plot(x_axis, b_sup_tmp');
        legend("Ax", "inf b", "sup b")
        ylabel("value")
        xlabel("i")
        title("Первое решение tolsolvty");
    end
    
    if(tolmax < 0)
        %проверяем решение и расширяем границы b
        iter = iter + 1;
        disp(strcat("iteration = ", num2str(iter)));
        shift = -tolmax;
        disp(strcat("delta b = ", num2str(shift)));
        b_sup_tmp = b_sup_tmp + shift;
        b_inf_tmp = b_inf_tmp - shift;
        [tolmax, argmax, envs, ccode] =  tolsolvty(A_inf, A_sup, b_inf_tmp,  b_sup_tmp);
        disp(strcat("tolmax = ", num2str(tolmax)));
        % графики Ax, sup_b, inf_b для каждой итерации
        if(is_need_plot(14))
           figure()
            hold on;
            grid on;    
            % т.к. сейчас А - не интервальная матрица, то (A_inf == A_sup)
            % и для вычисления значения можем взять любую
            plot_A = A_sup * argmax;
            x_axis = 1:length(b_sup_tmp);
            plot(x_axis, plot_A');
            plot(x_axis, b_inf_tmp')
            plot(x_axis, b_sup_tmp');
            legend("Ax", "inf b", "sup b")
            ylabel("value")
            xlabel("i")
            title("Второе( с расширенным b) решение tolsolvty");
        end
    end
    %графтк полученного решения x_i от i
    if(is_need_plot(15))
            figure()
            hold on;
            grid on;
            plot(argmax, "o");
            ylabel("x_i")
            xlabel("i")
            title("решение tolsolvty")
            figure()
            title("Гистограмма решения tolsolvty")
            hold on;
            grid on;
           	hist(argmax);
            ylabel("numder")
            xlabel("x_i")
    end
    
    cond_A = my_HeurMinCond(A_inf, A_sup);
    A_IVE = my_IVE(A_inf, A_sup,  b_inf_tmp,  b_sup_tmp, tolmax, argmax, length(argmax));
end

   

%%
%оценка числа обусловленности матрицы А
if(is_need_plot(16))
    %радиус элементов матрицы А - 10% от их величины
    matrix_radius = 0.1;
    A_inf_1 = A * (1 -  matrix_radius);
    A_sup_1 = A * (1 +  matrix_radius);
    b_sup_tmp_1 = sup_b;
    b_inf_tmp_1 = inf_b;
    
    for i = 10:10:100
        
        cond_A = my_HeurMinCond(A_inf_1, A_sup_1, i);  
        disp(strcat("rad = ", num2str(matrix_radius)," : HeurMinCond(A, ", num2str(i), ") = ", num2str(cond_A)));
    end
i = 1000;
cond_A = my_HeurMinCond(A_inf_1, A_sup_1, i);  
        disp(strcat("rad = ", num2str(matrix_radius)," : HeurMinCond(A, ", num2str(i), ") = ", num2str(cond_A)));
    iter = 100;
    for rad = 0.1:0.05:0.5
        A_inf_1 = A * (1 -  matrix_radius);
        A_sup_1 = A * (1 +  matrix_radius);
        b_sup_tmp_1 = sup_b;
        b_inf_tmp_1 = inf_b;
        cond_A = my_HeurMinCond(A_inf_1, A_sup_1, iter);        
        disp(strcat("rad = ", num2str(rad), " : HeurMinCond(A, ", num2str(iter), ") = ", num2str(cond_A)));
    end
   
end
%%
%вычисление IVE интервальной матрицы А
if(is_need_plot(17))
    %радиус элементов матрицы А - 10% от их величины
    A_inf_1 = A * 0.9;
    A_sup_1 = A * 1.1;
    b_sup_tmp_1 = sup_b;
    b_inf_tmp_1 = inf_b;
    
    cond_A = my_HeurMinCond(A_inf_1, A_sup_1);
    
    [tolmax_1, argmax_1, envs_1, ccode_1] =  tolsolvty(A_inf_1, A_sup_1, b_inf_tmp_1,  b_sup_tmp_1);
    if(tolmax_1 < 0)
        shift_1 = -tolmax_1;
        b_sup_tmp_1 = b_sup_tmp_1 + shift_1;
        b_inf_tmp_1 = b_inf_tmp_1 - shift_1;
        [tolmax_1, argmax_1, envs_1, ccode_1] =  tolsolvty(A_inf_1, A_sup_1, b_inf_tmp_1,  b_sup_tmp_1);
    end
    
    A_IVE_1 = my_IVE(A_inf_1, A_sup_1,  b_inf_tmp_1,  b_sup_tmp_1, tolmax_1, argmax_1, length(argmax_1));
    
end
%%
%--------------------------------------------------------------------------
% lab 6
%--------------------------------------------------------------------------
%%

inf_b;
sup_b;
b = (sup_b + inf_b)/2;
A = hord_matrix;


% N - число строк (256)
% M - число столбцов (174)
sizes_A = size(A);
N = sizes_A(1);
M = sizes_A(2);

% подготавливаем материал для ЗЛП
% создаём столбец где первые num_of_row нулей, далее num__of_col едениц
e = [zeros(M, 1); ones(N, 1)];
rad = ones(1, N);
diag_rad_b = diag(rad);
C = [A, -diag_rad_b
    -A, -diag_rad_b];
d = [b; -b];
lb = zeros(N+M, 1);
xw = linprog(e,C,d,[],[],lb);
x = xw(1:M, :);
w = xw(M+1:M+N, :);

%гистограмма решения x ЗЛП
if(is_need_plot(18))
    figure();
    hold on;
    grid on;
    histogram(x);
    title("Гистограмма решения ЗЛП");
    ylabel("numder")
    xlabel('x');
end

%график x решения ЗЛП
if(is_need_plot(19))
    x_axis = 1:M;
    figure();
    hold on;
    grid on;
    plot(x_axis, x, 'bo');
    title("Значение X");
    xlabel("i")
    ylabel('x');
end

%график w решения ЗЛП
if(is_need_plot(20))
    x_axis = 1:N;
    figure();
    hold on;
    grid on;
    plot(x_axis, w, 'bo');
    title("Значение W");
    xlabel("i")
    ylabel('w');
end

%сумма значений w
e = ones(1, N);
sum_w = e * w;
disp(strcat("функция цели: f(z) = ", num2str(sum_w)));

%график значений решения ЗЛП
if(is_need_plot(21))
    figure()
    hold on;
    grid on;    
    plot_A = A * x;
    x_axis = 1:length(plot_A);
    plot(x_axis, plot_A');
    plot(x_axis, inf_b');
    plot(x_axis, sup_b');
    legend("Ax", "inf b", "sup b")
    ylabel("value")
    xlabel("i")
    title("Значений решения ЗЛП");
end

