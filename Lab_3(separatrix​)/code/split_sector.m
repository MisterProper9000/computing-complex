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