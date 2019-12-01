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