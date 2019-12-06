function [] = Element_Draw(elem, color, num, vertex_flag)

N = length(elem.points);

x = zeros(N+1, 1);
y = zeros(N+1, 1);
for i=1:N
    x(i) = elem.points(i).x;
    y(i) = elem.points(i).y;
end
txt_pos_x = sum(x) / (N);
txt_pos_y = sum(y) / (N);
x(N+1) = elem.points(1).x;
y(N+1) = elem.points(1).y;

n=nargin;
if(n == 1)
    plot(x, y);
elseif(n == 2)
    plot(x, y, color);
elseif(n == 3)
    plot(x, y, color);
    %str_num = strcat(" (", num2str(num), ")");
    str_num = "";
    text(txt_pos_x, txt_pos_y, strcat(num2str(elem.index),str_num));
elseif(n == 4)
    plot(x, y, color);
    %str_num = strcat(" (", num2str(num), ")");
    str_num = "";
    text(txt_pos_x, txt_pos_y, strcat(num2str(elem.index),str_num));
    if(vertex_flag)
        for i=1:N
            text(x(i), y(i), num2str(i))
        end
    end
end

end