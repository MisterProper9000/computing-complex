function [hord_dist, intersection] = Element_GetHord(elem, k, b)
POINT = Point_get_interface();
ELEMENT = Element_get_interface();
N = length(elem.points);
    
intersection=[];
tmp = [elem.points, elem.points(1)];
for i = 1:N
    A = tmp(i);
    B = tmp(i+1);
       
    if((k * A.x + b - A.y) * (k * B.x + b - B.y) > 0)
        %����������� ���
        hord_dist = 0;
        continue;
    end
    
    if( (k * A.x + b - A.y) * (k * B.x + b - B.y) == 0)
        if(k * A.x + b - A.y == 0)
            %������ ����� ����� ����� �
            intersection = [intersection, A];
        end
        %������ ����� �������� ��� ������ ��������� �������
%         if(k * B.x + b - B.y == 0)
%             %������ ����� ����� ����� B
%             intersection = [intersection, B];
%         end
        continue
    end
    koef = abs((k * A.x + b - A.y)/(k * B.x + b - B.y));
    P = POINT.create((A.x + koef * B.x) / (1 + koef), (A.y + koef * B.y) / (1 + koef));
    intersection = [intersection, P];
end

if(length(intersection) <= 1)
    %����������� ��� (���� �� ���������)
    hord_dist = 0;
    %intersection = [];
    return 
end
    

%����������� ���� (������ 2 ����)
%FIXME �������� �� ������, ����� ����� 3 (��� ��� ����)? 
%����� ��������� ������������ ���������� ����� ���������???
hord_dist =  POINT.dist(intersection(1), intersection(2));
if(length(intersection) > 2)
    ELEMENT.draw(elem, "b", -1, true)
    disp(length(intersection));
end
end