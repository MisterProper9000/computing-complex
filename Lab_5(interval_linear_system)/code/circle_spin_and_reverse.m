function [RBDRY, ZBDRY] = circle_spin_and_reverse(RBDRY,ZBDRY,NBDRY, split_index)

%������� "�����"
new_RBDRY = RBDRY(split_index:NBDRY-1);
new_ZBDRY = ZBDRY(split_index:NBDRY-1);

%������� "������"
new_RBDRY = [new_RBDRY, RBDRY(1:split_index - 1)];
new_ZBDRY = [new_ZBDRY, ZBDRY(1:split_index - 1)];

%�������� ������ � �����
new_RBDRY = [new_RBDRY, RBDRY(split_index)];
new_ZBDRY = [new_ZBDRY, ZBDRY(split_index)];

% ������� ����� �� ���������������
RBDRY = rot90(new_RBDRY, 2);
ZBDRY = rot90(new_ZBDRY, 2);
end