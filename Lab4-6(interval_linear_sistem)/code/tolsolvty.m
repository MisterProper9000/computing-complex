function [tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb,varargin) 
%  
%   ���������� ��������� ������������� ����������� ����������� ��������� 
%   ������� ��� ������������ ������� �������� �������������� ���������. 
% 
%   TOLSOLVTY(infA, supA, infb, supb) ����� �������� ��������� ������������� 
%   ����������� ����������� ��������� ������� ��� ������������ ������� �������� 
%   ��������� Ax = b, � ������� ������� ������ � ������� ������ ��������� A  
%   ����� infA � supA, � ������� ������ � ������� ������ ������ ����� b  ����� 
%   infb  �  supb  ��������������. �������������  ���������  �������  �������� 
%   ��������� - ���������� ������� ��� ������������� ������������ �������� 
%   ������� Ax = b, ������� ���������� ���� ������������, � �����  ����������  
%   � �������/��������� ����������� ��������� ������� � ����������� ������. 
%  
%   ��������� ������:
%       [tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb, ... 
%                                           iprn,weight,epsf,epsx,epsg,maxitn) 
%  
%   ������������ ������� ��������� �������: 
%        infA, supA - ������� ����� � ������ ������ ������������ ������������� 
%                     ���  �����������  ���  ������������  �������  �������� 
%                     �������������� ���������; ��� ����� ���� ��������������, 
%                     �� ������ ����� ���������� �������; 
%        infb, supb - ������� ����� � ������ ������ ����������  ������ ����� 
%                     ������������ ������� �������� �������������� ���������. 
%  
%   �������������� ������� ��������� �������:
%              iprn - ������ ��������� ������; ���� iprn > 0 - ���������� 
%                     � ���� �������� ���������� ����� ������ iprn-��������;
%                     ���� iprn <= 0 (�������� �� ���������), ������ ���;
%            weight - ������������� ������ ������� ������������� ��� ���������� 
%                     ������������� �����������, �� ��������� ������ ������ 
%                     ������� �� ����� ���������� ������������; 
%              epsf - ������ �� �������� �� �������� �������� �����������,
%                     �� ��������� ��������������� 1.e-6;
%              epsx - ������ �� �������� �� ��������� �������� �����������,
%                     �� ��������� ��������������� 1.e-6;
%              epsg - ������ �� ������� ����� �������������� �����������,
%                     �� ��������� ��������������� 1.e-6;
%            maxitn - ����������� �� ���������� ����� ���������, 
%                     �� ��������� ��������������� 2000.
%  
%   �������� ��������� �������: 
%            tolmax - �������� ��������� ������������� �����������;
%            argmax - ������������ ��� ������ �������� ���������,������� 
%                     ����� � ���������� ��������� ������� ��� tolmax>=0;
%              envs - �������� ���������� ������������� ����������� � ����� 
%                     ��� ���������, ��������������� �� �����������; 
%             ccode - ��� ���������� ��������� (1 - �� ������� epsf �� 
%                     ��������� �������� �����������, 2 - �� ������� epsg 
%                     �� �������������, 3 - �� ������� epsx �� �������� 
%                     ���������, 4 - �� ����� ��������, 5 - �� ������ 
%                     �������� �� �����������). 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%    ���  ���������  ���������  ������������ ����������� ��������� ������� 
%    ��� ������������ �������  ��������  ��������������  ���������  Ax = b 
%    � ������������ �������� A = [infA, supA] � ������������ �������� ������ 
%    �����  b = [infb, supb] � ������� ������������ ������������� ����������� 
%    ����������� ��������� ������� ���� �������. ��. ����������� � 
%   
%       ����� �.�. ������������� ������������ ������. - �����������: XYZ, 
%       2017. - ����������� �����, ��������� �� http://www.nsc.ru/interval,
%       �������� 6.4;
%       Shary S.P. Solving the linear interval tolerance problem //
%       Mathematics and Computers in Simulation. - 1995. - Vol. 39.
%       - P. 53-85.
%   
%   ������ ������� ������������� ��� ���������� ������������� ����������� 
%   (����������, ��� ��������� ��������� �������) ��������� ��������� ������ 
%   �������� ��������� ���������, ��������������� ������������� ���������� 
%   � ������ �������������� ������������ � �.�. 
%
%   ���  ������������  ���������  �������������  �����������  ������������ 
%   ������� ��������� ����������������� ������� � ����������� ������������ 
%   � �����������  ��������  ���������������� ���������������, ������������
%   (��� ������ �����������) � ������ 
%       ��� �.�., �������� �.�. ����� �����������, ������������ ��������
%       ���������� ������������ � ����������� �������� ���� ���������������� 
%       ���������� // �����������. - 1971. - �3. - �. 51-59. 
%   
%   � �������� ������ ���� ����� ��������� ������������ ��������� ��������� 
%   ����������� ralgb5, ������������� � ������������� �.�.�������� (�������� 
%   ����������� ��� �������, ����). �������� ���� �������� ������ � ������ 
% 
%       ������ �.�. �������������� ������ ralgb5 � ralgb4 ��� ����������� 
%       �������� �������� ������� // �������������� ����������. - 2017. - 
%       �. 22, � 2. - �. 127-149. 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%   �.�. �����, ��� �� ��� & ���, 2007-2019 ��. 
%   �.�. ���������, ������, 2019 �.  
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%   
%   �������� ������������ ������� ������ 
%   
  
mi = size(infA,1);  ni = size(infA,2);  
ms = size(supA,1);  ns = size(supA,2); 
if mi==ms   %   m - ���������� ��������� � ������� 
    m = ms; 
else 
    error('���������� ����� � �������� ����� � ������ ������ �����������')
end
if ni==ns 
    n = ns; %   n - ���������� ����������� ���������� � ������� 
else 
    error('���������� �������� � �������� ����� � ������ ������ �����������')
end 
  
ki = size(infb,1); 
ks = size(supb,1); 
if ki==ks 
    k = ks; 
else 
    error('���������� ��������� � �������� ����� � ������ ������ �����������')
end
if k~=m 
    error('������� ������� ������� �� ������������� �������� ������ �����') 
end
  
if ~all(all(infA <= supA)) 
    error('� ������� ������� ����� ������������ ������������ �������') 
end 
  
if ~all(infb <= supb) 
    error('� ������� ������ ����� ������ ������������ ������������ ����������') 
end 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%   ������� ���������� ��������� ����������������� ������� � ������ 
%
maxitn = 2000;          %   ����������� �� ���������� ����� ���������
nsims  = 30;            %   ���������� ���������� ���������� �����
epsf = 1.e-6;           %   ������ �� ��������� �������� ����������� 
epsx = 1.e-6;           %   ������ �� ��������� ��������� �����������
epsg = 1.e-6;           %   ������ �� ����� �������������� �����������
  
alpha = 2.3;            %   ����������� ���������� ������������ � ���������
hs = 1.;                %   ��������� �������� ���� ����������� ������
nh = 3;                 %   ����� ���������� ����� ����������� ������ 
q1 = 0.9;               %   q1, q2 - ��������� ���������� �����������
q2 = 1.1;               %       �������� ���������
  
iprn = 0;               %   ������ � ���� �������� ����� ������ iprn-��������
                        %   (���� iprn < 0, �� ������ �����������) 
weight = ones(m,1);     %   ������� ������� ������� ������������� ��� ���������� 
format short g;         %   ������ ������ ������  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ������������ ��������� �������� ��� ���������� ��������� ������ 
%
HorLine = '-------------------------------------------------------------';
TitLine = '�������� ������������ ������������� ����������� Tol';
TabLine = '���        Tol(x)         Tol(xx)   ������/���  ������';
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   �������������� ���������� ���������, �������� ������������� 
% 
if nargin >= 5 
    iprn = ceil(varargin{1}); 
    if nargin >= 6 
        weight = varargin{2}; 
        if size(weight,1)~=m 
            error('������ ������� ������� ������������� ����� �����������') 
        end 
        if any( weight <= 0 ) 
            error(' ������ ������� ������������� ������ ���� �������������') 
        end 
        if nargin >= 7 
            epsf = varargin{3}; 
            if nargin >= 8 
                epsx = varargin{4}; 
                if nargin >= 9 
                    epsg = varargin{5}; 
                    if nargin >= 10 
                        maxitn = varargin{6}; 
                    end
                end
            end
        end 
    end 
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function [f,g,tt] = calcfg(x)
% 
%   �������, ������� ��������� �������� f ���������������� ������������� 
%   ����������� � ��� ������������� g;  ����� ����, ��� ����� ������ tt 
%   �� �������� ���������� ����������� � ������ ����� ��������� 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ��������������� ���������� ������� �������� 
  
infs = zeros(m,1); 
sups = zeros(m,1); 
tt = zeros(m,1); 
dl = zeros(n,1); 
ds = zeros(n,1); 
dd = zeros(n,m); 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   ��������� �������� ������������� ����������� � ������� dd, 
    %   ������������ �� ��������������� ��� ���������� 
    
    for i = 1:m 
        los = 0.5*(infb(i)+supb(i));
        sus = los; 
        for j = 1:n
            if x(j) >= 0  
                los = los - supA(i,j)*x(j); 
                sus = sus - infA(i,j)*x(j); 
            else 
                los = los - infA(i,j)*x(j); 
                sus = sus - supA(i,j)*x(j); 
            end 
        end 
        %   ������ � ������� ����� ���������, ������� ����������  
        %   ��� ������� � ��������� ��� i-�� ���������� 
        infs(i) = los; 
        sups(i) = sus; 
    end 
    
    %   ���������� �������� i-�� ���������� ������������� 
    %   ����������� � � �������������� 
    for i = 1:m 
        alos = abs(infs(i));  asus = abs(sups(i)); 
        %   ���������� �������������� dl ������� �����
        for j = 1:n 
            dm = infA(i,j);
            dp = supA(i,j);
            if x(j) < 0
                dl(j) = dm;
            else 
                dl(j) = dp;
            end
        end
        %   ���������� �������������� ds �������� ����� 
        for j = 1:n 
            dm = supA(i,j);
            dp = infA(i,j);
            if x(j) < 0 
                ds(j) = dm;
            else
                ds(j) = dp;
            end
        end
        %   ������ ������� �������������� i-�� ����������
        if alos ~= asus
            if alos < asus 
                mags = asus;  dd(:,i) = weight(i)*ds; 
            else 
                mags = alos;  dd(:,i) = -weight(i)*dl; 
            end 
        else
            mags = alos; 
            if sups(i) > 0 
                dd(:,i) = weight(i)*ds; 
            else 
                dd(:,i) = -weight(i)*dl; 
            end
        end 
        %   ���������� � ����������� �������� i-�� ���������� 
        tt(i) = weight(i)*(0.5*(supb(i) - infb(i)) - mags);
    end  
   
    %   �������� ����������� �� �������� ����������
    %   � ������������ ����� ������������� 
    [f,mc] = min(tt);
    g = dd(:,mc);
  
end 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   ��������� ��������� ����������� x ��� ������� ���� ������������� 
%   '�������' �������� �������, ���� ��� �� ������� ����� �����������,
%   ����� ���� ��������� ������������ ������� ������ 
% 
Ac = 0.5*(infA + supA); 
bc = 0.5*(infb + supb);
sv = svd(Ac);
minsv = min(sv);
maxsv = max(sv);
  
if ( minsv~=0 && maxsv/minsv < 1.e+15 ) 
    x = Ac\bc; 
else
    x = zeros(n,1);
end  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ������� �������:
%       B - ������� ��������� �������������� ������������
%       vf - ������ ���������� ����������� �� ��������� ����� ���������
%       g, g0, g1 - ������������ ��� �������� ��������������� ��������,
%           �������������� ��������������� ����������� � ��.

B = eye(n,n);                   %   �������������� ��������� �������� 
vf = realmax*ones(nsims,1);     %   �������������� ������ �������� ������� 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   ��������� ��������� ����������
  
w = 1./alpha - 1.;
lp = iprn; 
  
[f, g0, tt] = calcfg(x); 
ff = f ;  xx = x;
cal = 1;  ncals = 1; 
  
if iprn > 0 
    fprintf('\n\t%52s\n',TitLine); 
    fprintf('%65s\n',HorLine); 
    fprintf('\t%50s\n',TabLine); 
    fprintf('%65s\n',HorLine); 
    fprintf('\t%d\t%f\t%f\t%d\t%d\n',0,f,ff,cal,ncals); 
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   �������� ���� ���������: 
%       itn - ������� ����� ��������
%       xx  - ����������� � ��������� ��������� �����������
%       ff  - ����������� � ��������� ����������� 
%       cal - ���������� ���������� ����������� �� ������� ���� 
%     ncals - ����� ���������� ���������� �������� �����������
%  
for itn = 1:maxitn;
    vf(nsims) = ff;
    %   �������� �������� �� ����� ��������������
    if  norm(g0) < epsg
        ccode = 2;  
        break
    end
    %   ��������� ������������� � ��������������� ������������,
    %   ���������� ����������� ������� 
    g1 = B' * g0;
    g = B * g1/norm(g1); 
    normg = norm(g);
    %   ���������� ������ �� ����������� g:
    %       cal - ������� ����� ����������� ������,
    %       deltax - �������� ��������� � �������� ������
    r = 1;  
    cal = 0;  
    deltax = 0;
    while ( r > 0. && cal <= 500 )
        cal = cal + 1; 
        x = x + hs*g; 
        deltax = deltax + hs*normg; 
        [f, g1, tt] = calcfg(x); 
        if f > ff 
            ff = f; 
            xx = x; 
        end 
        %   ���� ������ nh ����� ����������� �������, 
        %   �� ����������� �������� ���� hs 
        if mod(cal,nh) == 0  
            hs = hs*q2; 
        end 
        r = g'*g1; 
    end 
    %   ���� �������� ����� ����� ����� ����������� �������, �� �����
    if cal > 500 
        ccode = 5; 
        break; 
    end 
    %   ���� ���������� ������ ����� ���� ���, 
    %   �� ��������� �������� ���� hs 
    if cal == 1
        hs = hs*q1;
    end 
    %   �������� ���������� � ��� ������������� ������� �
    ncals = ncals + cal;
    if itn==lp
        fprintf('\t%d\t%f\t%f\t%d\t%d\n',itn,f,ff,cal,ncals); 
        lp = lp + iprn;
    end
    %   ���� �������� ��������� � ���������� ������ ����, �� �����
    if deltax < epsx 
        ccode = 3;  
        break; 
    end
    %   ������������� ������� �������������� ������������ 
    dg = B' * (g1 - g0);
    xi = dg / norm(dg);
    B = B + w*B*(xi*xi');
    g0 = g1;
    %   �������� ��������� �������� �����������, �������������� 
    %   ���� �����������, �� ��������� nsims ����� ���������
    vf = circshift(vf,1);
    vf(1) = abs(ff - vf(1)); 
    if abs(ff) > 1
        deltaf = sum(vf)/abs(ff);
    else 
        deltaf = sum(vf);
    end
    if deltaf < epsf 
        ccode = 1;  
        break
    end 
    ccode = 4; 
end
   
tolmax = ff;
argmax = xx; 
  
%   ��������� ���������� ������������� ����������� �� ����������� 
tt = [(1:m)', tt];
[z,ind] = sort(tt(:,2));
envs = tt(ind,:);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ����� ����������� ������ 
   
if iprn > 0
    if rem(itn,iprn)~=0
        fprintf('\t%d\t%f\t%f\t%d\t%d\n',itn,f,ff,cal,ncals); 
    end
    fprintf('%65s\n',HorLine); 
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
disp(' ');
if tolmax >= 0
    disp(' ���������� ��������� ������� ������������ �������� ������� ������� ')
else 
    disp(' ���������� ��������� ������� ������������ �������� ������� ����� ')
end 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
disp(' ');
if ( tolmax < 0. && abs(tolmax/epsx) < 10 ) 
    disp(' ���������� �������� ������������ ���������');
    disp('                          ��������� � �������� �������� ��������'); 
    disp(' ������������� ���������  � �������� ����������  epsf �/��� epsx');
    disp(' ��� ��������� ������� ���������� � ������������ ���������������'); 
    disp(' ������ � ��������');
    disp(' ');
end 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  