%  ����ÿһ�ε��� ÿһ������Ӧ����ֵ 
%����20�ε��� 10������ ��Th Ϊ10*20  �����е�ֵΪ��ֵ
  load('part.mat');
  load('Threshold.mat')
  Th = zeros(10,13);
  for par  = 1:10
    for iter = 1:13
         ParIndex =  find(part(iter,:)==par);
         Thresh = T(iter,ParIndex);
         for j = 1:15 %��ֵ����
             num = find (Thresh ==(280+(j-1)*5));
             s(j) = size(num,2);
         end
         [value,index] = max(s);
         Th(par,iter) = 280+(index-1)*5;
    end
  end
 
%  Th(1,:) =190;
  save('Th.mat','Th');