%  ��ÿһ�ε���ÿһ������ԭͼ���Critical 
%�ҵ�ÿһ�ε����Ǹÿ���ԭͼ�����С������ֵT
%�õ�T������Ϊ IterNum * patchnum  �����иôε����ÿ��������ֵ
load('CostCritical.mat','Total_Critical_T');
  IterNum = 13;
  PatchNum = 63*63;
  T = zeros(IterNum,PatchNum);
  
  for IterNums = 1:IterNum
    for patchIndex = 1:PatchNum

        [value,j] =min(Total_Critical_T(patchIndex,IterNums,:));
        T(IterNums,patchIndex) =  280+(j -1)*5;
    
    end     
  end
  save('Threshold.mat','T');