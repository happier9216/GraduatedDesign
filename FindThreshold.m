%  将每一次迭代每一个块与原图差距Critical 
%找到每一次迭代是该块与原图差距最小差距的阈值T
%得到T矩阵中为 IterNum * patchnum  矩阵中该次迭代该块的最优阈值
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