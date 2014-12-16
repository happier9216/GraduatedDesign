%part 将每一次迭代每一个块不同差别Cost分成10个part
%得到part矩阵中为 IterNum * patchnum  矩阵中为1-10
load('CostCritical.mat','Total_Cost_T');
temp = zeros(1,10);
point = zeros(13,10);
part = zeros(3969,13); 
for Iter = 1:13
    Top = max(Total_Cost_T(:,Iter));
    Low = min(Total_Cost_T(:,Iter));
    
    step = (Top - Low)/9;
    temp = [Low:step:Top];
    point(Iter,: ) = temp(:);
    
    for patch = 1:3969
        if Total_Cost_T(patch,Iter)<= point(Iter,2)
            part(patch,Iter) = 1;
        elseif Total_Cost_T(patch,Iter)<=point(Iter,3)
            part(patch,Iter) = 2;
        elseif Total_Cost_T(patch,Iter)<=point(Iter,4)
            part(patch,Iter) = 3;
        elseif Total_Cost_T(patch,Iter)<=point(Iter,5)
            part(patch,Iter) = 4;
        elseif Total_Cost_T(patch,Iter)<=point(Iter,6)
            part(patch,Iter) = 5;
        elseif Total_Cost_T(patch,Iter)<=point(Iter,7)
            part(patch,Iter) = 6;
        elseif Total_Cost_T(patch,Iter)<=point(Iter,8)
            part(patch,Iter) = 7;
        elseif Total_Cost_T(patch,Iter)<=point(Iter,9)
            part(patch,Iter) = 8;
        elseif Total_Cost_T(patch,Iter)<=point(Iter,10)
            part(patch,Iter) = 9;
        else
            part(patch,Iter) = 10;      
        end
        
    end
       
end
part = part';
save('part.mat','part');



