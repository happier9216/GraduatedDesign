function [ImgRec,diff_Cost,Critical] = GSR_Solver_Denoising(ImgInput,ImgNoise, Opts,iter)
%%
if ~isfield(Opts,'PatchSize')
    Opts.PatchSize = 7;
end

if ~isfield(Opts,'SlidingDis')
    Opts.SlidingDis = 4;
    Opts.Factor = 240;
end

% if ~isfield(Opts,'ArrayNo')
%   %  Opts.ArrayNo = 60;
%   Opts.ArrayNo = 60;
% end

if ~isfield(Opts,'SearchWin')
    Opts.SearchWin = 20;
end


[Hight Width]   =   size(ImgInput);
SearchWin = Opts.SearchWin;
PatchSize    =    Opts.PatchSize;
PatchSize2    =   PatchSize*PatchSize;
ArrayNo   =   Opts.ArrayNo;
SlidingDis = Opts.step;

%tau =  Opts.lambda*Opts.Factor/Opts.mu;
%Threshold = sqrt(2*tau);
%Threshold = Opts.Threshold;

N     =  Hight-PatchSize+1;
M     =  Width-PatchSize+1;
L     =  N*M;

Row     =  [1:SlidingDis:N];
Row     =  [Row Row(end)+1:N];
Col     =  [1:SlidingDis:M];
Col    =  [Col Col(end)+1:M];
%%
PatchSet     =  zeros(PatchSize2, L, 'single');
%---add---
Patch_orgSet     =  zeros(PatchSize2, L, 'single');
%---------
x_org = Opts.org;
Count     =  0;
for i  = 1:PatchSize
    for j  = 1:PatchSize
        Count    =  Count+1;
        Patch  =  ImgInput(i:Hight-PatchSize+i,j:Width-PatchSize+j);
        Patch  =  Patch(:);
        PatchSet(Count,:) =  Patch';
        
        NoiPatch  =  ImgNoise(i:Hight-PatchSize+i,j:Width-PatchSize+j);
        NoiPatch  =  NoiPatch(:);
        NoiPatch_Set(Count,:) =  NoiPatch';
        
        Patch_org  =  x_org(i:Hight-PatchSize+i,j:Width-PatchSize+j);
        Patch_org  =  Patch_org(:);
        Patch_orgSet(Count,:) =  Patch_org';
        
    end
end

PatchSetT  =   PatchSet';
NoiseSetT  = NoiPatch_Set';
Patch_orgSetT =  Patch_orgSet';


%%
%-----------------------------¼ÆËãÔëÉù-------------------------------------------------------------
Opts.lam = 0.56;
Opts.sigma = 35;
SigmaNoi = zeros(1,(Hight-PatchSize+1)*(Hight-PatchSize+1));
if iter ==1
     SigmaNoi = Opts.sigma * ones(size(SigmaNoi));  
else
    Noi= sqrt(abs(repmat(Opts.sigma^2,1,size(PatchSetT,1))-mean((NoiseSetT-PatchSetT).^2,2)'));
    SigmaNoi = Opts.lam * Noi;
end
%--------------------------------------------------------------------------------------------------
%%
I        =   (1:L);
I        =   reshape(I, N, M);
NN       =   length(Row);
MM       =   length(Col);

ImgTemp     =  zeros(Hight, Width);
ImgWeight   =  zeros(Hight, Width);
IndcMatrix  =  zeros(NN, MM, ArrayNo);
PatchArray  =  zeros(PatchSize, PatchSize, ArrayNo);
%%
%tic;
diff_Cost = zeros(NN*MM,1);
Critical = zeros(NN*MM,1);
for  i  =  1 : NN
    for  j  =  1 : MM
        
        CurRow      =   Row(i);
        CurCol      =   Col(j);
        Off      =   (CurCol-1)*N + CurRow;
        
        
        CurPatchIndx  =  PatchSearch(PatchSetT, CurRow, CurCol, Off, ArrayNo, SearchWin, I);
        %------add-------------------------------------------------------------------------
        CurPatch_OrgIndx  =  PatchSearch(Patch_orgSetT, CurRow, CurCol, Off, ArrayNo, SearchWin, I);
        %-------------------------------------------------------------------------------
        
        
        %IndcMatrix(i,j,:) = CurPatchIndx;
        
        CurArray = PatchSet(:, CurPatchIndx);
        CurArray_Org = Patch_orgSet(:, CurPatch_OrgIndx);
        
        %----add------
        diff_Cost((i-1)*NN+j) = sum(sum((CurArray -repmat(PatchSetT(Off,:)',1,size(CurArray,2))).^2));
        %----------
        
        CurArray_Mean  =   repmat(mean( CurArray, 2 ),1,ArrayNo);
        CurArray_res    =   CurArray  -   CurArray_Mean;     
        
       [U,SigmaY,V]  = svd(CurArray_res,'econ');
    %  [U,SigmaY,V]  = svd(CurArray,'econ');
        NSig =  SigmaNoi(Off);
        Temp   =   sqrt(max( diag(SigmaY).^2 - ArrayNo*NSig^2, 0 ));
        C = 2.8284;
        for ite =1:3
            W_Vec    =  double( (C*sqrt(ArrayNo)*NSig^2)./( Temp + eps ));               % Weight vector
            SigmaX   =  soft(SigmaY,diag(W_Vec));
            Temp     = diag(SigmaX);
        end
          %    CurArray_post =  U*SigmaX*V';
                CurArray_post =  U*SigmaX*V' +  CurArray_Mean; 
        
        
        
%         SG_Z = SG_V.*(abs(SG_V)>Threshold);
%         non_zero = length(find(SG_Z>0));
%         CurArray = SG_S*SG_Z*SG_D';
        
        %--add-----------------------------------------
        Critical((i-1)*NN+j) = sum(sum((CurArray_post - CurArray_Org).^2));
        %-------------------------------------------        
        Critical_1((i-1)*NN+j) =  sum(sum((CurArray - CurArray_Org).^2));
        
        for k = 1:ArrayNo
            PatchArray(:,:,k) = reshape(CurArray_post(:,k),PatchSize,PatchSize);
        end
        
        for k = 1:length(CurPatchIndx)
            RowIndx  =  ComputeRowNo((CurPatchIndx(k)), N);
            ColIndx  =  ComputeColNo((CurPatchIndx(k)), N);
            ImgTemp(RowIndx:RowIndx+PatchSize-1, ColIndx:ColIndx+PatchSize-1)    =   ImgTemp(RowIndx:RowIndx+PatchSize-1, ColIndx:ColIndx+PatchSize-1) + PatchArray(:,:,k)';
            ImgWeight(RowIndx:RowIndx+PatchSize-1, ColIndx:ColIndx+PatchSize-1)  =   ImgWeight(RowIndx:RowIndx+PatchSize-1, ColIndx:ColIndx+PatchSize-1) + 1;
        end
        
    end
end
%%
%save ('IndcMatrix.mat', 'IndcMatrix');
ImgRec = ImgTemp./(ImgWeight+eps);

%toc;

return;



