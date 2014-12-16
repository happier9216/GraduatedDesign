clear;
clc;
cur = cd;
addpath(genpath(cur));


for ImgNo = 1    
    for Type = 5
        switch ImgNo
            case 1
                ImgName = 'house.png';
           
        end
         x_org = double(imread(ImgName));
        
         
        switch Type
            case 1
                Img_After_BM3D = 'sig15_house34.9447.png';
                x_noi = double(imread(Img_After_BM3D));                   
            case 2
                Img_After_BM3D = 'sig25_house32.8646.png';
                x_noi = double(imread(Img_After_BM3D));   
            case 3
                Img_After_BM3D = 'sig30_house_32.087.png';
                x_noi = double(imread(Img_After_BM3D));   
            case 4
                Img_After_BM3D = 'sig35_house31.3762.png';
                x_noi = double(imread(Img_After_BM3D)); 
                
            case 5
                sigma = 35;
                randn('state', 1); % initialization
                x_noi = x_org + randn(size(x_org)) * sigma; 
               
        end
        Opts   = ParSet(sigma,x_org);  
        
        
        %Opts = [];
       % Opts.org = x_org;

        %Opts.mu = 2.5e-3*3;
        %Opts.lambda = 0.5532;
        
%         if ~isfield(Opts,'ArrayNo')
%             Opts.ArrayNo = 140;
%         end
        
      
      %  Opts.mu = 0.12;
      %  Opts.lambda = 3.8283;
     %   Opts.block_size = 8;
        %Opts.Threshold = 376;

        IterNum = Opts.Iter;
        MSE = zeros(1,IterNum+1);
        
        MSE(1) = sum(sum((x_noi-x_org).^2))/numel(x_org);
        fprintf('Initial PSNR = %f\n',csnr(x_noi,x_org,0,0));

%         mu = Opts.mu;
%         invmu = 1/(mu+1);
%         muinv = mu/(mu+1);
        
%         Threshold = 150;
%         Threshold = 192;
   
        
%         Total_Cost_T = zeros(3969,IterNum,size(Threshold,2));
%         Total_Critical_T = zeros(3969,IterNum,size(Threshold,2));
%         
%         Total_Cost_T = zeros(3969,IterNum,size(Threshold,2));
%         Total_Critical_T = zeros(3969,IterNum,size(Threshold,2));
        
               

            x = x_noi;
            w = zeros(size(x_org));
            c = zeros(size(x_org));
     
           
            for iter = 1:IterNum
                 if (mod(iter-1,Opts.Innerloop)==0)
                    Opts.ArrayNo =  Opts.ArrayNo - 10;
                 end
                
                 est = x + Opts.delta*(x_noi - x);
        
                [x,diff_Cost,Critical] = GSR_Solver_Denoising(est,x_noi,Opts,iter);
%                 Total_Cost_T(:,iter,Outloop) = diff_Cost;
%                 Total_Critical_T(:,iter,Outloop)= Critical;
                
%                 delta = 0.1;
%                 x = w + delta*(x_noi - w);
                
                
%                 x = muinv .*(w+c) + invmu .* x;
%                 c = c + (w-x);

                x_resid = est - x_org;
                MSE(iter+1) =  (x_resid(:)'*x_resid(:))/numel(x_org);
                %MSE(iter+1) = mean(x_resid2(:).^2);
                fprintf('iter number = %d, PSNR x= %f,  PSNR est= %f\n',iter,csnr(x_org,x,0,0),csnr(x_org,est,0,0));

               % x_aft  =  w;
            end
            
            PSNR_seq = 10*log10((255.^2)./(MSE));
  %    end
        
%         save('T_house_sig15_1207.mat','Total_Cost_T','Total_Critical_T');
%         save('MSE.mat','MSE');
%         save('PSNR_seq.mat','PSNR_seq');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display Picture and Results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(4);imshow(uint8(x_noi));title(['Degraded Image, PSNR = ',num2str(csnr(x_noi,x_org,0,0))]);
        figure(5);imshow(uint8(x));title(['update Image via GSR, PSNR = ',num2str(csnr(x_org,x,0,0))]);
        figure; plot(1:length(MSE), PSNR_seq, 'LineWidth',2.0),
        title('Evolution of PSNR (dB)','FontName','Times','FontSize',15),
        set(gca,'FontName','Times'),
        set(gca,'FontSize',14),
        xlabel('Iterative Numbers ');
        ylabel('PSNR');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % Compute PSNR
        PSNR = csnr(x_org,x,0,0);

        Deblurred_Name = strcat(ImgName,'_GSR','_PSNR_',num2str(PSNR),'dB','.png');
       % imwrite(uint8(est),strcat('GSR_Results\',Deblurred_Name));
        imwrite(uint8(x),strcat('GSR_Results\',Deblurred_Name));
        pause(1);
        fprintf('GSR Finished...\n')
        fprintf('***************************************************************\n')
        fprintf('***************************************************************\n')

    end    
end
