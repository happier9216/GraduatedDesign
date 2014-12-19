A=imread('house.png');
imshow(A) %A is unit8(0,255)

E = double(A) + randn(size(A)) * 50; 
figure;imshow(E./255);
E = dct2(F);
E=log(abs(F));
figure;imshow(E);
colormap(jet(64)); %显示为64级灰度
colorbar; %显示颜色条，显示变换后的系数分布


C=dct2(A); %进行余弦变换
figure;
B=log(abs(C));
imshow(B)
colormap(jet(64)); %显示为64级灰度
colorbar; %显示颜色条，显示变换后的系数分布
C(abs(C)<10)=0; %将DCT变换后的系数值小于10的元素设为0
%E=idct2(C);
D=idct2(C)./255; %对DCT变换值归一化，进行余弦反变换???
figure;
imshow(D) ;
% imshow(uint8(E)); is the same as D=idct2(C)./255
% imshow(E,[]); is the same as D=idct2(C)./255