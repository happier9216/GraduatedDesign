A=imread('house.png');
imshow(A) %A is unit8(0,255)

E = double(A) + randn(size(A)) * 50; 
figure;imshow(E./255);
E = dct2(F);
E=log(abs(F));
figure;imshow(E);
colormap(jet(64)); %��ʾΪ64���Ҷ�
colorbar; %��ʾ��ɫ������ʾ�任���ϵ���ֲ�


C=dct2(A); %�������ұ任
figure;
B=log(abs(C));
imshow(B)
colormap(jet(64)); %��ʾΪ64���Ҷ�
colorbar; %��ʾ��ɫ������ʾ�任���ϵ���ֲ�
C(abs(C)<10)=0; %��DCT�任���ϵ��ֵС��10��Ԫ����Ϊ0
%E=idct2(C);
D=idct2(C)./255; %��DCT�任ֵ��һ�����������ҷ��任???
figure;
imshow(D) ;
% imshow(uint8(E)); is the same as D=idct2(C)./255
% imshow(E,[]); is the same as D=idct2(C)./255