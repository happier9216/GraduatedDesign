x=[1.2 3.4 3.1;1.3 4.2 3.6;1.9 2.6 0.5];
sigma = 10;
y = x +sigma*randn(3,3);
E= svd(y);
% Threshold = 
% E1 = E .*(abs(E)>Threshold);
% 
%  y=svd(diag([[1.7 2.5] zeros(1,98)])+randn(100)/sqrt(100));