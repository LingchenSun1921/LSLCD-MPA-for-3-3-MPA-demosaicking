function [val] = cpsnr1(a,b,border)
[n,m,~]=size(a);
X1(:,:,1)=a;
X2(:,:,1)=b;
% X1(:,:,3)=c;
% X1(:,:,4)=d;
% X2(:,:,1)=a1;
% X2(:,:,2)=b1;
% X2(:,:,3)=c1;
% X2(:,:,4)=d1;

MES=sum(sum((X1(1+border:n-border,1+border:m-border,:)-X2(1+border:n-border,1+border:m-border,:)).^2))*1^2/((n-2*border)*(m-2*border));     %¾ù·½²î
averageMES=sum(MES(:));

val=10*log10(1*1/averageMES);