function [val] = cpsnr(a,c,e,a1,c1,e1,border)


% [n,m,ch]=size(X1);
[n,m,~]=size(a);
X1(:,:,1)=a;
X1(:,:,2)=c;
X1(:,:,3)=e;

X2(:,:,1)=a1;
X2(:,:,2)=c1;
X2(:,:,3)=e1;


% MSE_R = sum(sum((R(frame+1:N1-frame,frame+1:N2-frame)-...
%     R_ESTM(frame+1:N1-frame,frame+1:N2-frame)).^2))*255^2/((N1-2*frame)*(N2-2*frame));
% MSE_G = sum(sum((G(frame+1:N1-frame,frame+1:N2-frame)-...
%     G_ESTM(frame+1:N1-frame,frame+1:N2-frame)).^2))*255^2/((N1-2*frame)*(N2-2*frame));
% MSE_B = sum(sum((B(frame+1:N1-frame,frame+1:N2-frame)-...
%     B_ESTM(frame+1:N1-frame,frame+1:N2-frame)).^2))*255^2/((N1-2*frame)*(N2-2*frame));
% diff=X1(1+border:n-border,1+border:m-border,:)-X2(1+border:n-border,1+border:m-border,:);
MES=sum(sum((X1(1+border:n-border,1+border:m-border,:)-X2(1+border:n-border,1+border:m-border,:)).^2))*1^2/((n-2*border)*(m-2*border));     %¾ù·½²î
averageMES=sum(MES(:))/3;
% mse=mean(diff(:).^2);
% val=10*log10(255*255/mse);
val=10*log10(1*1/averageMES);
