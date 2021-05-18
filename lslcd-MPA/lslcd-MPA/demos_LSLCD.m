% Least-Squares Luma-Chroma Demultiplexing Algorithm for Demosaicking
% IEEE Trans. Image Processing (submitted)
% This software is for provided for non-commercial and research purposes only
% Copyright Eric Dubois (edubois@site.uottawa.ca), University of Ottawa 2010 
% **************Q1=imfilter(A,h,'conv','same');
% a=1;
% A1_1=[];
% for i=1:5
%     for j=1:5
%         B=A1(i:i+2,j:j+2);
%         A1_1=[A1_1;(B(:))'];
%     end
% end
%       [RGB] = demos_freq_adapt(CFA, fpkg)
%       CFA(:,:)   input Bayer CFA mosaicked image, double according to the
%       pattern       G R
%                     B G
%       fpkg       path of the filter package that contains the follwing
%                  filter: h1, h2a, h2b, hG1, hG2
%       RGB(:,:,3) output RGB image, double

function [S0,S1,S2,a1,b1,c1]= demos_LSLCD(CFA,fpkg)

% Load the filter package

load('fpkg.mat');
% Use the filters to demosaic the CFA image
S = size(CFA); N1 = S(1); N2 = S(2);
yc = 0:N1-1; xc = 0:N2-1;
[XC,YC] = meshgrid(xc,yc);
% 
% F=fft2(CFA);
% % F=abs(F);
% 
% PQ=paddedsize(size(F));
% Fp = fft2(CFA,PQ(1),PQ(2));
% % gausFilter1=fspecial('gaussian',[445,493],0.7);
% % gausFilter1=fspecial('gaussian',[PQ(1),PQ(2)],0.7);
% % x1=fft2(gausFilter1);
%  Hp =lpfilter('gaussian',PQ(1),PQ(2),(PQ(1))/6);%生成一个PQ大小的高斯滤波器
% % x1=circshift(x1,floor(1/2*size(x1)));
% Lhat=Fp.* Hp;
% Lhat=ifft2(Lhat);
% Lhat = Lhat(1:size(CFA,1),1:size(CFA,2));
% % Lhat=ifftshift(Lhat);
% imshow(Lhat,[]);
% % gausFilter1=ifft2(x1);
% % Lhat=imfilter(CFA,gausFilter1,'replicate','same');
% % gausFilter2=fspecial('gaussian',[445,493],1.0);
% F=fftshift(Fp);
% 
% gausFilter2 =lpfilter('gaussian',PQ(1),PQ(2),50);
% x1=circshift(gausFilter2,floor(1/3*size(gausFilter2)));
% C1hat=F.*x1;
% % C1hat=ifftshift(C1hat);
% C1hat=ifft2(C1hat);
% C1hat = C1hat(1:size(CFA,1),1:size(CFA,2));
% % gausFilter2=ifft2(x1);
% % C1hat=imfilter(CFA,gausFilter2,'replicate','same');
% gausFilter3 =lpfilter('gaussian',PQ(1),PQ(2),50);
% x1=circshift(gausFilter3,floor(-1/3*size(gausFilter3)));
% C2hat=F.*x1;
% % C2hat=ifftshift(C2hat);
% C2hat=ifft2(C2hat);
% C2hat = C2hat(1:size(CFA,1),1:size(CFA,2));
% Extract chrominance in corners using h1
% C1hatR= imfilter(CFA,h1R,'replicate','same');
% C1hatI= imfilter(CFA,h1I,'replicate','same');
% imshow(C1ahatR,[]);
C1hatR= imfilter(CFA,h1R,'replicate','same');

C1hatI= imfilter(CFA,h1I,'replicate','same');
C2hatR= imfilter(CFA,h2R,'replicate','same');
% imshow(C1ahatR,[]);
C2hatI= imfilter(CFA,h2I,'replicate','same');
% imshow(C1ahatI,[]);


% imshow(C2ahatI,[]);
% Extract chrominance on sides at f_y = 0
Lhat = CFA -complex(C1hatR,C1hatI)-complex(C2hatR,C2hatI);

C1ahat=complex(C1hatR,C1hatI);
C2ahat=complex(C2hatR,C2hatI);

X1=exp(complex(0,2*pi*((1/3)*XC+(1/3)*YC)));
X2=exp(complex(0,2*pi*(-(1/3)*XC-(1/3)*YC)));
% imshow(C1ahat,[]);

% C1hat=C1ahat.*conj(X1);
C1hat=complex((real(C1ahat.*conj(X1))+real(C2ahat.*conj(X2)))/2,(imag(C1ahat.*conj(X1))-imag(C2ahat.*conj(X2)))/2);
C2hat=conj(C1hat);

a1=1.0003*Lhat+ 1*C1hat+1*C2hat;
b1=  1.*Lhat+complex(-0.5, - 0.866)*C1hat+complex(-0.5, 0.866)*C2hat;
c1 = 1.*Lhat+complex(-0.5, 0.866)*C1hat+complex(-0.5, - 0.866)*C2hat;
S0=2*Lhat;
S2=complex(0,-2)*C1hat+complex(0,2)*C2hat;
S1=2*C1hat+2*C2hat;