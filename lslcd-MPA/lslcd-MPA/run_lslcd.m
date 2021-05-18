% Perform demosaicking on the given experimental dataset using one single filter 

% IEEE Trans. Image Processing (submitted)
% This software is provided for non-commercial and research purposes only
% This software is modified based on the original version given by:
%  Eric Dubois (edubois@site.uottawa.ca), University of Ottawa 2010

% Copyright Li Lin(lin li@nwpu.edu.cn),NWPU (Northwestern Polytechnical University) 
% *****************

% Initialize MATLAB Workspace
clc; clear all; close all;
N_sum = zeros(8,40); 
% for SNR = 11:1:40

    for z=38:40
% Setup frame size when computing MSE/PSNR
frame = 13; 
y3331S0=[];
y3331DoLP=[];
% Select and load filter
% LSLCD algorithm using $f_{ST}$

% load('Result12.2.Test332.Test1.mat');
% name = ['Filternew\',num2str(z)];
% load(name);
% name = ['mat\',num2str(z)];
% load(name');

name = ['polar_data\',num2str(z)];
file=dir(name);
% load(name');
% II=Im3Dg1;
% II=I2;
name1 = [name,'\',file(3).name];name4 = [name,'\',file(4).name];
name2 = [name,'\',file(5).name];name3 = [name,'\',file(6).name];

% R1=RGB_0(:,:,3);R3=RGB_45(:,:,3);R5=RGB_90(:,:,3);R7=RGB_135(:,:,3); %%%%0,45,90,135
R1=imread(name1);R3=imread(name2);R5=imread(name3);R7=imread(name4); %%%%0,45,90,13
% load('filters\ST\filter_ST_01.mat');
name = ['Filternew2\',num2str(z)];
load(name);

% II=I;
% load('filters\ST\filter_ST_01.mat');
N=13;
% h111=zeros(N,N);
% h112=zeros(N,N);
h113=zeros(N,N);
h114=zeros(N,N);
h115=zeros(N,N);
h116=zeros(N,N);


% for i=1:N
% h111(1:N,i)=h1R(i*N-N+1:i*N,1);
% end
% for i=1:N
% h112(1:N,i)=h1I(i*N-N+1:i*N,1);
% end
for i=1:N
h113(1:N,i)=h1R(i*N-N+1:i*N,1);
end
for i=1:N
h114(1:N,i)=h1I(i*N-N+1:i*N,1);
end
for i=1:N
h115(1:N,i)=h2R(i*N-N+1:i*N,1);
end
for i=1:N
h116(1:N,i)=h2I(i*N-N+1:i*N,1);
end


% h1R=h111;
% h1I=h112;
h1R=h113;
h1I=h114;
h2R=h115;
h2I=h116;




% % LSLCD algorithm using $f_{TO}$
% load('filters\TO\filter_TO_01.mat');
% % LSLCD algorithm using $f_{ST_RE}$
% load('filters\ST_RE\filter_ST_RE_01.mat');
% % LSLCD algorithm using $f_{TO_RE}$
% load('filters\TO_RE\filter_TO_RE_01.mat');

% Load Gaussian filter
% load('filters\GF\GF.mat');

% Prepare the filter package
fpkg = 'fpkg.mat';
save(fpkg,'h1R','h1I','h2R','h2I');
clear h2R h2I h22R h22I;
% save(fpkg,'h1R','h1I','h2R','h2I');
% clear h1R h1I h2R h2I;



R1=im2double(R1(:,:,3));R3=im2double(R3(:,:,3));
R5=im2double(R5(:,:,3));R7=im2double(R7(:,:,3)); %%%%0,45,90,135

[R2 R4 R6 R8 ORIGS0 ORIGS1 ORIGS2]=SSCalculate(R1,R3,R5,R7);
I_disp = ORIGS0 / 2;
I_disp = imagesc(I_disp);caxis([0 1]);
colorbar;


clear Im3Dg1 II I;

ORIGDOLP=sqrt(ORIGS1.^2+ORIGS2.^2)./ORIGS0;

Ao = lin2rgb(ORIGDOLP);
fn = ['./'  num2str(z) '_Dolp_orig_.png'];
name = ['Dolp_orig_',num2str(z)];
imwrite(Ao, fn, 'png');
% ORIGDOLP = lin2rgb(ORIGDOLP);
% imshow(ORIGDOLP);
ORIG1=R1;
ORIG3=R4;
ORIG5=R6;




CFA = create_CFA(ORIG1,ORIG3,ORIG5);

[S0,S1,S2,a1,c1,e1]= demos_LSLCD(CFA,fpkg);

S0=real(S0);
I_disp = S0 / 2;
imagesc(I_disp);caxis([0 1]);
colorbar;

% fn = ['./'  num2str(z) '_S0_42_.png'];
% imwrite(I_disp, fn, 'png');

% imshow(S0);
S1=real(S1);
S2=real(S2);
a1=real(a1);
c1=real(c1);
e1=real(e1);
% [S0,S1,S2,DOLP]=S0Calcalate(a1,c1,e1);
DOLP=sqrt(S1.^2+S2.^2)./S0;
I_disp = ORIGDOLP;
I_disp = imagesc(I_disp);caxis([0 1]);
colorbar;
Co = lin2rgb(DOLP);
fn = ['./'  num2str(z) '_Dolp_4_.png'];
imwrite(Co, fn, 'png');

k1= cpsnr(ORIG1,ORIG3,ORIG5,a1,c1,e1,0);
k2= cpsnr1(ORIGS0,S0,0);
k3= cpsnr1(ORIGS1,S1,0);
k4= cpsnr1(ORIGS2,S2,0);
k5= cpsnr1(ORIGDOLP,DOLP,0);
k6 = (ssim(ORIG1,a1)+ssim(ORIG3,c1)+ssim(ORIG5,e1))/3;
k7 = ssim(ORIGS0,S0);
k8 = ssim(ORIGDOLP,DOLP);
    if k5<0
        print('error')
        break
    end
    N_sum(:,z) = [k1;k2;k3;k4;k5;k6;k7;k8];
%     N_sum(SNR-10,z) = abs(k1)+abs(k2)+abs(k5);
    end
N = sum(N_sum,2)/40;
% end
% N_sum_four = sum(N_sum,2)/size(N_sum,2);
% j=11:1:40;
% plot(j,N_sum_four);
% xlabel('SNR');ylabel('PSNR');