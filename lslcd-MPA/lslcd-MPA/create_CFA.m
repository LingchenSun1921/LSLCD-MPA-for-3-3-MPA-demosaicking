function [CFA] = create_CFA(a,b,c)
%%%%%%%[0 120 60;120 60 0;60 0 120]
S = size(a); N1 = S(1); N2 = S(2);


a0=[1 0 0;0 0 1;0 1 0];

a60=[0 0 1;0 1 0;1 0 0];

a120=[0 1 0;1 0 0;0 0 1];

n1=size(a0,1);
n2=size(a0,2);


%Write four cosets of RGB to CFA
a0=repmat(a0,ceil(N1/n1),ceil(N2/n2));
a0(N1+1:size(a0,1),:)=[];
a0(:,N2+1:size(a0,2))=[];

a60=repmat(a60,ceil(N1/n1),ceil(N2/n2));
a60(N1+1:size(a60,1),:)=[];
a60(:,N2+1:size(a60,2))=[];

a120=repmat(a120,ceil(N1/n1),ceil(N2/n2));
a120(N1+1:size(a120,1),:)=[];
a120(:,N2+1:size(a120,2))=[];


CFA=a0.*a+a60.*b+a120.*c;