clear;clc;
ACC=0;NMI=0;PUR=0;
acc = [];nmi = [];pur = [];
type=1;
cir_num=5;
for cir=1:cir_num
addpath('mtv_data');
addpath('Utilities');

% ORL dataset
data=cell(1,3);
load('COIL20_Gaussian_Random_0.02_view1.mat');data{1}=X1;
load('COIL20_Gaussian_Random_0.02_view2.mat');data{2}=X2;
load('COIL20_Gaussian_Random_0.02_view3.mat');data{3}=X3;

X=data;
viewnum = size(X,2);
y = gt;
classNum = length(unique(y));
distX1 = cell(1,viewnum);
idx = cell(1,viewnum);
S = cell(1,viewnum);
m = zeros(1,viewnum);
n = zeros(1,viewnum);
for v = 1:viewnum
    [m(v),n(v)]=size(X{v});
end
for v = 1:viewnum
    S(v) = {constructW_PKN(X{v},2)};
end

% Parameter Setting
beta =0.1;
iter =6;
mu =0.1;
delta=0.1;
k = 6;
% Initialization 
Y1 = cell(1,viewnum);
Y2 = cell(1,viewnum);
Z = cell(1,viewnum);
E = cell(1,viewnum);
J = cell(1,viewnum);
N = size(X{1},2);
V = length(X);
sX = [N, N, V];
for v=1:viewnum
    Y1{v} =zeros(m(v),n(v));
    Y2{v} =zeros(n(v),n(v));
    Z{v} = (X{v}'*X{v}+0.001*eye(n(v)))\(X{v}'*X{v});
    E{v} = zeros(m(v),n(v));
    S{v} = zeros(n(v),n(v));
    H{v} = zeros(n(v),n(v));
end

for i = 1:iter
     for v = 1:viewnum
           J{v} = updateJ(Z{v},Y2{v},mu);            
           distX = L2_distance_1(X{v}*Z{v},X{v}*Z{v});
           [distX1, idx] = sort(distX,2);
           [gamma] = cal_gamma(X{v},distX1,beta,k);
           S{v} = updateS(X{v}, distX1,idx,k,gamma,beta,H{v},S{v});
           S{v} = (S{v}+S{v}')/2;
           L{v} = diag(sum(S{v}))-S{v};
           Z{v} = updateZ(X{v},beta,mu,Y2{v},Y1{v},L{v},E{v},J{v});
           E{v} = updateE(X{v},1./sqrt(max(m(v),n(v))),Z{v},Y1{v},mu,type);
           Y1{v} = Y1{v}+mu*(X{v}*Z{v}+E{v}-X{v});
           Y2{v} = Y2{v}+mu*(J{v}-Z{v});
     end
         S_tensor = cat(3, S{:,:});
         s = S_tensor(:);
         H = updateH(s,sX,delta,viewnum);
         
         D = zeros(size(H{1}));  
        for j=1:viewnum
           D = D + (abs(H{j})+(abs(H{j}))')/2;
        end
                A = D/viewnum; 
                [actual_ids] = spectral_clustering(A, classNum);
                result = ClusteringMeasure(actual_ids ,y);
                acc = [result(1) acc];
                nmi = [result(2) nmi];
                pur = [result(3) pur];
end
ACC = ACC + max(acc);
NMI = NMI + max(nmi);
PUR = PUR + max(pur);
acc = []; nmi = []; pur = [];
end
disp(['ACC: ' num2str(ACC/cir_num)]);
disp(['NMI: ' num2str(NMI/cir_num)]);
disp(['PUR: ' num2str(PUR/cir_num)]);
