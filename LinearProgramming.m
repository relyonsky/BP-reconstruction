clear;close all;clc; 
% signal length  
N = 2450;  
% number of sp ikes in the signal  
T = 147;  
% number of observations to make  
K = 441;  



%生成x向量
x=zeros(N,1); %生成0矩阵
p = randperm(N); %随机打乱序列
x1=unidrnd(1,T,1); %产生从1到N所指定的最大数数之间的离散均匀随机整数
x(p(1:T))=x1;

%[x]=xlsread('C:\MATLAB7\work\X矩阵2500.xlsx')
yi=sum(~~x(:));
% measurement matrix  
disp('Creating measurment matrix...');  
  A=randsrc(K,2450,[[0 1];[0.99 0.01]])  ; 
  % 生成K行 2k列矩阵 0出现概率0.99,1出现概率0.01
% observations  


y = A*x;


figure(1);
plot(y,'r*');
legend('y(i)');%画出观测信号y

xp = 2 * ones(N,1);%将恢复信号元素均设置成2


%步骤一 
for i=1:K
    if y(i) == 0
        for j = 1:N
            if A(i,j) ~= 0
                xp(j) = 0;
            end
        end
    end
end
%将所有y(i)=0情况考虑进去 得出n个确定的x为0


p = sum(xp()==0);
fprintf('\n第一步x(j)=0的个数为%d',p);

m = sum(y()==0); n = sum(y()==1); k = sum(y()==2);
fprintf('\n');
fprintf('y中0，1，2的个数为%d,%d,%d',m,n,k);


%步骤二
for j = 1:N
    if xp(j) == 0     %消除所有x为0的情况
        index(j) = 1;
    else
        index(j) = 0;
    end
end
xn = xp(~index,:);%xn是去除所有x(i)=0的列之后的x
Q(:,:) = A(:,~index);%Q是去除所有已经确定x(i)=0的列后的观测矩阵
 
%抽离y=1的情况
 n = 1;
 for i=1:K
     if y(i) == 1
         yn(n,:) = y(i);
         M(n,:) = Q(i,:); % 现在考虑的是M * xn = T
         n=n+1;           % M大小(n-1)*(2450-p)
     end
 end

xr= M\yn;
t1=length(xr);
for i=1:t1
    T(i)=sum(M(:,i));
end
% for i=1:t1
%     if xr(i)~=0
%         xr(i) = 1;
%     end
% end
% t2=sum(xr);
 
 
%求解抽离后的方程
%for i=1:K
%    if xp(i) ~= 0 || A(n-1,i) == 1;
%        xn(i) = 0;
%    end
%end
% format rat;
% [S_H, S_P]= solveLS(M,T); 
%这种方法解出来的情况是没有考虑步骤一的0解
%难点是如何在使得前面解出来x(i)=0的情况，带入步骤二进行再次求解

        
                              
%p = sum(y() == 1);

                
            
            
