% l1eq_example.m  
%  
% Test out l1eq code (l1 minimization with equality constraints).  
%  
% Written by: Justin Romberg, Caltech  
% Email: jrom@acm.caltech.edu  
% Created: October 2005  
%  
  
% put key subdirectories in path if not already there  
% path(path, './Optimization');  
% path(path, './Data');  
  
% To reproduce the example in the documentation, uncomment the   
% two lines below  
%load RandomStates  
%rand('state', rand_state);  
%randn('state', randn_state);  
%function tim = BP_no_noise (K)
%clear all;close all;clc; 
% signal length  
N = 16467; %2450;  
% number of spikes in the signal  
T = 330;  
% number of observations to make  
K = 1320;  
  
% random +/- 1 signal  
%x = zeros(N,1);  
%q = randperm(N);  
%x(q(1:T)) = sign(randn(T,1));  
%下面只产生0和1的信号
%
x=zeros(N,1); %生成0矩阵
p = randperm(N); %随机打乱序列
x1=unidrnd(1,T,1); %产生从1到N所指定的最大数数之间的离散均匀随机整数
x(p(1:T))=x1;

%[x]=xlsread('C:\MATLAB7\work\X矩阵2500.xlsx')
yi=sum(~~x(:));
% measurement matrix  
disp('Creating measurment matrix...');  
%p=0.01;
%A = rand(K,N);  
%A =A<p
%[A]=xlsread('C:\MATLAB7\work\X矩阵1_1.xlsx')
disp('Done.');  
  A=randsrc(K,N,[[0 1];[0.96 0.04]])  ; 
  % 生成K行 2k列矩阵 0出现概率0.98,1出现概率0.02
% observations  



y = A*x;    
% initial guess = min energy  
%x0 = A'*y;  

% take epsilon a little bigger than sigma*sqrt(K)  

% solve the LP  
 tic  %计时开始
 xp = l1eq_pd(x, A, [], y, 1e-3);  
 tim = toc;  %计时结束
  
% large scale  
%Afun = @(z) A*z;  
%Atfun = @(z) A'*z;  
%tic  
 

%xp = l1eq_pd(x0, Afun, Atfun, y, 1e-3, 30, 1e-8, 200);toc 

if(~isempty(xp))
    % 对恢复信号进行取整操作
    xp = round(xp);
    
    N = length(xp);
    for k = 1:N
        %因为只有0和1 所以对非0的信号全取1
        %if xp(k) == 0
        %else xp(k) = 1
        % 对恢复信号中大于 1 的全取为 1
        if xp(k) >= 1
            xp(k) = 1;
        end
    end

end

%对未恢复信号标注
%figure(1);
xn = xp - x;
xn = abs(xn);
%plot(xn,'r'); %绘出未恢复的x信号
%hold on;
%plot (x,'k.');
%hold off;
%legend('Recovery failed signal','Original signal');

%figure(2);
%plot(xp,'r');%绘出x的恢复信号    
%hold on;    
%plot(x,'k');%绘出原信号x    
%hold off;    
%legend('Recovery','Original')    
%fprintf('\n恢复残差：');    
%norm(xp-x)%恢复残差 

%对原始信号和恢复信号进行比较
figure;    
plot(xp,'k');%绘出x的恢复信号    
hold on;    
plot(x,'r');%绘出原信号x    
hold off;    
legend('Recovery','Original')    
%fprintf('\n恢复残差：');    
%norm(xp-x)%恢复残差  

figure;
plot(y,'b');
legend('y(i)');


fprintf('\n未恢复x的个数为%d',norm(xn,1));
fprintf('\n运行的时间为%d',tim);
s = norm(xn,1);


%f=xp-x;
%s=sum(~~f(:))
Cur = (N-s)/N;
fprintf('\n恢复正确率：');
fprintf('%2.2f%%',((N-s)/N)*100)