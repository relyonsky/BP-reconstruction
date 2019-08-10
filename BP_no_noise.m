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
%����ֻ����0��1���ź�
%
x=zeros(N,1); %����0����
p = randperm(N); %�����������
x1=unidrnd(1,T,1); %������1��N��ָ�����������֮�����ɢ�����������
x(p(1:T))=x1;

%[x]=xlsread('C:\MATLAB7\work\X����2500.xlsx')
yi=sum(~~x(:));
% measurement matrix  
disp('Creating measurment matrix...');  
%p=0.01;
%A = rand(K,N);  
%A =A<p
%[A]=xlsread('C:\MATLAB7\work\X����1_1.xlsx')
disp('Done.');  
  A=randsrc(K,N,[[0 1];[0.96 0.04]])  ; 
  % ����K�� 2k�о��� 0���ָ���0.98,1���ָ���0.02
% observations  



y = A*x;    
% initial guess = min energy  
%x0 = A'*y;  

% take epsilon a little bigger than sigma*sqrt(K)  

% solve the LP  
 tic  %��ʱ��ʼ
 xp = l1eq_pd(x, A, [], y, 1e-3);  
 tim = toc;  %��ʱ����
  
% large scale  
%Afun = @(z) A*z;  
%Atfun = @(z) A'*z;  
%tic  
 

%xp = l1eq_pd(x0, Afun, Atfun, y, 1e-3, 30, 1e-8, 200);toc 

if(~isempty(xp))
    % �Իָ��źŽ���ȡ������
    xp = round(xp);
    
    N = length(xp);
    for k = 1:N
        %��Ϊֻ��0��1 ���ԶԷ�0���ź�ȫȡ1
        %if xp(k) == 0
        %else xp(k) = 1
        % �Իָ��ź��д��� 1 ��ȫȡΪ 1
        if xp(k) >= 1
            xp(k) = 1;
        end
    end

end

%��δ�ָ��źű�ע
%figure(1);
xn = xp - x;
xn = abs(xn);
%plot(xn,'r'); %���δ�ָ���x�ź�
%hold on;
%plot (x,'k.');
%hold off;
%legend('Recovery failed signal','Original signal');

%figure(2);
%plot(xp,'r');%���x�Ļָ��ź�    
%hold on;    
%plot(x,'k');%���ԭ�ź�x    
%hold off;    
%legend('Recovery','Original')    
%fprintf('\n�ָ��в');    
%norm(xp-x)%�ָ��в� 

%��ԭʼ�źźͻָ��źŽ��бȽ�
figure;    
plot(xp,'k');%���x�Ļָ��ź�    
hold on;    
plot(x,'r');%���ԭ�ź�x    
hold off;    
legend('Recovery','Original')    
%fprintf('\n�ָ��в');    
%norm(xp-x)%�ָ��в�  

figure;
plot(y,'b');
legend('y(i)');


fprintf('\nδ�ָ�x�ĸ���Ϊ%d',norm(xn,1));
fprintf('\n���е�ʱ��Ϊ%d',tim);
s = norm(xn,1);


%f=xp-x;
%s=sum(~~f(:))
Cur = (N-s)/N;
fprintf('\n�ָ���ȷ�ʣ�');
fprintf('%2.2f%%',((N-s)/N)*100)