%直接运行时将 fuction和 K的赋值命令注释取消掉
%function tim = BPDN (K)

clear all;
clc;
 
% signal length  
% N = 2450;  
% number of spikes in the signal  
% T = 147;  
% number of observations to make  
% K = 588;  

N = 16467;
T = 330;
K = 1318;


  
% random +/- 1 signal  
%x = zeros(N,1);  
%q = randperm(N);  
%x(q(1:T)) = sign(randn(T,1));  
%下面只产生0和1的信号
%[x]=xlsread('C:\MATLAB7\work\X矩阵.xlsx')
x=zeros(N,1);
p = randperm(N);
x1=unidrnd(1,T,1);
x(p(1:T))=x1;
%yi=sum(~~x(:));

% measurement matrix  
disp('Creating measurment matrix...');  
%p=0.04;
%A = rand(K,N);  
%A =A<p
%[A]=xlsread('C:\MATLAB7\work\A矩阵.xlsx')
A = randsrc(K,N,[[0 1];[0.94 0.06]]);%这里可以改变矩阵A的密度以达到较高的准确率
disp('Done.');  
      
%e = round (1 * wgn(K,1,1)); %产生高斯白噪声并且取整
 e = zeros(K,1);
for i=1:length(e)
%     if rem(i,15) == 0  %0.3的几率
%         e(i) = round (-2+ 4*rand);
%     elseif rem(i,15) == 1 
%         e(i) = round (-2+ 4*rand);   
%     elseif rem(i,15) == 7
%         e(i) = round (-2+ 4*rand);
%     else
%         e(i) = 0; %产生-2到2的随机噪声
    if rem(i,10) ~= 0  %将4/5的位置取0
        e(i) = 0;
    else
        e(i) = round (-3+ 6*rand); %产生-2到2的随机噪声
     end
end
%   

% 生成噪声信号e(i)
% e = zeros(K,1);
% for i=1:K 
% e(i) = round (-2+ 4*rand);
% end
% randomIndex = 1+floor(rand(1,floor(K*0.20))*K); %噪声中0.8概率的数值为0
% 
% e(randomIndex) = 0;
   
%e = random('poisson',0.08,K,1);  

    %为啥要变为负值？
    %for i=1:length(e)
    %if(rand > 0.5)
    %    e(i) = -e(i);
    %end
    %end
    %上边的泊松分布噪声有缺陷，用另外的方法产生噪声

    %随机噪声，并且用程序使其向量的噪声个数以及0的个数确定
    %e=round(UNIFRND(-2,2,256,1));
    %randomIndex=1+floor(rand(1,floor(256*0.9))*256);
    %e(randomIndex)=0;
    %这个方法有问题，重复的点太多，导致0的个数不够


    %现在用这个方法来产生随机噪声，570是为了多找部分点，来避免重复点带来的干扰
    %e=[]
    %for i =1:514
     %   e(i,1) = round(-2+ 4*rand); -2到2取整
    %end
    %B = randint(1,900,[1,514]); 1到900个[1到514范围内随机数]
    %for i=1 :900
    %e(B(i))=0; 514个元素随机赋0，重复900次
    %end
    %e
    %这个方法有点绕远路，生成了噪声还给它部分位置零，用下边的新方法


    %e=zeros(K,1);
    %p = randperm(K);
    %z=unidrnd(3,7,1);%生成7*1的列向量取值为1,2,3
    %z1=randi(7,1,6); % z1=randint(1,6,[1 7])
    %for i=1:6
    %    z(z1(i))=-z(z1(i));
    %end
    %z
    %e(p(1:7))=z
    %产生固定个数噪声

y0 = A*x;

% initial guess = min energy  
 

% noise_amount=length(nonzeros(e));
% noise_total=sum(abs(e));


%t(i)为噪声与原信号比值
t = zeros(1,length(e));%給t预分配内存
for i=1:length(e)    
    if y0(i) ~= 0
        t(i) = e(i)/y0(i);
        if abs(t(i)) > 0.6
            e(i) = fix(0.5 * y0(i)); %筛除噪声相对于原信号过大的
            t(i) = e(i)/y0(i); %更新修改过后的t(i)的数值 
            %e(i) = 0; y(i) = 0;
        end
    else
        t(i) = 0;
    end
    if y0(i)-e(i)<0
        e(i) = y0(i);
    end
end
%t = t';

y = A*x + e;  %更新e值给y
x0 = A'* y; 

% take epsilon a little bigger than sigma*sqrt(K)  
sigma = 0.005;  
epsilon =  sigma*sqrt(K)*sqrt(1 + 2*sqrt(2)/sqrt(K)); %参数yi pu sei le
% solve the LP  
% tic  
% xp = l1eq_pd(x0, A, [], y, 1e-3);  
% toc  
 
tic  
% xp = l1qc_logbarrier(x0, A, [], y, epsilon, 1e-3);  
tim = toc;

% large scale  
%Afun = @(z) A*z;  
%Atfun = @(z) A'*z;  
%tic  
%xp = l1qc_logbarrier(x0, Afun, Atfun, y, epsilon, 1e-3, 50, 1e-8, 500);   

%xp = l1eq_pd(x0, Afun, Atfun, y, 1e-3, 30, 1e-8, 200);toc 

% if(~isempty(xp))
%      %对恢复信号进行取整操作
%     xp = round(xp);
%     % 对恢复信号中大于 1 的全取为 1
%    
%    for k = 1:length(xp)
%         if xp(k) > 1
%             xp(k) = 1;
%         end
%    end
% 
% end

%figure(1);    
%plot(xp,'r');%绘出x的恢复信号    
%hold on;    
%plot(x,'k');%绘出原信号x    
%hold off;    
%legend('Recovery','Original')    
%fprintf('\n恢复残差：');    
%norm(xp-x)%恢复残差  

% figure(2);
% plot(x,'k');%绘出原信号x
% xn = xp - x;
% xn = abs(xn); %未能恢复的信号
% plot(t,'k.-'); %画出信噪比
% legend('signal noise ratio');

% figure(3);%画出噪声
% plot(e);

figure(4);%画出原观测信号
plot(y0);
% 
% figure(5);%画出含噪声的观测信号
% plot(y);

%f=xp-x;
s=sum(~~f(:));
Cur = (N-s)/N;
fprintf('\n恢复正确率：');
fprintf('%2.2f%%\n',((N-s)/N)*100);
%end

%fprintf('time is %d',tim);
