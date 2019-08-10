%ֱ������ʱ�� fuction�� K�ĸ�ֵ����ע��ȡ����
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
%����ֻ����0��1���ź�
%[x]=xlsread('C:\MATLAB7\work\X����.xlsx')
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
%[A]=xlsread('C:\MATLAB7\work\A����.xlsx')
A = randsrc(K,N,[[0 1];[0.94 0.06]]);%������Ըı����A���ܶ��Դﵽ�ϸߵ�׼ȷ��
disp('Done.');  
      
%e = round (1 * wgn(K,1,1)); %������˹����������ȡ��
 e = zeros(K,1);
for i=1:length(e)
%     if rem(i,15) == 0  %0.3�ļ���
%         e(i) = round (-2+ 4*rand);
%     elseif rem(i,15) == 1 
%         e(i) = round (-2+ 4*rand);   
%     elseif rem(i,15) == 7
%         e(i) = round (-2+ 4*rand);
%     else
%         e(i) = 0; %����-2��2���������
    if rem(i,10) ~= 0  %��4/5��λ��ȡ0
        e(i) = 0;
    else
        e(i) = round (-3+ 6*rand); %����-2��2���������
     end
end
%   

% ���������ź�e(i)
% e = zeros(K,1);
% for i=1:K 
% e(i) = round (-2+ 4*rand);
% end
% randomIndex = 1+floor(rand(1,floor(K*0.20))*K); %������0.8���ʵ���ֵΪ0
% 
% e(randomIndex) = 0;
   
%e = random('poisson',0.08,K,1);  

    %ΪɶҪ��Ϊ��ֵ��
    %for i=1:length(e)
    %if(rand > 0.5)
    %    e(i) = -e(i);
    %end
    %end
    %�ϱߵĲ��ɷֲ�������ȱ�ݣ�������ķ�����������

    %��������������ó���ʹ�����������������Լ�0�ĸ���ȷ��
    %e=round(UNIFRND(-2,2,256,1));
    %randomIndex=1+floor(rand(1,floor(256*0.9))*256);
    %e(randomIndex)=0;
    %������������⣬�ظ��ĵ�̫�࣬����0�ĸ�������


    %����������������������������570��Ϊ�˶��Ҳ��ֵ㣬�������ظ�������ĸ���
    %e=[]
    %for i =1:514
     %   e(i,1) = round(-2+ 4*rand); -2��2ȡ��
    %end
    %B = randint(1,900,[1,514]); 1��900��[1��514��Χ�������]
    %for i=1 :900
    %e(B(i))=0; 514��Ԫ�������0���ظ�900��
    %end
    %e
    %��������е���Զ·����������������������λ���㣬���±ߵ��·���


    %e=zeros(K,1);
    %p = randperm(K);
    %z=unidrnd(3,7,1);%����7*1��������ȡֵΪ1,2,3
    %z1=randi(7,1,6); % z1=randint(1,6,[1 7])
    %for i=1:6
    %    z(z1(i))=-z(z1(i));
    %end
    %z
    %e(p(1:7))=z
    %�����̶���������

y0 = A*x;

% initial guess = min energy  
 

% noise_amount=length(nonzeros(e));
% noise_total=sum(abs(e));


%t(i)Ϊ������ԭ�źű�ֵ
t = zeros(1,length(e));%�otԤ�����ڴ�
for i=1:length(e)    
    if y0(i) ~= 0
        t(i) = e(i)/y0(i);
        if abs(t(i)) > 0.6
            e(i) = fix(0.5 * y0(i)); %ɸ�����������ԭ�źŹ����
            t(i) = e(i)/y0(i); %�����޸Ĺ����t(i)����ֵ 
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

y = A*x + e;  %����eֵ��y
x0 = A'* y; 

% take epsilon a little bigger than sigma*sqrt(K)  
sigma = 0.005;  
epsilon =  sigma*sqrt(K)*sqrt(1 + 2*sqrt(2)/sqrt(K)); %����yi pu sei le
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
%      %�Իָ��źŽ���ȡ������
%     xp = round(xp);
%     % �Իָ��ź��д��� 1 ��ȫȡΪ 1
%    
%    for k = 1:length(xp)
%         if xp(k) > 1
%             xp(k) = 1;
%         end
%    end
% 
% end

%figure(1);    
%plot(xp,'r');%���x�Ļָ��ź�    
%hold on;    
%plot(x,'k');%���ԭ�ź�x    
%hold off;    
%legend('Recovery','Original')    
%fprintf('\n�ָ��в');    
%norm(xp-x)%�ָ��в�  

% figure(2);
% plot(x,'k');%���ԭ�ź�x
% xn = xp - x;
% xn = abs(xn); %δ�ָܻ����ź�
% plot(t,'k.-'); %���������
% legend('signal noise ratio');

% figure(3);%��������
% plot(e);

figure(4);%����ԭ�۲��ź�
plot(y0);
% 
% figure(5);%�����������Ĺ۲��ź�
% plot(y);

%f=xp-x;
s=sum(~~f(:));
Cur = (N-s)/N;
fprintf('\n�ָ���ȷ�ʣ�');
fprintf('%2.2f%%\n',((N-s)/N)*100);
%end

%fprintf('time is %d',tim);
