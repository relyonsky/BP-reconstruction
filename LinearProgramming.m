clear;close all;clc; 
% signal length  
N = 2450;  
% number of sp ikes in the signal  
T = 147;  
% number of observations to make  
K = 441;  



%����x����
x=zeros(N,1); %����0����
p = randperm(N); %�����������
x1=unidrnd(1,T,1); %������1��N��ָ�����������֮�����ɢ�����������
x(p(1:T))=x1;

%[x]=xlsread('C:\MATLAB7\work\X����2500.xlsx')
yi=sum(~~x(:));
% measurement matrix  
disp('Creating measurment matrix...');  
  A=randsrc(K,2450,[[0 1];[0.99 0.01]])  ; 
  % ����K�� 2k�о��� 0���ָ���0.99,1���ָ���0.01
% observations  


y = A*x;


figure(1);
plot(y,'r*');
legend('y(i)');%�����۲��ź�y

xp = 2 * ones(N,1);%���ָ��ź�Ԫ�ؾ����ó�2


%����һ 
for i=1:K
    if y(i) == 0
        for j = 1:N
            if A(i,j) ~= 0
                xp(j) = 0;
            end
        end
    end
end
%������y(i)=0������ǽ�ȥ �ó�n��ȷ����xΪ0


p = sum(xp()==0);
fprintf('\n��һ��x(j)=0�ĸ���Ϊ%d',p);

m = sum(y()==0); n = sum(y()==1); k = sum(y()==2);
fprintf('\n');
fprintf('y��0��1��2�ĸ���Ϊ%d,%d,%d',m,n,k);


%�����
for j = 1:N
    if xp(j) == 0     %��������xΪ0�����
        index(j) = 1;
    else
        index(j) = 0;
    end
end
xn = xp(~index,:);%xn��ȥ������x(i)=0����֮���x
Q(:,:) = A(:,~index);%Q��ȥ�������Ѿ�ȷ��x(i)=0���к�Ĺ۲����
 
%����y=1�����
 n = 1;
 for i=1:K
     if y(i) == 1
         yn(n,:) = y(i);
         M(n,:) = Q(i,:); % ���ڿ��ǵ���M * xn = T
         n=n+1;           % M��С(n-1)*(2450-p)
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
 
 
%�������ķ���
%for i=1:K
%    if xp(i) ~= 0 || A(n-1,i) == 1;
%        xn(i) = 0;
%    end
%end
% format rat;
% [S_H, S_P]= solveLS(M,T); 
%���ַ���������������û�п��ǲ���һ��0��
%�ѵ��������ʹ��ǰ������x(i)=0����������벽��������ٴ����

        
                              
%p = sum(y() == 1);

                
            
            
