function [S_H, S_P] = solveLS(A,b) 
% �������A��ϵ������ 
% �������b��Ax=b�ĳ�����������b 
% S_H��������Է�����Ļ�����ϵ 
% S_P����������Է�������ؽ� 
if size(A,1) ~= length(b) %size(A,1)���������� 
    error('�������ݴ������������룡'); 
    return; 
else
    B = [A,b]; %������� 
    rank_A = rank(A); %��ϵ��������� 
    rank_B = rank(B); %������������ 
    if rank_A ~= rank_B %�޽���� 
        disp('���Է������޽⣡'); 
        S_H = []; S_P = []; 
    else if rank_B == size(A,2) %������������ = δ֪������ 
            %size(A,2)�������������൱��length(A) 
            disp('���Է�������Ψһ�⣡'); 
            S_P = A\b; %��Ψһ�� 
            S_H = []; 
        else
            disp('���Է�����������⣡'); 
            S_H = null(A,'r');%�����η�����Ļ�����ϵ 
            S_P = A\b; %�����η�������ؽ� 
        end
    end
end





