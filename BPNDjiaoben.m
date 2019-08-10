clear; 
T = 147;
for i=7:8
K = i * T;
s = 0;
for j=1:30        %循环20次取均值
    s0 = BPDN(K);    
    s = s+s0;
end
s = s/30;
%fprintf('i = %d , s = %d \n',i,s)
fid = fopen('tim030.xls','a');
fprintf(fid,'%d \t',s);
fprintf(fid,'\r\n');
fclose(fid);
end
 