clear; 
T = 147;
for mi = 0.04:0.02:0.1
for i=5:8
K = i * T;
s = 0;
for j=1:50        %循环30次取均值
    s0 = BPDN(K,mi);
    s = s+s0;
end
s = s/50;
%fprintf('i = %d , s = %d \n',i,s)
fid = fopen('tim.xls','a');
fprintf(fid,'%d \t',s);
fprintf(fid,'\r\n');
fclose(fid);
end
end