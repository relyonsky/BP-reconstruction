clear; 
T = 147;
for i=4:4
K = i * T;
s = 0;
for j=1:50        %循环20次取均值
    s0 = BP_no_noise(K);    
    s = s+s0;
end
s = s/50;
%fprintf('i = %d , s = %d \n',i,s)
fid = fopen('tim_p01.xls','a');
fprintf(fid,'%d \t',s);
fprintf(fid,'\r\n');
fclose(fid);
end
 