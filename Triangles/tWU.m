% rule UFE > prct. of RAND 

function t_WU = tWU(v)

global stat EE zz Wei

det1 = 0;
FP1 = 0;

[N] = find(stat(:,8) > v);
       
    c = isempty(N);
 if c == 0
     
for n = 1 : size(N,1) 
       
             if EE(stat(N(n),1), stat(N(n),2)) == 1
                 
                det1 = det1+1;
             else
                 FP1 = FP1+1;
             end
             
end
 
         det1 = - det1/zz;   %% detection (0-1) 
          FP1 = FP1/n;    %% False pos. (0-1)
          
          t_WU = det1+Wei*FP1;  %% minimise min(-det + FP)
 else
          t_WU = 1+rand;
end

if t_WU == 1 
     t_WU = 1+rand;
 else
 end 
 
end