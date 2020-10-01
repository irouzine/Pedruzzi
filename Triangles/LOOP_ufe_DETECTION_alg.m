clc
close all
clear all %#ok<CLALL>
clearvars 

global stat zz Wei EE

tic



% Parameters
                                              Wei = 1;
                                              L = 50;  
                                              N = 10000;  N1 = N;
                                              
  s0 = .1; 
E_st = .75;
  f0 = .45*L;                                                
   U = 0.07; 
  td = 20+1;  
 
                                                                           
muL = U;  
mu=muL/L; 

c_1 = N * (mu/s0)^2;
SUlogNs = s0/(mu*L)*log(N*s0);

EE = zeros(L,L);
det_mat = zeros(L,L,td);


%create a random matrix:

a = 1;
b = L;

r = ceil(-a + (b-a)*rand(1,L));


for i = 2:2:L   
 
%     if r(i) ~= 0
    
% EE(i,r(i)) = 1; 
% EE(r(i),i) = 1;
EE(i,i-1) = 1; 
EE(i-1,i) = 1;

%     else
%     end
    
end

Ee = tril(EE,-1);   
                                                                
 [DAt, ww] = epi('binary',  s0,  0, L, N,  td,  2,   U,  f0, 1 , E_st, Ee);   

 
 
 
                                   tmin = 1;
                                  tmax = td-1;
 
                                                                           for t = tmin:tmax % time of analysis

  A = 6; 

I = find(Ee); 
    
  tic                                                                          
     
DATA = DAt(:,:,t); 
 
zz = length(I);
W = ww(:,t); 

% add detection code


    
    Ko = DATA; N = N1; 
    
    
i_div = (1:L); ip = nchoosek(i_div,2); np = size(ip,1);
fhap = zeros(np,4); UFE = zeros(np,1); WU = zeros(np,1);

for i=1:np
f00 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(N,1)*[0 0],2)); f01 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(N,1)*[0 1],2));
f10 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(N,1)*[1 0],2)); f11 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(N,1)*[1 1],2));
                 
fhap(i,:) = [f00 f01 f10 f11]; UFE(i) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); 
end  
   

UFE(UFE<0) = 0; 
UFE(UFE>1) = 0; 

stat = horzcat(ip,fhap,UFE);
     
            %%
 
m = 1; mm = m; tt = 1; 
 zz = length(I); n = 7;  
 
pick = ceil(rand(1,200)*np);
Pi = sort(pick);
stat_r = stat(Pi, :); 
 
      [I,J] = find(Ee); 
      
  % --------------------------------------------      

  % Threeshold
  
%    t00 = prctile(stat_r(:,3),80);
%    t01 = prctile(stat_r(:,4),80);
%    t10 = prctile(stat_r(:,5),80);     
%    t11 = prctile(stat_r(:,6),80);

   tUF = prctile(stat_r(:,7),80); 
   
  % ----------------------------------------
  
  v = tUF; [tUF_opt] = fminsearch(@tUFE,v);
  
[N2] = find(stat(:,7) > tUF_opt);

 det = 0;FP = 0;  c4 = isempty(N2);
    if c4 == 0

for n2 = 1 : size(N2,1) 
    
    % reconstruct detected matrix for comparison with Ee (EE)
    
    
    
             if EE(stat(N2(n2),1), stat(N2(n2),2)) == 1
                 
                 det_mat(stat(N2(n2),1), stat(N2(n2),2),t) = -1; % detected
                 
                 
                det = det+1;
             else
                 FP = FP+1;
                 
                 det_mat(stat(N2(n2),1), stat(N2(n2),2),t) = 1; % false positive
                 
             end    
end

         det = 100*det/zz; FP = 100*FP/n2;    
         out22 = [det;FP];
          else
      out22 = [det;FP];
    end
       
    
  % -------------------------
 

  % -------------------------
    
 
% PLOT

Y = horzcat(out22);

   x_str = {'UFE'};
                 
          
     OUT_det_UF(1,t) = Y(1,1); %#ok<SAGROW>
     
     OUT_fp_UF(1,t) = Y(2,1);  %#ok<SAGROW>


                                                                           end


%%

     Det_UF = mean(OUT_det_UF,3); 
     
     Fp_UF = mean(OUT_fp_UF,3); 
      
% ----------------------------
figure(1)
subplot(1,2,1), plot(Det_UF','LineWidth',2), title('Detection UFE'),  hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax])
subplot(1,2,2), plot(Fp_UF','LineWidth',2),  title('False Pos. UFE'), hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax])

    legend off




%----------------------------
gm = figure(1);

fig_n2 = 'FIGURE';
set(gm, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

set(gm,'Units','Inches');pos = get(gm,'Position');
set(gm,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gm,fig_n2,'-dpdf','-r0')



%%

figure(2)

Z = peaks;
surf(Z); 
axis tight manual 
set(gca,'nextplot','replacechildren'); 


v = VideoWriter('DET.avi', 'Uncompressed AVI');
open(v);

E_k = -EE;

for k = 1:20

    
% tot = det_mat(:,:,k) + E_k ;

% tot = tril(E_k) + triu(E_k).*det_mat(:,:,k);
% FP = triu(-E_k)+det_mat(:,:,k);

% tot(tot == )

  a =  surf(det_mat(:,:,k),'LineStyle','none'); colormap(redblue), view([0 90]), axis square, 
  ylabel('KNOWN Epis. pairs'), xlabel('DETECTED Epis. pairs'), colorbar, caxis ([-1 1]), 
  
   frame = getframe(gcf);
   writeVideo(v,frame);
   
end

close(v);

%%
implay('DET.avi');



