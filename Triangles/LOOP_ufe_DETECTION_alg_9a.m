clc
close all
clear all %#ok<CLALL>
clearvars 


tic

% Parameters
ii = 0;

                                             % L = 10;
                                             for L = [10]
                                             
%                                                  if ii == 0 
                                                     th = .5;
%                                                  elseif ii == 2
% %                                                      th = .6;
%                                                  else
%                                                  end
                                                 
                                           runn = 200;  
                                              N = 500;  N1 = N; 
                                              
                                              
                                              for tp = 2%1:4
                                              
  s0 = .1; 
E_st = .75;
  f0 = .45*L;                                                
   U = 0.07; 
  td = 19+1;  
 
        r = 0; 
      M = 0;  
                                                                      
muL = U;  
mu=muL/L; 

c_1 = N * (mu/s0)^2;
SUlogNs = s0/(mu*L)*log(N*s0);

EE = zeros(L,L);
det_mat = zeros(L,L,td);


%create a random matrix:

if tp == 1
    str = 'simple arches';
    topology = 'arches'; st = 2; % 1-2
    
    
elseif tp == 2
    str = 'double arches';
    topology = 'con_arches'; st = 3; % 1-2-3
    
elseif tp == 3
    str = 'triple arches';
    topology = 'con_arches+'; st = 3; % 1-2-3, 1-3
    
elseif tp == 3
    str = '1-2-3-4';
    topology = '1-2-3-4'; st = 4;
end

% 


for i = 2:st:L
    
switch topology
      
    case 'arches' 
        
      EE(i,i-1) = 1; 
      EE(i-1,i) = 1;
  
    case 'con_arches' 
        
      EE(i,i-1) = 1;  EE(i,i+1) = 1; 
      EE(i-1,i) = 1;  EE(i+1,i) = 1; 
      
       case 'con_arches+' 
        
      EE(i,i-1) = 1;  EE(i,i+1) = 1;  EE(i-1,i+1) = 1; 
      EE(i-1,i) = 1;  EE(i+1,i) = 1;  EE(i+1,i-1) = 1; 
      
          case '1-2-3-4' 
        
      EE(i,i-1) = 1;  EE(i,i+1) = 1;  EE(i+1,i+2) = 1; 
      EE(i-1,i) = 1;  EE(i+1,i) = 1;  EE(i+2,i+1) = 1; 
      
      
    case 'rand' 

a = i+1;
b = L-a;

r = ceil(a + (b-a)*rand(1));

     EE(i,r) = 1; 
     EE(r,i) = 1;

end
    
end

EE = EE(1:L, 1:L);
Ee = tril(EE,-1);   
                                                                 
Dat_run = zeros(N*runn,L,td);
   
for nn = 1:runn

 [DAt, ww] = epi('binary',  s0,  0, L, N,  td,  2,   U,  f0, 1 , E_st, Ee);   
%  [DAt, ww] = epi('half',  s0,  0, L, N,  td,  2,   U,  f0, 1 , E_st, Ee, coef);   

a = 1 + (nn-1)*N;
b = N + (nn-1)*N;

disp(nn)

Dat_run(a:b,:,:) = DAt(:,:,:);

end

 %%
                                 
 
   for AVG = [10 25 50 100 runn]     
       
         tmin = 1;
         tmax = td;
         
                                                                           for t = td %tmin:tmax % time of analysis

ss = [1:L];
I = find(Ee);                                                                        
     
 DATA = Dat_run(1:AVG*N,:,t);
 
% add detection code    
    
    Ko = DATA; Ng = N*AVG; 
      
i_div = (1:L); ip = nchoosek(i_div,2); np = size(ip,1);
fhap = zeros(np,4); UFE = zeros(np,1); %WU = zeros(np,1);
UFEnot = zeros(np, L-2);

for i=1:np
    
f00 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 0],2)); f01 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 1],2));
f10 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 0],2)); f11 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 1],2));
                 
fhap(i,:) = [f00 f01 f10 f11]; UFE(i) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); %WU(i) = log((f11*f00)/(f01*f10));

end  
    
UFE(UFE<0) = 0; 
UFE(UFE>=1) = 0; 

                                                                           end
                                                                           

                                                                           
                                                                           
                                                                           
%%

ufe = UFE(UFE > 0.15);
[ac,bc] = hist(ufe);


str = sprintf('Th=%g, E=%g, muL=%g, s_0=%g, f_0=%g  \n N=%d, L=%g, t=%g,  \n DOUBLE ARCHEs', th, E_st, U, s0, f0/L, N, L, td); %disp('a=0, pos epi.'), 


figure(1)
subplot(1,3,ii+1), plot(bc,ac,'LineWidth',2), title(str),  hold on, axis square, set(gca,'FontSize',16),xlabel('UFE'),ylabel('# of sites')
legend('10 POPs', '25 POPs', '50 POPs', '100 POPs', '200 POPs')

   end

                                              end

                                              %%
                                              

stat = horzcat(ip,UFE);     
 
m = 1; mm = m; tt = 1; 
 zz = length(I); n = 7;  
 
pick = ceil(rand(1,200)*np);
Pi = sort(pick);
stat_r = stat(Pi, :); 
 
      [I,J] = find(Ee); 
      
  % --------------------------------------------      

  % Threeshold

   tUF = th;  % prctile(stat_r(:,7),80); 
   
  % ----------------------------------------
  
%   v = tUF; [tUF_opt] = fminsearch(@tUFE,v);
%   
[N2] = find(stat(:,3) > tUF);

raw_det(t).Pairs = N2; % memorize all detected pairs "RAW"
% 
%  det = 0;FP = 0;  c4 = isempty(N2);
%  
%      if c4 == 0
% 
% for n2 = 1 : size(N2,1) 
%     
%     % reconstruct detected matrix for comparison with Ee (EE)
%     
%              if EE(stat(N2(n2),1), stat(N2(n2),2)) == 1
%                  
%                  det_mat(stat(N2(n2),1), stat(N2(n2),2),t) = -1; % detected
%                  
%                  
%                 det = det+1;
%              else
%                  FP = FP+1;
%                  
%                  det_mat(stat(N2(n2),1), stat(N2(n2),2),t) = 1; % false positive
%                  
%              end    
% end
% 
%          det = 100*det/zz; FP = 100*FP/n2;    
%          out22 = [det;FP];
%          
%     else
%               
%       out22 = [det;FP];
%       
%     end
%     
%     % -------------------------
%          
% for ji = 1 : zz
%     
%    det_mat(I(ji), J(ji),t) = -1; % known pairs
%     
% end
%  
%   % -------------------------
%  
% % x PLOT
% 
% Y = horzcat(out22);
%                     
%      OUT_det_UF(1,t) = Y(1,1); %#ok<SAGROW>
%      OUT_fp_UF(1,t) = Y(2,1);  %#ok<SAGROW>
%                                                                      end
% 
% 
% %%
% 
%      Det_UF = mean(OUT_det_UF,3); 
%      
%      Fp_UF = mean(OUT_fp_UF,3); 
%       
% % ----------------------------
% figure(1)
% subplot(1,2,1), plot(Det_UF','LineWidth',2), title('Detection UFE'),  hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax])
% subplot(1,2,2), plot(Fp_UF','LineWidth',2),  title('False Pos. UFE'), hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax])
% 
%     legend off10 25 50 runn]     
       
% 
% 
% 
% 
% % %----------------------------
% % gm = figure(1);
% % 
% % fig_n2 = 'FIGURE';
% % set(gm, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % 
% % set(gm,'Units','Inches');pos = get(gm,'Position');
% % set(gm,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % print(gm,fig_n2,'-dpdf','-r0')
% % 


  %%     
  % Analysis of False-Bonds by triplets f_{110}


Nn =  raw_det(td).Pairs;


DATA = Dat_run(1:AVG*N,:,t);
 

% add detection code    
    
%     Ko = DATA; Ng = N*AVG; 
%      
% DATA = DAt(:,:,tmin); 
%   Ko = DATA;   

Raw = stat(Nn,:);


for k = 1 : L
    chr{k} = num2str(k);
end

chr = chr';

% all sites involved in pairs:

sites = unique(Raw(:,1:2)); %size(sites)
pairs = Raw(:,1:2);

G = graph(pairs(:,1),pairs(:,2),Raw(:,end),chr);

% G = graph(pairs(:,1),pairs(:,2),[],L); %   weights: Raw(:,7),L);

% figure(1)
% subplot(2,2,2), plot(G,'NodeColor','k','EdgeAlpha',0.9,'EdgeLabel',G.Edges.Weight), axis square



                                             end
% %%
% 
% bin = conncomp(G,'Type','weak');
% 
% %[za,zb] = hist(bin,unique(bin))  % unique(bin) number of cluster (including single elements)
% [bin,~,binsize] = histcounts(bin);
% 
% Cls = bin(bin >2)
% 
% figure(4)
% 
% for cl = 1:length(Cls)
%     
% 
% idx = bin(binsize) == Cls(cl);
% SG = subgraph(G, idx); hold on
% subplot(2,2,cl), plot(SG) % separate each subclusters (greater than two sites)
% 
% % identify the sites within sub-clusters > 2 sites
% 
% Sub_Cl = [1:L]; Sub_Cl = Sub_Cl(idx)
% 5
% [~,b1] = intersect(Raw(:,1), Sub_Cl)
% [~,b2] = intersect(Raw(:,2), Sub_Cl)
% 
% %b = mod(bb,size(Raw,1))
% %b(b > size(Raw,1)) = mod(b,size(Raw,1))
% 
% b = union(b1,b2)
% 
% sub_Raw = Rawmin(b',1:2)
% 
% 
% 
% %consider all triplets:
% 
% ii_div = Sub_Cl;
% iip = nchoosek(ii_div,3)
% nt = size(iip,1);
% 
% fhap_iii = zeros(nt,4);
%     UFEiii = zeros(nt,3);
%     
% for i=1:nt
%     
%     for i0 = 1:3
%     
%         if i0 == 1
%         a = 1; b = 2; c = 3;
%         elseif i0 == 2
%         a = 1; b = 3; c = 2;
%         elseif i0 == 3
%         a = 3; b = 2; c = 1;
%         else
%         end
%     
% 
% f00 = mean(all(Ko(:,[iip(i,a), iip(i,b), iip(i,filename = sprintf('%s_r=%g_s0=%g_%g_%g_%g_%g_%g_%g',...
% f01 = mean(all(Ko(:,[iip(i,a), iip(i,b), iip(i,c)]) == ones(N,1)*[0 1 0],2));
% f10 = mean(all(Ko(:,[iip(i,a), iip(i,b), iip(i,c)]) == ones(N,1)*[1 0 0],2)); 
% f11 = mean(all(Ko(:,[iip(5i,a), iip(i,b), iip(i,c)]) == ones(N,1)*[1 1 0],2));
%                
% % fhap(i,:) = [f00 f01 f10 f11]; 
% UFEiii(i,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); 10 25 50 runn]     
       
% 
%     end
% 
% end
%    
% % 
% % UFE(UFE<0) = 0; 
% % UFE(UFE>1) = 0; 
% % 
% % stat = horzcat(ip,fhap,UFE);
% 
% 
% Results = horzcat(iip,UFEiii)
% 
% 
% end
% 
% %%
%%
clc
format short

Bin = conncomp(G,'Type','weak');

%[za,zb] = hist(bin,unique(bin))  % unique(bin) number of cluster (including single elements)
[bin,~,binsize] = histcounts(Bin);

Cls = bin(bin > 2); disp(Cls)

% issue: how to separate small cluster of same size



sb = length(Cls);

for cl = 1%1:sb%length(Cls)
    
    
   clear UFEiii 
    
idx = bin(binsize) == Cls(cl);

ss = [1:L];
ss = ss(idx);

Gsub = subgraph(G, ss); hold on
% subplot(2,sb,cl), 
subplot(1,3,2), sub = plot(Gsub,'NodeColor','k','EdgeAlpha',0.5, 'EdgeLabel',Gsub.Edges.Weight,'LineWidth',2); axis square% separate each subclusters (greater than two sites)


% identify the sites within sub-clusters > 2 sites

Sub_Cl = [1:L]; Sub_Cl = Sub_Cl(idx)
disp(' ')
disp(' ')

[~,b1] = intersect(Raw(:,1), Sub_Cl);
[~,b2] = intersect(Raw(:,2), Sub_Cl);

%b = mod(bb,size(Raw,1))
%b(b > size(Raw,1)) = mod(b,size(Raw,1))

b = union(b1,b2);

sub_Raw = Raw(b',:); 
disp('sub_RAW'), disp(sub_Raw(:,[1 2 3])) % i,j and UFE (3-6 fhap)

Ko = DATA;



%consider all triplets:             %% instead of ALL triplets, target specific links





    nt = size(sub_Raw,1);
%     UFEiii = zeros(nt,3);
        
   test = size(Gsub.Edges,1); % # of bonds to be tested
    
for it = 1 : test
    
    te = setdiff(Sub_Cl,sub_Raw(it,[1 2]));
    test2 = size(te,2);
     
    for i0 = 1:test2  

        a = sub_Raw(it,1); b = sub_Raw(it,2); c = te(i0); 

        %disp(a), disp(b), disp(c)      
        
f00 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 0 0],2)); 
f01 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 1 0],2));
f10 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 0 0],2)); 
f11 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 1 0],2));
               
% fhap(i,:) = [f00 f01 f10 f11]; 
UFEiii(it,i0,cl) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));  %#ok<SAGROW> % UFE for triplets



% 
% te = setdiff(i_div,ip(i,:), 'stable'); % A that is not in B
%     test2 = size(te,2);
%      
%     for i0 = 1:test2
%    
%         a = ip(i,1); b = ip(i,2); c = te(i0);       
%         
% f00 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 0 0],2)); 
% f01 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 1 0],2));
% f10 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 0 0],2)); 
% f11 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 1 0],2));
%                
% UFEiii(i,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));  %#ok<SAGROW> % UFE for triplets

    end

end
   

clear UFEiii_avg tag


% UFE(UFE<0) = 0; 
% UFE(UFE>1) = 0; 
UFEiii_avg  = mean(UFEiii(:,:,cl),2);

disp(' ')
disp(' ')

Results = horzcat(sub_Raw(:,[1 2 3]),UFEiii_avg);
ratio = (Results(:,3)./Results(:,4))-1;

disp('    A         B         UFE_AB   <UFE_ABc>  ratios ')
Results = horzcat(Results, ratio); disp(Results)


% threeshold must be otpimised, right now, quite a random 0.4 without much meaning
TRUE = abs(ratio) < .09; %some ratios are negative

Gsub.Edges.Weight = UFEiii_avg;
% 
% subplot(2,sb,sb+cl), 
subplot(1,3,3), h = plot(Gsub,'NodeColor','k','EdgeAlpha',0.5, 'EdgeLabel',Gsub.Edges.Weight,'LineWidth',2); axis square% separate each subclusters (greater than two sites)

     
tag(TRUE) = {'true'}; 
tag(~TRUE) = {'false'};

labeledge(h,1:numedges(Gsub),tag);

end                                               