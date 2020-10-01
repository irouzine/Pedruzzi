clc
close all
clear all %#ok<CLALL>
clearvars 


tic

% Parameters
ii = 0;

                                             % L = 10;
                                             for L = [50]
                                             
%                                                  if ii == 0 
                                                     th = .5;
%                                                  elseif ii == 2
% %                                                      th = .6;
%                                                  else
%                                                  end
                                                 
                                           runn = 100;  
                                              N = 1000;  N1 = N; 
                                              
                                              
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
indi_EE = EE;
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
      
      indi_EE(i-1,i+1) = 1;
      indi_EE(i+1,i-1) = 1;
      
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
                                 
 
   for AVG = runn%[10 25 50 100 runn]     
       
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

index_true = zeros(np,1);
index_indi = zeros(np,1);

for i=1:np   
    
    
     if EE(ip(i,1), ip(i,2)) == 1 % track true pairs
index_true(i,1) = 1;
     else
     end

     
     if indi_EE(ip(i,1), ip(i,2)) == 1 % track INDI pairs, for duble arch topo.
index_indi(i,1) = 1;
     else
     end
     
     
    
f00 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 0],2)); f01 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 1],2));
f10 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 0],2)); f11 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 1],2));
                 
fhap(i,:) = [f00 f01 f10 f11]; UFE(i) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); %WU(i) = log((f11*f00)/(f01*f10));


te = setdiff(i_div,ip(i,:), 'stable'); % A that is not in B
    test2 = size(te,2);
     
    for i0 = 1:test2
   
        a = ip(i,1); b = ip(i,2); c = te(i0);       
        
f00 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 0 0],2)); 
f01 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 1 0],2));
f10 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 0 0],2)); 
f11 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 1 0],2));
               
UFEiii(i,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));  % UFE for triplets

    end

end  
    
  UFE(UFE<0) = 0; 
UFE(UFE>=1) = 0; 

  UFEiii(UFEiii<0) = 0; 
UFEiii(UFEiii>=1) = 0;      
                        
UFEiii_avg  = mean(UFEiii(:,:),2);
UFEiii_min  = min(UFEiii(:,:),[],2);
                                                                           end
                                                                           
ind_true = find(index_true); UFE_true = UFE(ind_true,:);  UFEiii_true_avg = UFEiii_avg(ind_true,:); UFEiii_true_min = UFEiii_min(ind_true,:);                                                                           
ind_indi = find(index_indi); UFE_indi = UFE(ind_indi,:);  UFEiii_indi_avg = UFEiii_avg(ind_indi,:); UFEiii_indi_min = UFEiii_min(ind_indi,:); 

 %   ufe_t = UFE_true(UFE_t > 0.15);                                                                         
                                                                           
%%

ufe = UFE(UFE > 0.15);

[ac,bc] = hist(ufe);

[at,bt] = hist(UFE_true); [ati,bti] = hist(UFEiii_true_avg); [atm,btm] = hist(UFEiii_true_min);


[ai,bi] = hist(UFE_indi); [aii,bii] = hist(UFEiii_indi_avg); [aim,bim] = hist(UFEiii_indi_min); 


str = sprintf('Th=%g, E=%g, muL=%g, s_0=%g, f_0=%g  \n N=%d, L=%g, t=%g,  POPs=%g, \n DOUBLE ARCHEs', th, E_st, U, s0, f0/L, N, L, td, runn); %disp('a=0, pos epi.'), 

% figure(1)
% subplot(1,3,ii+1), plot(bc,ac,'LineWidth',2), title(str),  hold on, axis square, set(gca,'FontSize',16),xlabel('UFE'),ylabel('# of pairs')
% legend('10 POPs', '25 POPs', '50 POPs', '100 POPs', '200 POPs')

figure(2)
subplot(2,2,1), plot(bt,at,'r','LineWidth',2), title(str),  hold on, axis square, set(gca,'FontSize',16),xlabel('UFE_{ij}'),ylabel('# of pairs')
subplot(2,2,1), plot(bi,ai,'b','LineWidth',2), hold on, legend('red TRUE', 'blue INDI')

subplot(2,2,2), plot(bti,ati,'r','LineWidth',2), hold on, axis square, set(gca,'FontSize',16),xlabel('<UFE_{ij0}>'),ylabel('# of pairs')
subplot(2,2,2), plot(bii,aii,'b','LineWidth',2), hold on, 

subplot(2,2,3), plot(btm,atm,'r','LineWidth',2), hold on, axis square, set(gca,'FontSize',16),xlabel('min(UFE_{ij0})'),ylabel('# of pairs')
subplot(2,2,3), plot(bim,aim,'b','LineWidth',2), hold on,



%%

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
  
[N2] = find(stat(:,3) > tUF);

raw_det(t).Pairs = N2; % memorize all detected pairs "RAW"


%   %%     
%   % Analysis of False-Bonds by triplets f_{110}
% 
% Nn =  raw_det(td).Pairs;
% 
% DATA = Dat_run(1:AVG*N,:,t);  
% 
% Raw = stat(Nn,:);
% 
% 
% for k = 1 : L
%     chr{k} = num2str(k);
% end
% 
% chr = chr';
% 
% % all sites involved in pairs:
% 
% sites = unique(Raw(:,1:2)); %size(sites)
% pairs = Raw(:,1:2);
% 
% G = graph(pairs(:,1),pairs(:,2),Raw(:,3),chr);
% 
% % G = graph(pairs(:,1),pairs(:,2),[],L); %   weights: Raw(:,7),L);
% 
% % figure(1)
% % subplot(2,2,2), plot(G,'NodeColor','k','EdgeAlpha',0.9,'EdgeLabel',G.Edges.Weight), axis square



                                             end


%%
% 
% clc
% format short
% 
% Bin = conncomp(G,'Type','weak');
% 
% %[za,zb] = hist(bin,unique(bin))  % unique(bin) number of cluster (including single elements)
% [bin,~,binsize] = histcounts(Bin);
% 
% Cls = bin(bin > 2); disp(Cls)
% 
% % issue: how to separate small cluster of same size
% 
% 
% 
% sb = length(Cls);
% 
% for cl = 1%1:sb%length(Cls)
%     
%     
%    clear UFEiii 
%     
% idx = bin(binsize) == Cls(cl);
% 
% ss = [1:L];
% ss = ss(idx);
% 
% Gsub = subgraph(G, ss); hold on
% % subplot(2,sb,cl), 
% figure(1)
% subplot(1,3,2), sub = plot(Gsub,'NodeColor','k','EdgeAlpha',0.5, 'EdgeLabel',Gsub.Edges.Weight,'LineWidth',2); axis square% separate each subclusters (greater than two sites)
% 
% 
% % identify the sites within sub-clusters > 2 sites
% %
% Sub_Cl = [1:L]; Sub_Cl = Sub_Cl(idx)
% disp(' ')
% 
% [~,b1,~] = intersect(Raw(:,1), Sub_Cl)
% [~,b2,~] = intersect(Raw(:,2), Sub_Cl)
% 
% bU = union(b1,b2)
% 
% 
% 
% A = Sub_Cl
% B = Raw(:,1)'
% 
% c=ismember(B,A)
% % C=A(c,:)
% % D=A(~c,:)
% sub_Raw = Raw(c',:)
% 
% %
% % sub_Raw = Raw(bU',:)
% disp('sub_RAW'), disp(sub_Raw(:,[1 2 3])) % i,j and UFE (3-6 fhap)
% 
% Ko = DATA;
% 
% 
% 
% %consider all triplets:             %% instead of ALL triplets, target specific links
% 
% 
% 
% 
% 
%     nt = size(sub_Raw,1);
% %     UFEiii = zeros(nt,3);
%         
%    test = size(Gsub.Edges,1); % # of bonds to be tested
%     
% for it = 1 : test
%     
%     te = setdiff(Sub_Cl,sub_Raw(it,[1 2]));
%     test2 = size(te,2);
%      
%     for i0 = 1:test2  
% 
%         a = sub_Raw(it,1); b = sub_Raw(it,2); c = te(i0); 
% 
%         %disp(a), disp(b), disp(c)      
%         
% f00 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 0 0],2)); 
% f01 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 1 0],2));
% f10 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 0 0],2)); 
% f11 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 1 0],2));
%                
% % fhap(i,:) = [f00 f01 f10 f11]; 
% UFEiii(it,i0,cl) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));  %#ok<SAGROW> % UFE for triplets
% 
%     end
% 
% end
%    
% 
% clear UFEiii_avg tag
% 
% 
% % UFE(UFE<0) = 0; 
% % UFE(UFE>1) = 0; 
% UFEiii_avg  = mean(UFEiii(:,:,cl),2);
% 
% disp(' ')
% disp(' ')
% 
% Results = horzcat(sub_Raw(:,[1 2 3]),UFEiii_avg);
% ratio = (Results(:,3)./Results(:,4))-1;
% 
% disp('    A         B         UFE_AB   <UFE_ABc>  ratios ')
% Results = horzcat(Results, ratio); disp(Results)
% 
% 
% 
% 
% 
% % threeshold must be otpimised, right now, quite a random value without much meaning
% TRUE = abs(ratio) < .09; %some ratios are negative
% 
% 
% 
% 
% 
% GsuB = subgraph(G, ss);
% GsuB.Edges.Weight = UFEiii_avg;
% % 
% figure(1)
% % subplot(2,sb,sb+cl), 
% subplot(1,3,3), h = plot(GsuB,'NodeColor','k','EdgeColor','r','EdgeAlpha',0.5, 'EdgeLabel',GsuB.Edges.Weight,'LineWidth',2); axis square% separate each subclusters (greater than two sites)
% 
%      
% tag(TRUE) = {'true'}; 
% tag(~TRUE) = {'false'};
% 
% labeledge(h,1:numedges(GsuB),tag);
% 
% %%
% 
% T = find(TRUE);
% F = find(~TRUE);
% 
% T = sub_Raw(T, [1 2]); 
% F = sub_Raw(F, [1 2]); 
% 
% % This section recover the position of "ith" Node, to reorder the Nodes
% % for TRUE/FALSE coloring
% 
% Tt = T(1,:); Ff = F(1,:);
% 
% for sp = 2 : size(T,1)  
%     Tt = horzcat(Tt,T(sp,:));
% end
% for sp = 2 : size(F,1)  
%     Ff = horzcat(Ff,F(sp,:));
% end
% 
% Tt_u = unique(Tt,'Stable'); Tt_0 = Tt; Tt_0(:) = 0;
% % Ff_u = unique(Ff,'Stable'); Ff_0 = Ff; Ff_0(:) = 0;
% C = table2cell(GsuB.Nodes)'
% 
% for uu = 1 : length(Tt_u)
% %     num2str(Tt_u(uu))
% %     C(uu)
%     
% index = find(strcmp(num2str(Tt_u(uu)), C));
% ind = find(Tt == Tt_u(uu));
% Tt_0(1,ind) = index; 
% end
% %Tt_0
% 
% % for vv = 1 : length(Ff_u)
% %     num2str(Ff_u(vv))
% %     C
% %     
% % index = find(strcmp(num2str(Ff_u(vv)), C))
% % ind = find(Ff == Ff_u(vv))
% % Ff_0(1,ind) = index
% % end
% % Ff_0
% 
% % highlight(h,[Ff_0],'EdgeColor','r','LineWidth',1.5)
% highlight(h,[Tt_0],'EdgeColor','g','LineWidth',1.5)
% 
% %%
% 
% end                                               