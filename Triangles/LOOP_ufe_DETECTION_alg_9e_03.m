clc
close all
clear all %#ok<CLALL>
clearvars 



% Parameters
ii = 0;

                                           L = 50;
                                             
                                                     th = .6;

                                           runn = 200;  
                                              N = 1000;  N1 = N; 
                                              
                                              
                                              for tp = 2%1:4
                                              
  s0 = .1; 
E_st = .75;
  f0 = .45*L;                                                
   U = 0.07; 
  td = 29+1;  
 
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

                                              end
                                              
                                              str = sprintf('Th=%g, E=%g, muL=%g, s_0=%g, f_0=%g  \n N=%d, L=%g, t=%g,  POPs=%g, \n DOUBLE ARCHEs', th, E_st, U, s0, f0/L, N, L, td, runn); %disp('a=0, pos epi.'),

%
 
 
 % Epistatic analysis: UFE_{AB} and UFE_{AB0}
 
   for AVG = runn % @ POPs max
       
         tmin = 1;
         tmax = td;
                                                                           t = td; % @ one time-point

ss = 1:L;I = find(Ee);                                                                        
  
 DATA = Dat_run(1:AVG*N,:,t);
 
% initialise statistics and haplotypes 

Ko = DATA; Ng = N*AVG; 
i_div = (1:L); 
ip = nchoosek(i_div,2); 
np = size(ip,1);
UFE = zeros(np,1); 

for i=1:np   
       
f00 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 0],2)); f01 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 1],2));
f10 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 0],2)); f11 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 1],2));
                 
UFE(i) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); 

end  

   %Filtering (not rocket science)
   
  UFE(UFE <0) = 0; 
  UFE(UFE>=1.5) = 0; 
                                                                 
                                                                                                                                                  
%%

   end

                                             

                                              %%
                                              
                                               tic % Cluster inspection:
                                              
                                              
stat = horzcat(ip,UFE); m = 1; mm = m; tt = 1; zz = length(I); n = 7;  
pick = ceil(rand(1,200)*np); Pi = sort(pick); stat_r = stat(Pi, :); 
 
      [I,J] = find(Ee);     
  % --------------------------------------------      
  % Threeshold

   tUF = th;  % prctile(stat_r(:,7),80); 
 
  % ----------------------------------------
  
[N2] = find(stat(:,3) > tUF); raw_det(t).Pairs = N2; %#ok<SAGROW> % memorize all detected pairs "RAW"
  %%     
  % Analysis of False-Bonds by triplets f_{110}

Nn =  raw_det(td).Pairs;
DATA = Dat_run(1:AVG*N,:,t);  
Raw = stat(Nn,:);

for k = 1 : L
    chr{k} = num2str(k);
end

chr = chr';

% all sites involved in pairs:

sites = unique(Raw(:,1:2)); %size(sites)
pairs = Raw(:,1:2);
G = graph(pairs(:,1),pairs(:,2),Raw(:,3),chr); 
                                           

clc
format short

Bin = conncomp(G);

[bin,~,binsize] = histcounts(Bin);

Cls = bin(bin > 2); disp(Cls)

% note: does not separate cluster of same size

sb = length(unique(Cls));
sup = 0;

for cl = 1:sb
    
   clear UFEiii fhap Tt Ff T F
    
idx = bin(binsize) == Cls(cl);

ss = [1:L];
ss = ss(idx);

Gsub = subgraph(G, ss); hold on
subplot(sb,2,1+sup), sub = plot(Gsub,'NodeColor','k','EdgeAlpha',0.5, 'EdgeLabel',Gsub.Edges.Weight,'LineWidth',2); axis square, hold on, % separate each subclusters (greater than two sites)
xlabel('Raw UFE_{AB}')

% identify the sites within sub-clusters > 2 sites
Sub_Cl = [1:L]; Sub_Cl = Sub_Cl(idx);

[~,b1,~] = intersect(Raw(:,1), Sub_Cl);
[~,b2,~] = intersect(Raw(:,2), Sub_Cl);

bU = union(b1,b2);
A = Sub_Cl; B = Raw(:,1)';
c=ismember(B,A);
sub_Raw = Raw(c',:);


%consider all triplets:             %% instead of ALL triplets, target specific links


     nt = size(sub_Raw,1);        
   test = size(Gsub.Edges,1); % # of bonds to be tested 
  
   clc

for it = 1 : test
    
    te = setdiff(Sub_Cl,sub_Raw(it,[1 2]));
    test2 = size(te,2);
     
    for i0 = 1:test2  

        a = sub_Raw(it,1); b = sub_Raw(it,2); c = te(i0); 
sites = [a b c];
        disp(sites),      
        
       
        
f00 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 0 0],2)); 
f01 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 1 0],2));
f10 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 0 0],2)); 
f11 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 1 0],2));

 fhap(it,:) = [f00 f01 f10 f11];   
 
UFEiii(it,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));

disp(UFEiii(it,i0));

disp('---') 
% pause


% UFE for triplets

    end

end
   
%   UFEiii = abs(UFEiii); 
  
  UFEiii(UFEiii <0) = 0; 
  UFEiii(UFEiii>=1.5) = 0; 

  
  
UFEiii_min  = min(UFEiii,[],2);

disp(' ')

Results = horzcat(sub_Raw(:,[1 2 3]),UFEiii_min);
Res_ii0 = horzcat(sub_Raw(:,[1 2 ]),UFEiii);
Res_ii0_fhap = horzcat(sub_Raw(:,[1 2]),fhap,UFEiii);
prufe = (Results(:,3)./Results(:,4)); %product of ufe AB x AB0

prufe(isinf(prufe)) = 10; %arbitrarily set to 10 fold reduction (because UFE_{AB0} = 0)

disp('    A         B         UFE_AB   min(UFE_AB0)  ratio ')
Results = horzcat(Results, prufe); disp(Results) %#ok<AGROW>

% threeshold must be otpimised:
TRUE = abs(prufe) < 2; % UFE_{AB}/UFE_{AB0} > 2

GsuB = subgraph(G, ss);
GsuB.Edges.Weight = UFEiii_min;

subplot(sb,2,sup+2), h = plot(GsuB,'NodeColor','k','EdgeColor','r','EdgeAlpha',0.5, 'EdgeLabel',GsuB.Edges.Weight,'LineWidth',2); xlabel('min(UFE_{AB0})'), axis square, hold on % separate each subclusters (greater than two sites)

if cl == 1
    title(str)
else 
end

tag(TRUE) = {'true'}; 
tag(~TRUE) = {'false'};

sup = sup +2;
%%

T = TRUE;
F = ~TRUE;

T = sub_Raw(T, [1 2]); 
F = sub_Raw(F, [1 2]); 

% This section recover the position of "ith" Node, to reorder the Nodes
% for TRUE/FALSE coloring

if size(T,1) == 0  &&  size(F,1) == 0

disp('no TRUE pairs detected')

else
    
Tt = T(1,:); Ff = F(1,:);

for sp = 2 : size(T,1)  
    Tt = horzcat(Tt,T(sp,:));
end
for sp = 2 : size(F,1)  
    Ff = horzcat(Ff,F(sp,:));
end

Tt_u = unique(Tt,'Stable'); Tt_0 = Tt; Tt_0(:) = 0;
C = table2cell(GsuB.Nodes)';

for uu = 1 : length(Tt_u)
%     num2str(Tt_u(uu))
%     C(uu)
    
index = find(strcmp(num2str(Tt_u(uu)), C));
ind = find(Tt == Tt_u(uu));
Tt_0(1,ind) = index; 
end


% highlight(h,[Ff_0],'EdgeColor','r','LineWidth',1.5)
highlight(h,[Tt_0],'EdgeColor','g','LineWidth',1.5)
end
%%pause

end    


%%
% Res_ii0
% Res_ii0_fhap
toc


%%

% % % SAVE FIGs 
% gm = figure(1);
% 
%  fig_n2= ['fig_Cluster_UFE_AB0.pdf'];
% 
% set(gm, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% 
% set(gm,'Units','Inches');pos = get(gm,'Position');
% set(gm,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gm,fig_n2,'-dpdf','-r0')
% close(gm)