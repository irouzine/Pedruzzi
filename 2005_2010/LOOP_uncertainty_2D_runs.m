clc
clear all %#ok<CLALL>
                    
% =========================================================================

% Sensitivity of Results w.r.t. detection parameters:

clc; clear all; close all;  

[Headers, Sequences] = multialignread('full_align_na_2010.aln');
load Data_na_2010;
% [Headers, Sequences] = multialignread('full_align_ha.aln');
% load Data_ha;

mK = mean(mean(K)); 


%%
%Preliminary Anlysis

tic
  
  % --------------------------------------------  
  
  % Threesholds and Parameters

 
   %Dw = .2;  % downweighting coeff.
   %dv = 0.05; % lower-limit for sequence similarity
   
   
   
    tU = .6;   % raw UFE detection threeshold
    
    fmin = 0.005; % observably diverse sites
    Fmin = 0.01;  % observably diverse pairs
   
    RATIO = 2;    % UFE/UFE_{AB00 ratio}
    
    minUfe = -1;  % UFE filtering:
    maxUfe = 1.5;
    
  % ----------------------------------------
  
  in = 1;
  
  runs = 4;
  
  for run = 1:runs
      
  yi = 1;
  
for Dw = .1 : .01 : .5
    
    
 xi = 1;
 
  for dv = .01 : .001 : .05
     
       disp(in) 
  % downweighting sequences with less than dv% mutations
  
f_i = mean(K, 1);f_n = mean(K, 2); 
f_mut = sum(K, 2); uu_th = find(f_n >= dv);
Km = K(uu_th,:); 

u0 = find(f_n < dv); 
K0 = K(u0,:); 

upDw = 1+ceil(size(u0,1)*Dw);

r = randi([1 size(u0,1)],1,upDw);

K0 = K0(r,:);

Kdw = vertcat(K0,Km);



% Calculate Haplotypes f_{i,j} and Epis. Statistics (UFE, Wu)
         
Ko = Kdw(1:end,:); N =  size(Ko,1); %size(uu,1); 
L = pos; mapU = zeros(L,L); mapW = mapU; 
i_div = (1:L); ip = nchoosek(i_div,2);  np = size(ip,1);
fhap = zeros(np,4); UFE = zeros(np,1);  WU = zeros(np,1); f_i1 = zeros(np,1); f_i2 = zeros(np,1); fhap_m = zeros(np,1);
  
mKo = mean(mean(Ko));

parfor ij=1:np
%  disp(ij/np)
    
f00 = mean(all(Ko(:,[ip(ij,1), ip(ij,2)]) == ones(N,1)*[0 0],2)); 
f01 = mean(all(Ko(:,[ip(ij,1), ip(ij,2)]) == ones(N,1)*[0 1],2));
f10 = mean(all(Ko(:,[ip(ij,1), ip(ij,2)]) == ones(N,1)*[1 0],2)); 
f11 = mean(all(Ko(:,[ip(ij,1), ip(ij,2)]) == ones(N,1)*[1 1],2));
                 
fhap(ij,:) = [f00 f01 f10 f11]; fhap_m (ij) = min(fhap(ij,:));
f_i1(ij) = f_i(ip(ij,1));f_i2(ij) = f_i(ip(ij,2));

UFE(ij) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); 

end

TABLE = horzcat(ip,f_i1,f_i2,fhap,UFE,fhap_m);  

ii=TABLE(:,3) > fmin & TABLE(:,3) < 1-fmin;     % observably diverse sites
TAB_i = TABLE(ii,:);

iii=TAB_i(:,10) > Fmin & TAB_i(:,10) < 1-Fmin;  % observably diverse pairs
TAB_ij = TAB_i(iii,:);

                                              
%  % Epistatic analysis: UFE_{AB} and UFE_{AB0}

   %Filtering (not rocket science)

UFE = TAB_i(iii,9);
  
UFE(UFE < minUfe) = 0; 
UFE(UFE >=maxUfe) = 0; %Topology-dependent value                                                                                                     
                                              
                                               % Cluster inspection:                           

stat = horzcat(TAB_ij(:,[1 2]),UFE);
[N2] = find(stat(:,3) > tU); raw_det.Pairs = N2; % memorize all detected pairs "RAW"
       
  % Analysis of False-Bonds by triplets f_{110}

Nn =  raw_det.Pairs;
Raw = stat(Nn,:);

for k = 1 : L
    chr{k} = num2str(k);
end
chr = chr';

% all sites involved in pairs:

G = graph(Raw(:,1),Raw(:,2),Raw(:,3),chr); 
                                           
[Bin, Binsize] = conncomp(G);
[bin,~,binsize] = histcounts(Bin);
Cls = bin(bin > 2); %disp(Cls)

% note: does not separate cluster of same size

Ng = N;
sb = 1;
cl = 1;     % Analysis of FIRST Cluster of dim > 2
    
   clear UFEiii fhap Tt Ff T F idx ss Gsub sub_Raw
    
idx = bin(binsize) == Cls(cl);

ss = [1:L];
ss = ss(idx);

Gsub = subgraph(G, ss); hold on

% identify the sites within sub-clusters > 2 sites
Sub_Cl = 1:L; Sub_Cl = Sub_Cl(idx);

[~,b1,~] = intersect(Raw(:,1), Sub_Cl);
[~,b2,~] = intersect(Raw(:,2), Sub_Cl);

bU = union(b1,b2);
A = Sub_Cl; B = Raw(:,1)';
c=ismember(B,A);
sub_Raw = Raw(c',:);                %% extract Sub_RAW only for sub_Cluster

% -----------------------------------------------------------------------------------%
%                          %% instead of ALL triplets, target specific links


     nt = size(sub_Raw,1);        
   test = size(Gsub.Edges,1);       % # of bonds to be tested 
  fhap = zeros(test,4);
   
UFEiii_min = zeros(test,1);

for it = 1 : test
    
iio = ismember(sub_Raw(:,1),sub_Raw(it,[1 ]),'rows')'; %1st site in 1st column
ijo = ismember(sub_Raw(:,1 ),sub_Raw(it,[2 ]),'rows')'; %2nd site in 1st column
Ii = or(iio,ijo); % rows of Sub_raw (1st col) that contains nodes i or j

jio = ismember(sub_Raw(:,2),sub_Raw(it,[1 ]),'rows')'; %1st site in 2nd column
jjo = ismember(sub_Raw(:,2),sub_Raw(it,[2 ]),'rows')'; %2nd site in 2nd column
Jj = or(jio,jjo); % rows of Sub_raw (2st col) that contains nodes i or j

IJ = or(Ii, Jj); % set rows containing i and j nodes.

i_j = sub_Raw(it,[1 2]);
nodes = sub_Raw(IJ,[1 2]); Nodes = unique(reshape(nodes,1,[]));

te = setdiff(Nodes(:,:),i_j);
  test2 = size(te,2);   
    
    UFEiii = zeros(1,test2);
    a = sub_Raw(it,1); b = sub_Raw(it,2);
    
    parfor i0 = 1:test2  
     
         c = te(i0); 
        
% sites = [a b c];disp(sites)      
%      a = 13; b = 14; c = 15; %to check haplotype for specif triplet

f00 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 0 0],2));
f01 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 1 0],2));
f10 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 0 0],2));
f11 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 1 0],2));

%  fhap(it,:) = [f00 f01 f10 f11];
 %
UFEiii(1,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));

    end

     UFEiii(UFEiii <0) = 0; 
%   UFEiii(UFEiii>=1.5) = 0; 

 UFEiii_min(it,1)  = min([UFEiii(1,:), 777],[],2); % cause some pairs are "terminal"

end
  
prufe = (sub_Raw(:,[3])./UFEiii_min);

prufe(isinf(prufe)) = 10; %arbitrarily set to 10 fold reduction (because UFE_{AB0} = 0)
Results = horzcat(sub_Raw(:,[1 2 3]),UFEiii_min,prufe);

% threeshold must be otpimised:
                                                                               TRUE = abs(prufe) < RATIO; % UFE_{AB}/UFE_{AB0} > 2

GsuB = subgraph(G, ss);
GsuB.Edges.Weight = UFEiii_min;

tag(TRUE) = {'true'}; 
tag(~TRUE) = {'false'};

T = TRUE;
F = ~TRUE;

T = sub_Raw(T, [1 2]); 
F = sub_Raw(F, [1 2]); 

zz = GsuB.Edges;
zz = zz(TRUE,:); %% true pairs

edg = table2array(zz(:,[1])); Edg = str2double(edg);
Wei = table2array(zz(:,[2]));

TAB_ij = horzcat(Edg,Wei);

L = max(TABLE(:,2)); 
xU = zeros(L,L); 

i_div = (1:L); ip = nchoosek(i_div,2);  np = size(ip,1);

for i=1:size(TAB_ij,1)
     
xU(TAB_ij(i,1),TAB_ij(i,2)) = TAB_ij(i,3);
xU(TAB_ij(i,2),TAB_ij(i,1)) = TAB_ij(i,3);

end

                                                                                thresh = 0.2;
                                                                                
xU(xU >  -thresh & xU <  thresh) = 0;
xU(xU <  thresh) = 0;

xU(isnan(xU)) = 0;
xU(isinf(xU)) = 0;

toc

xU_i(:,:,in,run) = xU;



dv_i(xi) = dv;

Ufe(xi,yi,run) = size(find(xU>0),1)/2;

% dv
% Dw
% Ufe(xi,yi)


in = in+1; 
xi = xi+1;

  end
  Dw_i(yi) = Dw; %#ok<SAGROW>
 yi = yi+1;
end

  end
%%

Ufe = mean(Ufe,3); save('Plateau_na_2010_runs.mat', 'dv_i', 'Dw_i', 'Ufe', 'xU_i')
    if in > 2

%% 

load Plateau_na_2010_runs.mat
   colormap jet
    figure(33)

% subplot(1,2,1),
surf(Dw_i,dv_i,Ufe), xlabel('Dw'), ylabel('dv'), title('detection plateau (mean(runs))'),  set(gca,'FontSize', 12), 
colorbar, %xlim([min(dv_i) max(dv_i)]), ylim([min(Dw_i) max(Dw_i)])

% subplot(1,2,2),
% surf(Ufe), xlabel('dv'), ylabel('Dw'), title('detection plateau'),  set(gca,'FontSize', 12,'Ydir','reverse'); 
% colorbar, %xlim([min(dv_i) max(dv_i)]), ylim([min(Dw_i) max(Dw_i)])
%%




else;end

%%


% %%
% figure
% subplot(2,2,1), hist(f_i,30), axis square, title('mut. rate per site'), 
% subplot(2,2,2), hist(f_n,30), axis square, title('mut. rate per sequence'), 
% subplot(2,2,3), plot(f_i), axis square, title('mut. rate per site'), 
% subplot(2,2,4), plot(f_n), axis square, title('mut. rate per sequence'), 
% 
% %%
% ma = multialignread('full_align_na_2010.aln');
% % showalignment(ma)
% seqalignviewer(ma(1:100))









