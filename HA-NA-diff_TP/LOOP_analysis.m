close all;  
clear all

Col = 'w'; %tits = {'A','B','C','D','E'};
          tits = {'J','H','I','J','K'};
                                                                                runs = 1; 



for col = 1:9 
    
Dw_cal_1 = [.005 .005 .005 .01 .01 .01 .02 .02 .02];  %[.075  .275   .175   .075  .275]; 
dv_cal_1 = [.05 .1 .15 .05 .1 .15 .05 .1 .15];     %[.015 .015  .0235   .032   .032]; 

% Dw_cal_2 = [.05 .2 .25 .3 .3];  
% dv_cal_2 = [.5 .005 .01 .015 .005];  

for tim = 1
    %1:8

    clc
    clearvars -except runs tim Dw_cal_1 Dw_cal_2 dv_cal_1 dv_cal_2 col Col tits
    load 2D_runs_rev_a;
    
    if tim == 1    
     [Headers, Sequences] = multialignread('full_align_NA_05_08.aln');
    load Data_NA_05_08; str = ('NA H1N1 2005-2008'); fig = 50; Tim = 1; 
    
    elseif tim == 2
     [Headers, Sequences] = multialignread('full_align_NA_09_10.aln');
    load Data_NA_09_10; str = ('NA H1N1 2009-2010'); fig = 60; Tim = 2; 
    
    end
                    
mK = mean(mean(K)); 
%%

figure(100*tim)
subplot(2,1,1), histogram(mean(K),50,'Normalization','probability'),ylabel('distribution f_i (per site)'), title(str); 
subplot(2,1,2), histogram(mean(K,2)',50,'Normalization','probability'),ylabel('distribution f_i (per genome)'), title(str);

%% Plateau and AOI for dv & Dw

V = Ufe_(:,:,tim);
M = size(V,2); steps = 99; 
[X,Y] = meshgrid(1:M);

% Create the query grid with spacing "steps+1"
dvq = linspace(min(dv_),max(dv_),steps+1);
dwq = linspace(min(Dw_),max(Dw_),steps+1);

[Xq,Yq] = meshgrid(1:(M-1)/steps:M);

[X1,Y1] = meshgrid(dvq,dwq);

% Interpolate at the query points.
Vq = interp2(X,Y,V,Xq,Yq,'linear');
                                                            minDet = 0;%prctile(mean(Vq),05); 
           
V_0 = Vq; V_0(V_0 < minDet) = NaN; 


% normalize between [0-1]:
Vq = Vq./(max(max(Vq)));

colormap(parula)

% figure(fig)
% % subplot(2,3,6), 
% s1 = surf(dvq,dwq,Vq); grid on, view([180 90]), axis tight, hold on %set(gca, 'Zdir', 'reverse')
% s1.EdgeColor = 'none'; axis square, xlim([.01 .05])
% title(str), set(gca, 'FontSize',16), grid on, 
% colorbar
%%
figure(234),%subplot(1,2,tim),
contourf(dv_,Dw_,Ufe_(:,:,tim),'--'), ylabel('Dw'), xlabel('dv'), ylim([0.005 .1])
title(str), set(gca,'FontSize', 12), hold on;

%%

      if col == 1
strR = 'Num. detected links';
ht = text(0.05, 0.05,strR); hold on,
set(ht,'Rotation',90)

% set(as,'FontSize',18)
else; end

xlabel('Cut-off d_v')
ylabel('Downweighting D_w')
% colorbar
% zlabel('detected pairs')    

%%

%Preliminary Anlysis

tic
  
  % --------------------------------------------  
  
  % Threesholds and Parameters


    tU = .6;   % raw UFE detection threeshold
    fmin = 0.005; % observably diverse sites
    Fmin = 0.01;  % observably diverse pairs
                                                                             RATIO = 2.5;    % UFE/UFE_{AB00 ratio}
    minUfe = -1;  % UFE filtering:
    maxUfe = 1.5;
    
  % ----------------------------------------
  in = 1;
  %%
                                                                                                                                        
                                                                           
L = pos; mapU = zeros(L,L); mapW = mapU; 
i_div = (1:L); ip = nchoosek(i_div,2); np = size(ip,1);

fhap_ = zeros(np,4,runs); %UFE_ = zeros(np,runs);  
f_i1_ = zeros(np,runs); f_i2_ = zeros(np,runs); fhap_m_ = zeros(np,runs);
 UFE = zeros(np,1);  f_i1 = zeros(np,1); f_i2 = zeros(np,1); %fhap_m = zeros(np,1);
                                                                               
 for run = 1:runs
     
     
if tim == 1
    
        Dw = Dw_cal_1(col);%*(1+r1);
        dv = dv_cal_1(col);%*(1+r2);

elseif tim == 2
    
        Dw = Dw_cal_2(col);%*(1+r1);
        dv = dv_cal_2(col);%*(1+r2);
end       

%          xx = char(tits(col));

figure(fig)
% subplot(2,3,6), 


Li = strsplit(sprintf('%c\n','A':'E'));  % Letter Labels

% plot(x, y, '+w')

% axis([-0.5  10.5    -0.5  10.5 ])

% s4a = scatter3(dv,Dw,200,100,'xw','LineWidth', 2); hold on, set(gca,'FontSize', 18);

figure(234)
% subplot(1,2,tim), 
scatter3(dv,Dw,200,100,'xw','LineWidth', 2); hold on, set(gca,'FontSize', 20);
Te = text(dv, Dw, Li(col), 'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 20, 'Color','w'); Te.Position = [dv+0.002 Dw 200];



  % downweighting sequences with less than dv% mutations
  
f_i = mean(K, 1);f_n = mean(K, 2); f_mut = sum(K, 2); uu_th = find(f_n >= dv);
Km = K(uu_th,:); 

u0 = find(f_n < dv); 
K0 = K(u0,:); 

upDw = 1+ceil(size(u0,1)*Dw);

r = randi([1 size(u0,1)],1,upDw);

K0 = K0(r,:);

Kdw = vertcat(K0,Km);

uu = find(f_n >= dv); %disp(size(uu,1))
u0 = find(f_n < dv);  %disp(size(u0,1))

r = randi([1 size(u0,1)],1,size(uu,1));



%% Calculate Haplotypes f_{i,j} and Epis. Statistics (UFE, Wu)
         
Ko = Kdw(1:end,:); N =  size(Ko,1); %size(uu,1); 
  
mKo = mean(mean(Ko));

for ij=1:np
    
f00 = mean(all(Ko(:,[ip(ij,1), ip(ij,2)]) == ones(N,1)*[0 0],2)); 
f01 = mean(all(Ko(:,[ip(ij,1), ip(ij,2)]) == ones(N,1)*[0 1],2));
f10 = mean(all(Ko(:,[ip(ij,1), ip(ij,2)]) == ones(N,1)*[1 0],2)); 
f11 = mean(all(Ko(:,[ip(ij,1), ip(ij,2)]) == ones(N,1)*[1 1],2));
                 
fhap_(ij,:,run) = [f00 f01 f10 f11]; fhap_m_ (ij,run) = min(fhap_(ij,:,run));

f_i1_(ij,run) = f_i(ip(ij,1));
f_i2_(ij,run) = f_i(ip(ij,2));

%UFE_(ij,run) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); 

end

 end
        
% averaged after "runs"

fhap = mean(fhap_,3); fhap_m = mean(fhap_m_,2);

f_i1(:) = mean(f_i1_,2);
f_i2(:) = mean(f_i2_,2);

for ij=1:np % calculate UFE after <f_ij>_runs
    
f00 = fhap(ij,1); f01 = fhap(ij,2);  f10 = fhap(ij,3); f11 = fhap(ij,4);       
UFE(ij) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); 

end


% record into TABLE: > run UFE_{AB0} only once!

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

stat = horzcat(TAB_ij(:,[1 2]),UFE); m = 1; mm = m; tt = 1; n = 7;
[N2] = find(stat(:,3) > tU); raw_det.Pairs = N2; % memorize all detected pairs "RAW"
       
  % Analysis of False-Bonds by triplets f_{110}

Nn =  raw_det.Pairs;
Raw = stat(Nn,:);

for k2 = 1 : L
    chr{k2} = num2str(k2);
end
chr = chr';

% all sites involved in pairs:

G = graph(Raw(:,1),Raw(:,2),Raw(:,3),chr); 
                                           
[Bin, Binsize] = conncomp(G);
[bin,~,binsize] = histcounts(Bin);
Cls = bin(bin > 2); %disp(Cls)

% note: does not separate cluster of same size

Ng = N;
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

% -------------------------------------------------------------------------%
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
Jj = or(jio,jjo); % rows of Sub_raw (2st col) that cont2ains nodes i or j

IJ = or(Ii, Jj); % set rows containing i and j nodes.

i_j = sub_Raw(it,[1 2]);
nodes = sub_Raw(IJ,[1 2]); Nodes = unique(reshape(nodes,1,[]));

te = setdiff(Nodes(:,:),i_j);
  test2 = size(te,2);   
    
    UFEiii = zeros(1,test2);
    
    for i0 = 1:test2  
     
        a = sub_Raw(it,1); b = sub_Raw(it,2); c = te(i0); 
        
% sites = [a b c];disp(sites)      
%      a = 13; b = 14; c = 15; %to check haplotype for specif triplet

f00 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 0 0],2));
f01 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 1 0],2));
f10 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 0 0],2));
f11 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 1 0],2));

 fhap(it,:) = [f00 f01 f10 f11];
 %
UFEiii(1,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));

    end

     UFEiii(UFEiii <0) = 0; 

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


% start - end 

a = 1; b = pos;
xu = xU; %xu = median(xU_i,3);
xu = xu(a:b,a:b);

%Create custom node labels
myLabel = cell(length(xu));

for i = 1:length(xu)
    
    if sum(xu(:,i)) == 0
    
  myLabel{i} = ' ';%num2str(a-1+i);
    else
    myLabel{i} = num2str(a-1+i);
    end
    
end

figure(fig+col), %subplot(2,3,col)
cg = circularGraph(xu,'Label',myLabel); 

% close all

end 


end






%% SAVE FIGs
% 
%  fig_1 = ['fig_2ba_',num2str(tim)];
% 
% g1 = figure(fig);  
% set(g1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% set(g1,'Units','Inches');
% pos = get(g1,'Position');
% set(g1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(g1,fig_1,'-dpdf','-r0')
% 
% for iu = 1:5
%    
%      fig_2 = ['fig_2bb_',num2str(iu)];
% 
% g2 = figure(fig+iu);  
% set(g2, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% set(g2,'Units','Inches');
% pos = get(g2,'Position');
% set(g2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(g2,fig_2,'-dpdf','-r0')
%     
% end
