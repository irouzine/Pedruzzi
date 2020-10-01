clc 
close all
clear all %#ok<CLALL>
clearvars 

% format shortg

% Parameters
ii = 0;
fs = 22;
                                           L = 40;
                                             
                                                     th = .6;

                                           runn = 200;  
                                              N = 1000;  N1 = N; 
                                              
                                              
                                              for tp = 2 %1:4
                                              
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
    stra = 'simple arches';
    topology = 'arches'; st = 2; % 1-2
    
elseif tp == 2
    stra = 'double arches';
    topology = 'con_arches';                                st = 3; % 1-2-3
    
elseif tp == 3
    stra = 'triple arches';
    topology = 'con_arches+';                               st = 10+3; % 1-2-3, 1-3 %also detects asymettric?? what about RAND
    
elseif tp == 4
    stra = '1-2-3-4';
    topology = '1-2-3-4'; st = 4;
    
    
elseif tp == 5
    stra = '1-2-3-4-n-n+1';
    topology = '1-2-3-4-n+1'; st = 2;
else
end


 
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
      
      case '1-2-3-4-n+1' 
        
     EE(i,i-1) = 1;  EE(i,i+1) = 1; 
      EE(i-1,i) = 1;  EE(i+1,i) = 1; 
      
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
                                              
str = sprintf(' Double-arches topology \n Num. populations = %g \n Th = %g \n E = %g \n muL = %g \n s_0 = %g \n f_0 = %g \n N = %d \n L = %g \n t = %g  ', runn, th, E_st, U, s0, f0/L, N, L, td); %disp('a=0, pos epi.'),

%
 
 
 % Epistatic analysis: UFE_{AB} and UFE_{AB0}
 
   for AVG = [1 runn] % @ POPs max
       
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
  UFE(UFE >=1.5) = 0; %Topology-dependent value
                                                                 
                                                                                                                                                  
%%

 

                                             

                                              %%
                                              
                                               tic % Cluster inspection:
                                              
                                              
stat = horzcat(ip,UFE); m = 1; mm = m; tt = 1; zz = length(I); n = 7;  
pick = ceil(rand(1,200)*np); Pi = sort(pick); stat_r = stat(Pi, :); 
 
      [I,J] = find(Ee);     
  % --------------------------------------------      
  % Threeshold

   tUF = th;  % prctile(stat_r(:,7),80); 
 
  % ----------------------------------------
  
[N2] = find(stat(:,3) > tUF); raw_det(t).Pairs = N2; % memorize all detected pairs "RAW"
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

Links = Raw(:,3); %Links(Links>1) = NaN; 



G = graph(pairs(:,1),pairs(:,2),Links,chr); 
                                           

clc
format short

Bin = conncomp(G);

[bin,~,binsize] = histcounts(Bin);

Cls = bin(bin > 2); disp(Cls)

% note: does not separate cluster of same size

sb = length(unique(Cls));
sup = 0;

if sb == 0
    disp('only isolated pairs are detected')
else
end

for cl = 1%:sb
    
   clear UFEiii fhap Tt Ff T F
    
idx = bin(binsize) == Cls(cl);

ss = [1:L];
ss = ss(idx);

Gsub = subgraph(G, ss); hold on
 
if AVG == 1
figure(1), sub = plot(Gsub,'NodeColor','k','EdgeColor','k','EdgeAlpha',0.4, 'NodeLabel',[], 'EdgeLabel',[],'LineWidth',2); axis square, hold on, % separate each subclusters (greater than two sites)
xlabel('Raw UFE_{AB}')
set(gca,'visible','off')

set(gca, 'FontSize', fs)
xticks([]), xticklabels('')
yticks([]), yticklabels('')


else
    

 
 lbs = round(Gsub.Edges.Weight,2);   
 format bank
    
 
 
 
    
figure(2), sub = plot(Gsub,'NodeColor',[.1 .1 .1],'EdgeColor','k','EdgeAlpha',0.5, 'EdgeFontSize', 11, 'NodeLabel',[], 'EdgeLabel',[], 'LineWidth',3); axis square, hold on, % separate each subclusters (greater than two sites)
xlabel('Raw UFE_{AB}')
set(gca,'visible','off', 'FontSize', fs)
xticks([]), xticklabels('')
yticks([]), yticklabels('')  
   
%%

% tl = sub.EdgeLabel;
% 
%  itvec = 1:length(tl);
%  Newlbs = zeros(size(itvec));
% %  
% for i= itvec
%     textstr = string(tl(i));
%     Newlbs(i) = round(str2double(textstr), 2);
% end
% 
% sub = plot(Gsub,'NodeColor',[.1 .1 .1],'EdgeColor','k','EdgeAlpha',0.5, 'EdgeFontSize', 11, 'NodeLabel',[], 'EdgeLabel',Newlbs', 'LineWidth',2);

%%

end


if AVG > 1

% identify the sites within sub-clusters > 2 sites
Sub_Cl = [1:L]; Sub_Cl = Sub_Cl(idx);

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
%    clc
   
UFEiii_min = zeros(test,1);


for it = 1 : test
    
ii = ismember(sub_Raw(:,[1 ]),sub_Raw(it,[1 ]),'rows')'; %first site in first column
ij = ismember(sub_Raw(:,[1 ]),sub_Raw(it,[2 ]),'rows')'; %second site in first column
Ii = or(ii,ij); % rows of Sub_raw (1st col) that contains nodes i or j

ji = ismember(sub_Raw(:,[2 ]),sub_Raw(it,[1 ]),'rows')'; %first site in second column
jj = ismember(sub_Raw(:,[2 ]),sub_Raw(it,[2 ]),'rows')'; %second site in second column
Jj = or(ji,jj); % rows of Sub_raw (2st col) that contains nodes i or j

IJ = or(Ii, Jj); % set rows containing i and j nodes.

i_j = sub_Raw(it,[1 2]);
nodes = sub_Raw(IJ,[1 2]); Nodes = unique(reshape(nodes,1,[]));

te = setdiff(Nodes(:,:),i_j);
  test2 = size(te,2);   
    
    UFEiii = zeros(1,test2);
    
    for i0 = 1:test2  
     
        a = sub_Raw(it,1); b = sub_Raw(it,2); c = te(i0); 
        
%  sites = [a b c]     
    %% 
    
%       a = 13; b = 14; c = 15; %to check haplotype for specif triplet

f00 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 0 0],2));
f01 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 1 0],2));
f10 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 0 0],2));
f11 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 1 0],2));

 fhap(it,:) = [f00 f01 f10 f11];   
 %%
UFEiii(1,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); disp(UFEiii(1,i0));
%  pause

disp(UFEiii(1,i0))

    end

     UFEiii(UFEiii <0) = 0; 
%   UFEiii(UFEiii>=1.5) = 0; 
  

% UFEiii(1,:)


 UFEiii_min(it,1)  = min([UFEiii(1,:), 777],[],2); % cause some pairs are "terminal"

%         UFEiii_min(it,1)  = min(UFEiii(1,:),[],2); 

end
  

prufe = (sub_Raw(:,[3])./UFEiii_min);

prufe(isinf(prufe)) = 10; %arbitrarily set to 10 fold reduction (because UFE_{AB0} = 0)
Results = horzcat(sub_Raw(:,[1 2 3]),UFEiii_min,prufe);

disp('    A         B         UFE_AB   min(UFE_AB0)  ratio ')

disp(Results) 

                                                                                    figure(3), 
sca = scatter(Results(:,3),Results(:,4),200,'MarkerEdgeColor','w','MarkerFaceColor','g','LineWidth', 2); axis square, xlabel('UFE_{ij}'), ylabel('UFE_{ij0}'), hold on,

Res = Results(:,[3:5]);

idx = find(Res(:,3) < 2);
Res(idx,:) = [];

scatter(Res(:,1),Res(:,2),200,'MarkerEdgeColor','w','MarkerFaceColor','r','LineWidth', 2), hold on

set(gca,'FontSize', fs,'linewidth',3)
a = annotation('textbox', [0.35, 0.75, 0.1, 0.1], 'String', str, 'EdgeColor', 'w');%, 'Interpreter', 'latex');
a.FontSize = 20;
dim = [.2 .74 .25 .15];
xlim([0 1.5]), ylim([0 1.5])

% b = annotation('rectangle',dim,'Color','green', 'LineWidth', 2);

x = 0:0.01:1.5;
y1 = get(gca,'ylim');

cc = [.7 .7 .7];

hold on
plot([th th],y1,'Color',cc, 'LineWidth', 3), hold on

line(x,.5*x,'Color',cc,'LineStyle','--', 'LineWidth', 3), hold on
line(x,x,'Color',cc,'LineStyle','--', 'LineWidth', 3), hold on



% threeshold must be otpimised:
                                                                            TRUE = abs(prufe) < 2; % UFE_{AB}/UFE_{AB0} > 2

GsuB = subgraph(G, ss);
GsuB.Edges.Weight = UFEiii_min;
Lbs = round(GsuB.Edges.Weight,2);  

% Lbs(Lbs == 0) = [];

  cdarkred = [.54 0 0];  
  
  figure(4) 
h = plot(GsuB,'NodeColor',[.5 .5 .5],'EdgeColor','r','EdgeAlpha',0.6, 'EdgeFontSize', 14, 'NodeLabel',[], 'EdgeLabelColor',cdarkred, 'EdgeLabel',Lbs,'LineWidth',2); xlabel('min(UFE_{AB0})'), axis square, hold on % separate each subclusters (greater than two sites)
set(gca, 'FontSize', 16)
xticks([]), xticklabels('')
yticks([]), yticklabels('')
set(gca,'visible','off')

%%

% tl = sub.EdgeLabel;
% 
%  itvec = 1:length(tl);
%  Newlbs = zeros(size(itvec));
% %  
% for i= itvec
%     textstr = string(tl(i));
%     Newlbs(i) = round(str2double(textstr), 2);
% end
% 
% h = plot(GsuB,'NodeColor',[.5 .5 .5],'EdgeColor','r','EdgeAlpha',0.6, 'EdgeFontSize', 11, 'NodeLabel',[], 'EdgeLabel',Newlbs','LineWidth',2);
% %%

tag(TRUE) = {'true'}; 
tag(~TRUE) = {'false'};

sup = sup +2;
%

T = TRUE;
F = ~TRUE;

T = sub_Raw(T, [1 2]); 
F = sub_Raw(F, [1 2]); 

% This section recover the position of "ith" Node, to reorder the Nodes
% for TRUE/FALSE coloring

if size(T,1) == 0  &&  size(F,1) == 0

disp('no TRUE pairs detected')

else
    
Tt = T(1,:); %Ff = F(1,:);

for sp = 2 : size(T,1)  
    Tt = horzcat(Tt,T(sp,:));
end
% for sp = 2 : size(F,1)  
%     Ff = horzcat(Ff,F(sp,:));
% end

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
cdarkgreen = [0 100/256 0];  
% highlight(h,Tt_0)
highlight(h,Tt_0, 'EdgeColor','g','LineWidth',4, 'NodeColor', 'k', 'EdgeLabelColor',cdarkgreen)

end

end    

toc




end

   end


%%
% 
% % % % SAVE FIGs 
% g1 = figure(1); g2 = figure(2); g3 = figure(3); g4 = figure(4);
% 
% % close(g1), close(g2), close(g4)
%  
% 
% fig_n1 = 'fig_2a.pdf'; fig_n2 = 'fig_2b.pdf'; fig_n3 = 'fig_2c.pdf'; fig_n4 = 'fig_2d.pdf';
% 
% % set(g1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'Units','Inches'); pos = get(g1,'Position');
% % set(g1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % print(g1,fig_n1,'-dpdf','-r0')
% 
% % set(g2, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'Units','Inches'); pos = get(g2,'Position');
% % set(g2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % print(g2,fig_n2,'-dpdf','-r0')
% 
% set(g3, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'Units','Inches'); pos = get(g3,'Position');
% set(g3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(g3,fig_n3,'-dpdf','-r0')
% % 
% % set(g4, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'Units','Inches'); pos = get(g4,'Position');
% % set(g4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % print(g4,fig_n4,'-dpdf','-r0')










