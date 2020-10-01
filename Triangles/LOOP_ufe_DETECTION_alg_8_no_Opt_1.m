clc
close all
clear all %#ok<CLALL>
clearvars 


tic

% Parameters
                                              L = 10;  
                                              N = 2000;  N1 = N;
                                              
  s0 = .1; 
E_st = .75;
  f0 = .45*L;                                                
   U = 0.07; 
  td = 49+1;  
 
        r = 0; 
      M = 0;  
                                                                      
muL = U;  
mu=muL/L; 

c_1 = N * (mu/s0)^2;
SUlogNs = s0/(mu*L)*log(N*s0);

EE = zeros(L,L);
det_mat = zeros(L,L,td);


runn = 20;
 

%create a random matrix:
% 
% topology = 'arches'; st = 2; % 1-2
% topology = 'con_arches'; st = 3; % 1-2-3
topology = 'con_arches+'; st = 3; % 1-2-3, 1-3
% topology = '1-2-3-4'; st = 4;


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
                                                                
%  [DAt, ww] = epi('binary',  s0,  0, L, N,  td,  2,   U,  f0, 1 , E_st, Ee);   
% 
%  
%  
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
                                 
 
   for AVG = [runn]     
       
         tmin = 1;
         tmax = td;
         
                                                                           for t = tmin:tmax % time of analysis

ss = [1:L];
I = find(Ee);                                                                        
     
%DATA = DAt(:,:,t); 
 DATA = Dat_run(1:AVG*N,:,t);
 

% add detection code    
    
    Ko = DATA; Ng = N*AVG; 
      
i_div = (1:L); ip = nchoosek(i_div,2); np = size(ip,1);
fhap = zeros(np,4); UFE = zeros(np,1); %WU = zeros(np,1);
UFEnot = zeros(np, L-2);

for i=1:np
    
f00 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 0],2)); f01 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 1],2));
f10 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 0],2)); f11 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 1],2));
                 
fhap(i,:) = [f00 f01 f10 f11]; UFE(i) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); WU(i) = log((f11*f00)/(f01*f10));

te = setdiff(i_div,ip(i,:), 'stable'); % A that is not in B
    test2 = size(te,2);
     
    for i0 = 1:test2
   
        a = ip(i,1); b = ip(i,2); c = te(i0);       
        
f00 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 0 0],2)); 
f01 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[0 1 0],2));
f10 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 0 0],2)); 
f11 = mean(all(Ko(:,[a,b,c]) == ones(Ng,1)*[1 1 0],2));
               
UFEiii(i,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));  %#ok<SAGROW> % UFE for triplets

    end


end  
    
UFE(UFE<0) = 0; 
UFE(UFE>=1) = 0; 

UFEiii(UFEiii<0) = 0; 
UFEiii(UFEiii>=1) = 0; 

     
     stat0 = mean(UFEiii,2);
     
statp = UFE.*stat0;

stat_0 = horzcat(ip,UFE, stat0, statp);
 
            %%

[I,J] = find(Ee); 
     
 zz = length(I); 

  %% --------------------------------------------      

  % Threeshold

  [~,u2,u3] = kmeans(UFE,2)
  tUF =  .9*max(u2)%prctile(stat_r(:,3),80);%mean(UFE) + std(UFE);
  
  [~,v2,v3] = kmeans(stat0,2)
  tUF1 = .9*max(v2)   % <nonUFE>
  
  [~,x2,x3] = kmeans(statp,2)
  tUF2 = .9*max(x2)   % UFE * <nonUFE>
     
   tUF_t(t) = tUF;
   %%
  % ----------------------------------------
  
  
[N2] = find(stat_0(:,3) > tUF); [N3] = find(stat_0(:,4) > tUF1); [N4] = find(stat_0(:,5) > tUF2);





% N2 = C1;





% tUF
det = 0;FP = 0;  c4 = isempty(N2);
 
     if c4 == 0

for n2 = 1 : size(N2,1) 
             if EE(stat_0(N2(n2),1), stat_0(N2(n2),2)) == 1
                det_mat(stat_0(N2(n2),1), stat_0(N2(n2),2),t) = -1; % detected
                det = det+1;
             else
                 FP = FP+1;
                 det_mat(stat_0(N2(n2),1), stat_0(N2(n2),2),t) = 1; % false positive
             end    
end

         det = 100*det/zz; FP = 100*FP/n2;    
         out = [det;FP];
         
    else
              
      out = [det;FP];
      
     end
   
    
     % tUFnot
det = 0;FP = 0;  c5 = isempty(N3);
 
     if c5 == 0

for n3 = 1 : size(N3,1) 
             if EE(stat_0(N3(n3),1), stat_0(N3(n3),2)) == 1
                det_mat(stat_0(N3(n3),1), stat_0(N3(n3),2),t) = -1; % detected
                det = det+1;
             else
                 FP = FP+1;
                 det_mat(stat_0(N3(n3),1), stat_0(N3(n3),2),t) = 1; % false positive
             end    
end

         det = 100*det/zz; FP = 100*FP/n3;    
         out1 = [det;FP];
         
    else
              
      out1 = [det;FP];
      
     end
    
       % tUF * <tUFnot>
       
det = 0;FP = 0;  c6 = isempty(N4);
 
     if c6 == 0

for n4 = 1 : size(N4,1) 
             if EE(stat_0(N4(n4),1), stat_0(N4(n4),2)) == 1
                det_mat(stat_0(N4(n4),1), stat_0(N4(n4),2),t) = -1; % detected
                det = det+1;
             else
                 FP = FP+1;
                 det_mat(stat_0(N4(n4),1), stat_0(N4(n4),2),t) = 1; % false positive
             end    
end

         det = 100*det/zz; FP = 100*FP/n4;    
         out2 = [det;FP];
         
    else
              
      out2 = [det;FP];
      
     end  
     
     
    % -------------------------
         
for ji = 1 : zz
    
   det_mat(I(ji), J(ji),t) = -1; % known pairs
    
end
 
  % -------------------------
 
% x PLOT

Y = horzcat(out, out1, out2);
                    
     OUT_det_UF(1,t) = Y(1,1); 
     OUT_fp_UF(1,t) = Y(2,1);  

     OUT_det_UFnot(1,t) = Y(1,2); 
     OUT_fp_UFnot(1,t) = Y(2,2);  

     OUT_det_UFp(1,t) = Y(1,3); 
     OUT_fp_UFp(1,t) = Y(2,3);  

     
   %%  

                                                                           end


%%

     Det_UF = mean(OUT_det_UF,3); 
     Fp_UF = mean(OUT_fp_UF,3); 
      
      Det_UFnot = mean(OUT_det_UFnot,3); 
     Fp_UFnot = mean(OUT_fp_UFnot,3); 
     
      Det_UFp = mean(OUT_det_UFp,3); 
     Fp_UFp = mean(OUT_fp_UFp,3); 
%% ----------------------------

figure(1)
subplot(3,2,1), plot(Det_UF','LineWidth',2), title('Detection UFE'),  hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax]), title('UFE')
subplot(3,2,2), plot(Fp_UF','LineWidth',2),  title('False Pos. UFE'), hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax])

subplot(3,2,3), plot(Det_UFnot','LineWidth',2), title('Detection UFE'),  hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax]), title('UFEnot')
subplot(3,2,4), plot(Fp_UFnot','LineWidth',2),  title('False Pos. UFE'), hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax])

subplot(3,2,5), plot(Det_UFp','LineWidth',2), title('Detection UFE'),  hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax]), title('UFEp')
subplot(3,2,6), plot(Fp_UFp','LineWidth',2),  title('False Pos. UFE'), hold on, axis square, set(gca,'FontSize',16), ylim([0 100]), xlim([tmin tmax])


    legend off
    %%

   end


% %----------------------------
% gm = figure(1);
% 
% fig_n2 = 'FIGURE';
% set(gm, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% 
% set(gm,'Units','Inches');pos = get(gm,'Position');
% set(gm,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gm,fig_n2,'-dpdf','-r0')
% 
% 

%%
% 
% figure(2)
% 
% Z = peaks;
% surf(Z); 
% axis tight manual 
% set(gca,'nextplot','replacechildren'); 
% 
% 
% v = VideoWriter('DET.avi', 'Uncompressed AVI');
% open(v);
% 
% E_k = -EE;
% 
% for k = tmin:tmax
% 
%   a =  surf(det_mat(:,:,k),'LineStyle','none'); colormap(redblue), view([0 90]), axis square, 
%   ylabel('KNOWN Epis. pairs'), xlabel('DETECTED Epis. pairs'), colorbar, caxis ([-1 1]), 
%   
%    frame = getframe(gcf);
%    writeVideo(v,frame);
%    
% end
% 
% close(v);
% 
% %%
%  implay('DET.avi',1);
% 
% 
% toc
% 
% 
% % disp(Thress)
% % pause
% 
% 
% 
% %%
% %%
% %%
% %%
% 
% 
%  %%    
% % Analysis of False-Bonds by triplets f_{110}
% 
% tmin = 4;
% 
% Nn =  raw_det(tmin).Pairs;
% 
% Raw = stat(Nn,:);
% 
% % all sites involved in pairs:
% 
% sites = unique(Raw(:,1:2)); %size(sites)
% pairs = Raw(:,1:2);
% 
% % G = graph(pairs(:,1),pairs(:,2)); %   weights: Raw(:,7),L);
% % 
% figure(3)
% % subplot(1,2,1), plot(G,'NodeColor','k','EdgeAlpha',0.9);
% % 
% 
% %
% 
% for k = 1 : L
%     chr{k} = num2str(k);
% end
% 
% chr = chr';
% 
% %
% 
% G1 = graph(pairs(:,1),pairs(:,2),Raw(:,end),chr);
%      
% %subplot(2,2,1),
% plot(G1,'NodeColor','k','EdgeAlpha',0.9,'EdgeLabel',G1.Edges.Weight);
% 
% 
% % Nei = neighbors(G1,'3');
% % H = subgraph(G1, ['3'; Nei]);
% % subplot(1,2,2),
% % plot(H,'EdgeLabel',H.Edges.Weight)
% 
% %%
% clc
% format short
% 
% Bin = conncomp(G1,'Type','weak');
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
% figure(4)
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
% Gsub = subgraph(G1,ss); hold on
% % subplot(2,sb,cl), 
% subplot(1,2,1), sub = plot(Gsub,'NodeColor','k','EdgeAlpha',0.5, 'EdgeLabel',Gsub.Edges.Weight,'LineWidth',2); axis square% separate each subclusters (greater than two sites)
% 
% 
% % identify the sites within sub-clusters > 2 sites
% 
% Sub_Cl = [1:L]; Sub_Cl = Sub_Cl(idx)
% disp(' ')
% disp(' ')
% 
% [~,b1] = intersect(Raw(:,1), Sub_Cl);
% [~,b2] = intersect(Raw(:,2), Sub_Cl);
% 
% %b = mod(bb,size(Raw,1))
% %b(b > size(Raw,1)) = mod(b,size(Raw,1))
% 
% b = union(b1,b2);
% 
% sub_Raw = Raw(b',:); 
% disp('sub_RAW'), disp(sub_Raw(:,[1 2 7])) % i,j and UFE (3-6 fhap)
% 
% Ko = DATA;
% 
% %consider all triplets:
% 
%     nt = size(sub_Raw,1);
% %     UFEiii = zeros(nt,3);
%         
%    test = size(Gsub.Edges,1); % # of bonds to be tested
%     
%     z
% for it = 1 : test
%     
%     
%     te = setdiff(Sub_Cl,sub_Raw(it,[1 2]));
%     test2 = size(te,2);
%      
%     for i0 = 1:test2
%     
% 
%         a = sub_Raw(it,1); b = sub_Raw(it,2); c = te(i0); 
%  
% 
%         %disp(a), disp(b), disp(c)
%         
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
% Results = horzcat(sub_Raw(:,[1 2 7]),UFEiii_avg);
% ratio = (Results(:,3)./Results(:,4))-1;
% 
% disp('    A         B         UFE_AB   <UFE_ABc>  ratios ')
% Results = horzcat(Results, ratio); disp(Results)
% 
% 
% % threeshold must be otpimised, right now, quite a random 0.4 without much meaning
% TRUE = ratio < .09;
% 
% 
% 
% %
% 
% Gsub.Edges.Weight = UFEiii_avg;
% % 
% % subplot(2,sb,sb+cl), 
% subplot(1,2,2), h = plot(Gsub,'NodeColor','k','EdgeAlpha',0.5, 'EdgeLabel',Gsub.Edges.Weight,'LineWidth',2); axis square% separate each subclusters (greater than two sites)
% 
%      
% tag(TRUE) = {'true'}; 
% tag(~TRUE) = {'false'};
% 
% labeledge(h,1:numedges(Gsub),tag);
% 
% end
%%


Zz = horzcat(stat_0, UFEiii);