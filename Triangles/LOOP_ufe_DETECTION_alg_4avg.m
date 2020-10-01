% clc
% close all
% clear all %#ok<CLALL>
% clearvars 
% 
% global stat zz Wei EE
% 
% tic
% 
% 
% 
% % Parameters
%                                               Wei = 1;
%                                               L = 50;  
%                                               N = 5000;  N1 = N;
%                                               
%   s0 = .1; 
% E_st = .75;
%   f0 = .45*L;                                                
%    U = 0.07; 
%   td = 30+1;  
%  
%         r = 0; 
%       M = 0;  
%                                                                       
% muL = U;  
% mu=muL/L; 
% 
% c_1 = N * (mu/s0)^2;
% SUlogNs = s0/(mu*L)*log(N*s0);
% 
% EE = zeros(L,L);
% det_mat = zeros(L,L,td);
% 
%  
% runn = 20;
%  
%  
%  
%  
% %create a random matrix:
% 
% topology = 'arches';
% 
% for i = 2:2:L
%     
% switch topology
%       
%     case 'arches' 
%         
%       EE(i,i-1) = 1; 
%       EE(i-1,i) = 1;
%   
%         
%     case 'rand' 
% 
% a = i+1;
% b = L-a;
% 
% r = ceil(a + (b-a)*rand(1));
% 
%      EE(i,r) = 1; 
%      EE(r,i) = 1;
% 
% end
%     
% end
% 
% Ee = tril(EE,-1);   
%                                                                 
% %  [DAt, ww] = epi('binary',  s0,  0, L, N,  td,  2,   U,  f0, 1 , E_st, Ee);   
% % 
% %  
% %  
% Dat_run = zeros(N*runn,L,td);
%    
% 
% for nn = 1:runn
% 
%  [DAt, ww] = epi('binary',  s0,  0, L, N,  td,  2,   U,  f0, 1 , E_st, Ee);   
% %  [DAt, ww] = epi('half',  s0,  0, L, N,  td,  2,   U,  f0, 1 , E_st, Ee, coef);   
% 
% a = 1 + (nn-1)*N;
% b = N + (nn-1)*N;
% 
% disp(nn)
% 
% Dat_run(a:b,:,:) = DAt(:,:,:);
% 
% end
% 
%  %%
%                                  
%  
%    for AVG = [1 runn]     
%        
%          tmin = 1;
%          tmax = td;
%          
%                                                                            for t = tmin:tmax % time of analysis
% 
% 
% I = find(Ee);                                                                        
%      
% %DATA = DAt(:,:,t); 
%  DATA = Dat_run(1:AVG*N,:,t);
%  
% 
% % add detection code    
%     
%     Ko = DATA; Ng = N*AVG; 
%       
% i_div = (1:L); ip = nchoosek(i_div,2); np = size(ip,1);
% fhap = zeros(np,4); UFE = zeros(np,1); WU = zeros(np,1);
% 
% for i=1:np
% f00 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 0],2)); f01 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[0 1],2));
% f10 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 0],2)); f11 = mean(all(Ko(:,[ip(i,1), ip(i,2)]) == ones(Ng,1)*[1 1],2));
%                  
% fhap(i,:) = [f00 f01 f10 f11]; UFE(i) = 1 - (log(f11/f00))/((log(f01*f10/f00^2))); WU(i) = log((f11*f00)/(f01*f10));
% end  
%     
%    
% 
% UFE(UFE<0) = 0; 
% UFE(UFE>1) = 0; 
% 
% stat = horzcat(ip,fhap,UFE);
%      
%             %%
%  
% m = 1; mm = m; tt = 1; 
%  zz = length(I); n = 7;  
%  
% pick = ceil(rand(1,200)*np);
% Pi = sort(pick);
% stat_r = stat(Pi, :); 
%  
%       [I,J] = find(Ee); 
%       
%   % --------------------------------------------      
% 
%   % Threeshold
% 
%    tUF = prctile(stat_r(:,7),80); 
%    
%   % ----------------------------------------
%   
%   v = tUF; [tUF_opt] = fminsearch(@tUFE,v);
%   
% [N2] = find(stat(:,7) > tUF_opt);
% 
%                                                         raw_det(t).Pairs = N2; %#ok<SAGROW> % memorize all detected pairs "RAW"
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
% 
%      
% 
% 
% 
% 
% 
%    %%  
% 
%                                                                            end
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
%     legend off
% 
%    end
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
% % 
% 
% %%
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
% implay('DET.avi',1);
% 
% 
% toc







%%
%%
%%
%%


 %%    
% Analysis of False-Bonds by triplets f_{110}

tmin = 19;

Nn =  raw_det(tmin).Pairs;

Raw = stat(Nn,:);

% all sites involved in pairs:

sites = unique(Raw(:,1:2)); %size(sites)
pairs = Raw(:,1:2);

% G = graph(pairs(:,1),pairs(:,2)); %   weights: Raw(:,7),L);
% 
figure(3)
% subplot(1,2,1), plot(G,'NodeColor','k','EdgeAlpha',0.9);
% 

%

for k = 1 : L
    chr{k} = num2str(k);
end

chr = chr';

%

G1 = graph(pairs(:,1),pairs(:,2),Raw(:,end),chr);
     
%subplot(2,2,1),
plot(G1,'NodeColor','k','EdgeAlpha',0.9,'EdgeLabel',G1.Edges.Weight);


% Nei = neighbors(G1,'3');
% H = subgraph(G1, ['3'; Nei]);
% subplot(1,2,2),
% plot(H,'EdgeLabel',H.Edges.Weight)

%%
clc
format short

bin = conncomp(G1,'Type','weak');

%[za,zb] = hist(bin,unique(bin))  % unique(bin) number of cluster (including single elements)
[bin,~,binsize] = histcounts(bin);

Cls = bin(bin > 2); disp(Cls)

figure(4)

for cl = 1:length(Cls)
    
idx = bin(binsize) == Cls(cl);
Gsub = subgraph(G1, idx); hold on
subplot(2,2,cl), sub = plot(Gsub,'NodeColor','k','EdgeAlpha',0.5, 'EdgeLabel',Gsub.Edges.Weight,'LineWidth',2); % separate each subclusters (greater than two sites)


% identify the sites within sub-clusters > 2 sites

Sub_Cl = [1:L]; Sub_Cl = Sub_Cl(idx)

[~,b1] = intersect(Raw(:,1), Sub_Cl);
[~,b2] = intersect(Raw(:,2), Sub_Cl);

%b = mod(bb,size(Raw,1))
%b(b > size(Raw,1)) = mod(b,size(Raw,1))

b = union(b1,b2);

sub_Raw = Raw(b',:); disp(sub_Raw(:,[1 2 7])) % i,j and UFE (3-6 fhap)

Ko = DATA;

%consider all triplets:

ii_div = Sub_Cl;
iip = nchoosek(ii_div,3)
nt = size(iip,1);

fhap_iii = zeros(nt,4);
    UFEiii = zeros(nt,3);
    UFEii = zeros(nt,1); ratio = UFEiii;
    
    i_div = (1:L); ip = nchoosek(i_div,2); np = size(ip,1);
    
    
for i=1:nt
    
    for i0 = 1:3
    
        if i0 == 1
        a = 1; b = 2; c = 3; % AB-c
        elseif i0 == 2
        a = 1; b = 3; c = 2; % AC-b
        elseif i0 == 3
        a = 2; b = 3; c = 1; % BC-a
        else
        end
        
 
f00 = mean(all(Ko(:,[iip(i,a), iip(i,b)]) == ones(Ng,1)*[0 0],2)); 
f01 = mean(all(Ko(:,[iip(i,a), iip(i,b)]) == ones(Ng,1)*[0 1],2));
f10 = mean(all(Ko(:,[iip(i,a), iip(i,b)]) == ones(Ng,1)*[1 0],2)); 
f11 = mean(all(Ko(:,[iip(i,a), iip(i,b)]) == ones(Ng,1)*[1 1],2));
               
UFEii(i,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));  % UFE for triplets 



f00 = mean(all(Ko(:,[iip(i,a), iip(i,b), iip(i,c)]) == ones(Ng,1)*[0 0 0],2)); 
f01 = mean(all(Ko(:,[iip(i,a), iip(i,b), iip(i,c)]) == ones(Ng,1)*[0 1 0],2));
f10 = mean(all(Ko(:,[iip(i,a), iip(i,b), iip(i,c)]) == ones(Ng,1)*[1 0 0],2)); 
f11 = mean(all(Ko(:,[iip(i,a), iip(i,b), iip(i,c)]) == ones(Ng,1)*[1 1 0],2));
               
% fhap(i,:) = [f00 f01 f10 f11]; 
UFEiii(i,i0) = 1 - (log(f11/f00))/((log(f01*f10/f00^2)));  % UFE for triplets

ratio(i,i0) = UFEii(i,i0)/UFEiii(i,i0)-1; 


    end

end
   
% 
% UFE(UFE<0) = 0; 
% UFE(UFE>1) = 0; 
% 
% stat = horzcat(ip,fhap,UFE);
   disp('           ') 
disp('   A         B         C          /   UFE_AB  UFE_ABc  ratio    /  UFE_AC  UFE_ACb  ratio     /  UFE_BC  UFE_BCa  ratio ')
Results = horzcat(iip,UFEii(:,1),UFEiii(:,1),ratio(:,1),UFEii(:,2),UFEiii(:,2),ratio(:,2),UFEii(:,3),UFEiii(:,3),ratio(:,3))






end

%%
subplot(2,2,2+cl), h = plot(Gsub,'NodeColor','k','EdgeAlpha',0.5, 'EdgeLabel',Gsub.Edges.Weight,'LineWidth',2); % separate each subclusters (greater than two sites)
% highlight(sub_h,[1 3])


highlight(h,3,'NodeColor',[0 0.75 0])
highlight(h,n3,'NodeColor','red')
highlight(h,3,n3,'EdgeColor','red')












