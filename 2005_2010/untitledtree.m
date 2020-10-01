
clc
clear all
close all

load wheel_2010

[X,Y] = find(triu(xu) > 0);
coordinate = horzcat(X,Y);
x  = vertcat(X,Y);
detect_ij = unique(coordinate)';

[a,b]=sort(hist(x,unique(x)),'descend');
sort_detect_ij = detect_ij(b)
disp(a)

disp('HOT-SPOT')
disp(sort_detect_ij(1))

 [Headers, Sequences] = multialignread('full_align_na_2010.aln');

%%

 dv = 0.05;
 
f_i = mean(K, 1);
f_n = mean(K, 2); 
f_mut = sum(K, 2);

subplot(1,3,1), hist(f_n,50); axis square, hold on, title('dist. of genomes tot. mut. fq. (MFq)')

uu = find(f_n >= dv); disp(size(uu,1))
u0 = find(f_n < dv);  disp(size(u0,1))

r = randi([1 size(u0,1)],1,size(uu,1));

K_L = K(r,:); %disp(size(K_L))
K_R = K(uu,:); %disp(size(K_R))
K_tot = vertcat(K_L,K_R);

tot_ind = vertcat(uu,r');


subplot(1,3,2), hist(f_n(tot_ind'),50); axis square, hold on, title('dist. of genomes after D.W. genomes with MFq < dv (dv=0.05)')

%% plot only 50 sequences per group:

r1 = randi([1 size(uu,1)],1,50);

 K_L_ij = K_L(r1,sort_detect_ij);
 K_R_ij = K_R(r1,sort_detect_ij);
 K_tot_ij = vertcat(K_R_ij,K_L_ij);
 
subplot(1,3,3),
 imagesc(K_tot_ij), set(gca,'FontSize', 9,'xTickLabels',sort_detect_ij), %xticks(1:length(sort_detect_ij)), 
 xticks(1:16)
 xlabel('detected epistatic sites'), ylabel('LEFT peak      -       RIGHT peak'), axis square
 
 %%
S1 = Sequences(uu(1:50)); 
S2 = Sequences(r(1:50)); 
S = vertcat(S1,S2);

% names = Headers(r);

dist = seqpdist(S,'ScoringMatrix','GONNET','Alphabet','AA');
tree = seqlinkage(dist);%,'single',seqs)
% subplot(2,2,4), 
plot(tree)
% showalignment(ma)


