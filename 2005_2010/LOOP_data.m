clc
clear all %#ok<CLALL>

% Read DATA
tic 
seq = fastaread('804519368671-2005_2010.fasta'); % NA 2005-2010

Seq = seq; iu = 0;

for ii = 1 : length(seq) 
    
   disp(ii)
    seq(ii).aa_length = size(seq(ii).Sequence,2);
    seq(ii).ambiguous = strfind(seq(ii).Sequence,'J');
    
     
     if isempty(strfind(seq(ii).Sequence,'J'))    
     else
         Seq(ii-iu) = []   ;
         iu = iu+1;        % correct row index
     end
     
end  

%%
% 'Seq' dos NOT include sequences with ambiguous carachter (J)
S = Seq(1:end); % entire set.

for ih = 1 : length(S) 
        
    S(ih).Header = num2str(ih);
 
end 

tic
dist = seqpdist(S,'ScoringMatrix','GONNET','Alphabet','AA');
tree = seqlinkage(dist,'average',S);
ma = multialign(S,tree,'ScoringMatrix', {'pam150','pam200','pam250'},'UseParallel',true);
% showalignment(ma)
multialignwrite('full_align_na_2010.aln', ma)
toc
%%



%
                    [Headers, Sequences] = multialignread('full_align_na_2010.aln');
%

for iy = 1 : length(Sequences) 
    disp(iy)
%     Seq(iy).Header = Headers(iy); 
    Seq(iy).sequence_msa = char(Sequences(iy));
    Seq(iy).length = size(Seq(iy).sequence_msa,2);
end  

%
[CSeq,~] = seqconsensus(Sequences,'Gaps', 'all'); REF1 = CSeq; %disp(CSeq)

pos = length(CSeq); 
S = Seq(1:end); Reads = size(S,1);
K = zeros(Reads(1),pos);

% Convert to Binary sequences (0/1)

for i = 1 : size(S,1)
    
    disp(i)
    SEQ = S(i).sequence_msa;
  
    co = strcmp(REF1,SEQ);
    
    if co == 1
   
    K(i,:) = 0;
 
    else
        
   for j= 1:pos
        
      if size(setdiff(SEQ(j),REF1(j)),2) == 1         
        
        K(i,j) = 1;
      else
      end
            
    
   end
        
    end    
end


save('Data_na_2010.mat', 'K', 'pos', 'Seq', 'tree')


%% end of data
