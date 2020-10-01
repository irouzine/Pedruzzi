tic
 clc 
global SZ
 
%                                                                                SZ = 1000;

for tim = 5%:6% 1:2

%     clc
    clearvars -except tim CSeq REF1 SZ

    if tim == 1    
    seq = fastaread('NA_05_07.fasta');  tgseq = 470;% NA 2005-2007
    elseif tim == 2       
    seq = fastaread('NA_08_10.fasta');  tgseq = 469;% NA 2008-2010
    elseif tim == 3
    seq = fastaread('HA_05_07.fasta');  tgseq = 565;% HA 2005-2007
    elseif tim == 4
    seq = fastaread('HA_08_10.fasta');  tgseq = 566;% HA 2008-2010
     
    
    elseif tim == 5
    seq = fastaread('NA_05_08.fasta');  tgseq = 470;% NA 2005-2008
    elseif tim == 6
    seq = fastaread('NA_09_10.fasta');  tgseq = 470;% NA 2009-2010
    end
    
   
% Read DATA



Seq = seq; iu = 0;

for ii = 1 : length(seq) 
    
    disp(ii)
    seq(ii).aa_length = size(seq(ii).Sequence,2);
    seq(ii).ambiguous = strfind(seq(ii).Sequence,'J');
    
     
     if isempty(strfind(seq(ii).Sequence,'J'))  % &&  size(seq(ii).Sequence,2) == tgseq
     else
         Seq(ii-iu) = [];
         iu = iu+1;        % correct row index
     end
     
end  

%%
% 'Seq' dos NOT include sequences with ambiguous carachter (J)
                                                                                S = Seq(1:end); %entire set. 
%                                                                                 S = Seq(1:SZ); 

for ih = 1 : length(S) 
        
    S(ih).Header = num2str(ih);
 
end 

tic
dist = seqpdist(S,'ScoringMatrix','GONNET','Alphabet','AA');
tree = seqlinkage(dist,'average',S);
ma = multialign(S,tree,'ScoringMatrix', {'pam150','pam200','pam250'},'UseParallel',true);
% showalignment(ma)



    if tim == 1  
    multialignwrite('full_align_NA_05_07.aln', ma)
        [Headers, Sequences] = multialignread('full_align_NA_05_07.aln');
     
    elseif tim == 2      
    multialignwrite('full_align_NA_08_10.aln', ma)
        [Headers, Sequences] = multialignread('full_align_NA_08_10.aln');
     
    elseif tim == 3
    multialignwrite('full_align_HA_05_07.aln', ma)
        [Headers, Sequences] = multialignread('full_align_HA_05_07.aln');

    elseif tim == 4
    multialignwrite('full_align_HA_08_10.aln', ma)
        [Headers, Sequences] = multialignread('full_align_HA_08_10.aln');
    
        elseif tim == 5
    multialignwrite('full_align_NA_05_08.aln', ma)
        [Headers, Sequences] = multialignread('full_align_NA_05_08.aln');
              
        elseif tim == 6
    multialignwrite('full_align_NA_09_10.aln', ma)
        [Headers, Sequences] = multialignread('full_align_HA_09_10.aln');
        
        
    end


for iy = 1 : length(Sequences) 
    disp(iy)
%     Seq(iy).Header = Headers(iy); 
    Seq(iy).sequence_msa = char(Sequences(iy));
    Seq(iy).length = size(Seq(iy).sequence_msa,2);
end  

%
%      if tim == 1 || tim == 3
[CSeq,~] = seqconsensus(Sequences,'Gaps', 'all'); REF1 = CSeq; %disp(CSeq)
% else;end

pos = length(CSeq);

disp(tim)
 disp(pos)
 disp(iu)
  disp('--------')
 % pause


S = Seq(1:end); Reads = size(S,1);
K = zeros(Reads(1),pos);

% Convert to Binary sequences (0/1)

                                                                          for i = 1 : size(S,1) %%%%SZ
%                                                                           for i = 1 : SZ
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


    if tim == 1  
    save('Data_NA_05_07.mat', 'K', 'pos', 'Seq', 'tree')  % NA 4 time-points data
    elseif tim == 2
    save('Data_NA_08_10.mat', 'K', 'pos', 'Seq', 'tree')
    elseif tim == 3
    save('Data_HA_05_07.mat', 'K', 'pos', 'Seq', 'tree') 
    elseif tim == 4
    save('Data_HA_08_10.mat', 'K', 'pos', 'Seq', 'tree')
    
     elseif tim == 5
    save('Data_NA_05_08.mat', 'K', 'pos', 'Seq', 'tree')
     elseif tim == 6
    save('Data_NA_09_10.mat', 'K', 'pos', 'Seq', 'tree')
    end


%% end of data
end

toc

