function [Kk, Ww, track] = epi(distribution_s,s0,a,L,N,tf,run,muL,f0L,epi,E_st,Ee, coef)
format long;

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

DIM = N*L;
f0=f0L/L;
muesc=0;                      % standing variation f0(s): deleterious mutation rate during escape
fmin=0.5/N;                   % cutoff  for average sample of 1/fmin seq. 
T=0:tf';                      % times
mu=muL/L;
% SUlogNs = s0/(mu*L)*log(N*s0);

%% parameter for epistasis
E1 = E_st;   
E2 = -E1;

% Distributions of selection coefficients
switch distribution_s
    case 'binary' 
        s=s0*(-ones(1,L) + 2*(rand(1,L) > 1-a) + 1e-5*(rand(1,L)-0.5)); 

       case 'half'
 %%
% close all

% L = 1000;

% clc
% s=-s0*log(rand(1,L)); 

% coef=.25;
pd2 = makedist('HalfNormal','mu',0,'sigma',coef*sqrt(pi/2));

s = random(pd2,L,1)'; 
s = -s;

MM = mean(s);
St = std(s);

% x = 1:L;
% sh = pdf(pd2,x); 
% scr = sh(randperm(length(sh)));

% s = (s - min(s)) ./ (max(s) - min(s));

% [aa,bb] = hist(s,7,'Normalization','pdf')
% % 2/pi
% 
% s = sort(s,'descend');
% 
% % % figure;
% subplot(1,2,1),histogram(s);
% 
% subplot(1,2,2),histogram(s,'Normalization','pdf')

% legend({'mu = 0, sigma = 1'},'Location','NE');
%  s = scr;     
        
end

Knew=zeros(N,L); Kk = zeros(N,L,tf); Ww = zeros(N,tf);

if f0~=0
    K=(rand(N,L) < f0);
%     sum(sum(K))
%     pause
elseif muesc~=0
	f0=muesc*s.^(-1).*(1-exp(-s*tesc));       
	K=(rand(N,L) < ones(N,1)*f0);
else
    K=zeros(N,L);
end

%% Epistatic interactions matrix
E = zeros(L,L) ;  

if epi == 1
    
for i=1:L
    for j=1:L
          E(i,j) = (abs(s(i))+abs(s(j)))*E1*Ee(i,j);
    end  
end

elseif epi == 2
    
for i=1:L
    for j=1:L
          E(i,j) = (abs(s(i))+abs(s(j)))*E2*Ee(i,j);
    end        
end

else
end
     tE = tril(E,-1);
    % Evolution starts...
   
    
    % genome tags
    Gen_tag = cell(N,tf+1);
    gen_tag = split(int2str([1:N]));
    
    Gen_tag(:,1) = gen_tag;
    
    tracking = ones(N,tf+1);
    
    
for t=T(2:end)
    
      Kk(:,:,t)=K;
    K = xor(K,sprand(N,L,mu));
  
    if run == 1 
    w = K*(s');
    elseif run == 2
    DIA = diag(K*tE*K');    %% GAB using lower triangular matrix of E to avoid counting each epistatic pair twice
    w = K*(s') + DIA;
    end
    
    % Random drift: divide in 2 or die 
    % Same as Poisson, Var[n]=E[n], Neff=N for drift/tree; not for few copies of an allele
    
    nprogav1=min(2,exp(w)/sum(exp(w)));      % average progeny numbers <= 2N
    nprogav2=nprogav1/mean(nprogav1);         % renormalizing
    nprog=rand(N,1) < nprogav2/2 ;          % survive or die?
    

    % balancing survivals and deaths to keep population constant
    disbal=1;
    while disbal ~=0
%         disp('non N')
        isur=find(nprog); 
        idied=find(1-nprog);
        
        disbal=length(isur) - N/2;
        if disbal > 0
            iflip=isur(ceil(rand(1,disbal)*length(isur)));  
            nprog(iflip)=zeros(1,disbal); % replacing extra 1 by 0
        elseif disbal < 0
            iflip=idied(ceil(rand(1,-disbal)*length(idied)));
            nprog(iflip)=ones(1,-disbal);
        end
    end
    nprog=nprog*2 ;

    tracking(:,t+1) = nprog;
    
    % updating population
    is=[0;cumsum(nprog(1:(N-1)))];
    
    oc = 1;
    
    for i=1:N
        if nprog(i)
            Knew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*K(i,:);
                   Gen_tag(oc:oc+1,t+1) = Gen_tag(i,t);   oc = oc+2;
        end
    end
    K=Knew; 
   
    sK=size(Knew);sK=sK(1);
    if sK  ~= N, disp('N not conserved'),return;end
    
      Ww(:,t) = w;
    
%     Gen_tag
%     tracking

    % Evolution step ended
end

end %Evolution ended





