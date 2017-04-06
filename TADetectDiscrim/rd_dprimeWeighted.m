function [dprime, criterion] = rd_dprimeWeighted(nh,nfa,nsignal,nnoise,faWeight)

% [dprime, criterion] = rd_dprimeWeighted(nh,nfa,nsignal,nnoise)
%
% nh = number of hit trials
% nfa = number of false alarm trials
% nsignal = number of signal trials
% nnoise = number of noise trials

% inputs
if nargin<5
    faWeight = 1;
end

% adjust for ceiling or floor values
nh(nh==nsignal) = nsignal(nh==nsignal)-1;
nh(nh==0) = 1;
nfa(nfa==nnoise) = nnoise(nfa==nnoise)-1;
nfa(nfa==0) = 1;

% proportions
h = nh./nsignal;
fa = nfa./nnoise;

% dprime
zh = norminv(h,0,1); zfa = norminv(fa,0,1)*faWeight;

dprime = zh - zfa;
criterion = -0.5*(zh+zfa);
