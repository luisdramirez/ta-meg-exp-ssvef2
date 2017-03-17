function [blockOrder, attBlockOrder, targetBlockOrder, cueBlockOrder, targetTypeBlockOrder] = block_gen(blockNames, attBlockNames, targetBlockNames, cueBlockNames, run)
% Random block generator for makeTADetectStim: one run (32 trials in one 
% repetition = 4 trials for each target (4) x cue (2) condition ) with
% added blank trials every 4 target trials 
% Assume cue validity = 75% ( 3 trials valid and 1 trial invalid )
% One blockOrder condition: fast-left or slow-left
% One attBlockOrder condition: att-right

if nargin < 5
    run = 1;
end

%% define indices
blank = find(ismember(blockNames,'blank'));
fast_side = find(ismember(blockNames,{'fast-left','slow-left'}));

no_att = find(ismember(attBlockNames,'no-att'));
att_right = find(ismember(attBlockNames,'att-right'));

nt = find(ismember(targetBlockNames,'no-targ'));
pp = find(ismember(targetBlockNames,'pres-pres'));
pa = find(ismember(targetBlockNames,'pres-abs'));
ap = find(ismember(targetBlockNames,'abs-pres'));
aa = find(ismember(targetBlockNames,'abs-abs'));

nc   = find(ismember(cueBlockNames,'no-cue'));
c1c1 = find(ismember(cueBlockNames,'1-1'));
c1c2 = find(ismember(cueBlockNames,'1-2'));
c2c1 = find(ismember(cueBlockNames,'2-1'));
c2c2 = find(ismember(cueBlockNames,'2-2'));

%% cue order
% pre-cue = T1
A = [c1c1,c1c1,c1c1,c1c2]; % specify validity for pre-cue = T1
cueBlockOrder_cue1 = repmat(A,1,4); % for all target conditions
% 
% pre-cue = T2
B = [c2c2,c2c2,c2c2,c2c1]; % specify validity for cue = T2
cueBlockOrder_cue2 = repmat(B,1,4); % for all target conditions

cueBlockOrder = [cueBlockOrder_cue1 , cueBlockOrder_cue2];

% randomize cue order
indices = 1:length(cueBlockOrder);
% indices = randperm(length(cueBlockOrder));
cueBlockOrder = cueBlockOrder(indices);

%% target order
dummy = repmat([pp,pa,ap,aa],[4,1]);
targetBlockOrder_cue1 = dummy(:)';
% targetBlockOrder = repmat(targetBlockOrder_cue1,[1,2]);
targetBlockOrder = repmat(targetBlockOrder_cue1,[1,numel(cueBlockOrder)/numel(targetBlockOrder_cue1)]);

% randomize target order
targetBlockOrder = targetBlockOrder(indices);

%% target type order
% hard-coding this for lack of a better method
% zero means still unassigned, nan means target absent

% assigning types for post-cued targets
% attempting to evenly distribute target types across conditions, but need
% 2 runs to do this because there is an odd number of targets per condition
% in each run
if mod(run,2) % odd runs
    targetTypeOrder(:,1) = [1 1 2 0 2 2 1 0 NaN NaN NaN NaN NaN NaN NaN NaN ...
        0 0 0 2 0 0 0 1 NaN NaN NaN NaN NaN NaN NaN NaN]'; % T1
    targetTypeOrder(:,2) = [0 0 0 1 NaN NaN NaN NaN 0 0 0 2 NaN NaN NaN NaN ...
        2 2 1 0 NaN NaN NaN NaN 1 1 2 0 NaN NaN NaN NaN]'; % T2
else % even runs
    targetTypeOrder(:,1) = [2 2 1 0 1 1 2 0 NaN NaN NaN NaN NaN NaN NaN NaN ...
        0 0 0 1 0 0 0 2 NaN NaN NaN NaN NaN NaN NaN NaN]'; % T1
    targetTypeOrder(:,2) = [0 0 0 2 NaN NaN NaN NaN 0 0 0 1 NaN NaN NaN NaN ...
        1 1 2 0 NaN NaN NaN NaN 2 2 1 0 NaN NaN NaN NaN]'; % T2   
end

% assigning types for non-targets
% randomly distributing equal number of each target type among all
% nontarget positions. not feasible to have a full factorial design and
% want target types for T1 and T2 to be independent
nontargetTypes = Shuffle([ones(1,8) 2*ones(1,8)]');
targetTypeOrder(targetTypeOrder==0) = nontargetTypes;

% set NaN to 0 to match previous style
targetTypeOrder(isnan(targetTypeOrder)) = 0;

targetTypeBlockOrder = targetTypeOrder(indices,:)';

%% insert blank trials 
% (every 4 target trials for target and cue block order)
blank_tar = repmat(nt,[1,(length(targetBlockOrder)/4) + 1]);
ind = zeros(1, length(targetBlockOrder)+ length(blank_tar));
ind (1:5:length(ind)) = blank_tar;
ind2 = ind;
ind3 = repmat(ind,2,1);
ind3(ind3==1) = NaN;
ind2 (ind == nt ) = nc;
ind(ind == 0) = targetBlockOrder;
targetBlockOrder = ind;
ind2(ind2 == 0) = cueBlockOrder;
cueBlockOrder = ind2;
ind3(1, ind3(1,:) == 0) = targetTypeBlockOrder(1,:);
ind3(2, ind3(2,:) == 0) = targetTypeBlockOrder(2,:);
targetTypeBlockOrder = ind3;

%% block order (one condition) and attention order
block = [blank,repmat(fast_side,1,4)];
blockOrder = [repmat(block,1,8),blank];

att = [no_att,repmat(att_right,1,4)];
attBlockOrder = [repmat(att,1,8),no_att];

end

