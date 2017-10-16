function [blockOrder, attBlockOrder, targetBlockOrder, cueBlockOrder, targetTypeBlockOrder, targetPedestalBlockOrder] = ...
    block_gen2(blockNames, attBlockNames, targetBlockNames, cueBlockNames, run)

%% setup
blank = find(ismember(blockNames,'blank'));
fast_side = find(ismember(blockNames,{'fast-left','slow-left'}));

no_att = find(ismember(attBlockNames,'no-att'));
att_right = find(ismember(attBlockNames,'att-right'));

nt = find(ismember(targetBlockNames,'no-targ'));
nc = find(ismember(cueBlockNames,'no-cue'));

%% targetBlockOrder (no-targ, pres-pres)
targetBlockOrder = 2*ones(1,32);

%% cueBlockOrder
m = mod(run,4)-1;
b = [0 0 0 1];
bb = [];
for i = 1:length(b)
    bb = [bb circshift(b,[i+m,0])];
end

cueBlockOrder = [bb + 2, bb + 4];

%% targetTypeBlockOrder, targetPedestalBlockOrder
% precue, T1, T2
a = fliplr(fullfact([4 4 2]));

[t1TT, t1TP] = contrastToTypePedestal(a(:,2));
[t2TT, t2TP] = contrastToTypePedestal(a(:,3));

targetTypeBlockOrder = [t1TT' t2TT'];
targetPedestalBlockOrder = [t1TP' t2TP'];

%% randomize trial order
% indices = 1:length(cueBlockOrder); % for debugging
indices = randperm(length(cueBlockOrder));
cueBlockOrder = cueBlockOrder(indices);
targetBlockOrder = targetBlockOrder(indices);
targetTypeBlockOrder = targetTypeBlockOrder(indices,:)';
targetPedestalBlockOrder = targetPedestalBlockOrder(indices,:)';

%% block order (one condition) and attention order
blockOrder = repmat(fast_side,1,length(cueBlockOrder));
attBlockOrder = repmat(att_right,1,length(cueBlockOrder));

%% insert blank trials 
% (every 4 target trials for target and cue block order)
blank_tar = repmat(nt,[1,(length(targetBlockOrder)/4) + 1]);
ind = zeros(1, length(targetBlockOrder)+ length(blank_tar));
ind (1:5:length(ind)) = blank_tar;
ind2 = ind;
ind3 = repmat(ind,2,1);
ind3b = repmat(ind,2,1);
ind4 = ind;
ind5 = ind;
ind3(ind3==1) = NaN;
ind3b(ind3b==1) = NaN;
ind2 (ind == nt ) = nc;

ind(ind == 0) = targetBlockOrder;
targetBlockOrder = ind;

ind2(ind2 == 0) = cueBlockOrder;
cueBlockOrder = ind2;

ind3(1, ind3(1,:) == 0) = targetTypeBlockOrder(1,:);
ind3(2, ind3(2,:) == 0) = targetTypeBlockOrder(2,:);
targetTypeBlockOrder = ind3;

ind3b(1, ind3b(1,:) == 0) = targetPedestalBlockOrder(1,:);
ind3b(2, ind3b(2,:) == 0) = targetPedestalBlockOrder(2,:);
targetPedestalBlockOrder = ind3b;

ind4(ind4==1) = blank;
ind4(ind4 == 0) = blockOrder;
blockOrder = ind4;

ind5(ind5==1) = no_att;
ind5(ind5 == 0) = attBlockOrder;
attBlockOrder = ind5;
