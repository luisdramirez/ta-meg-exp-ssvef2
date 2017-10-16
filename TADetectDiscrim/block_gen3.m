% block_gen3.m

nContrastT1 = 4;
nContrastT2 = 4;
nV = 4;
nT = 2;

a = fullfact([nContrastT1 nContrastT2 nV nT]);

validity = a(:,3);
validity(validity~=4) = 1;
validity(validity==4) = 2;
a(:,3) = validity;

idx = randperm(size(a,1));
b = a(idx,:);





contrastT1 = a(:,1);
contrastT2 = a(:,2);
validity = a(:,3);
target = a(:,4);
