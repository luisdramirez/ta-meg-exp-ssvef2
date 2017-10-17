% codetest.m

cueValidity = b.cueValidity;
responseTarget = b.responseTarget;
targetContrast = (b.targetPedestal-1)*2 + b.targetType;

count = [];
vs = [1 -1]; 
ts = [1 2];
cs = 1:4;
for iT = 1:2
    target = ts(iT);
    wT = responseTarget==target;
    for iV = 1:2
        validity = vs(iV);
        wV = cueValidity==validity;
        for iC = 1:4
            contrast = cs(iC);
            wC = targetContrast==contrast;
            count(iT,iV,iC) = nnz(wT & wV & wC);
        end
    end
end

