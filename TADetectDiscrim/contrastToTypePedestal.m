function [tt, tp] = contrastToTypePedestal(contrastLevels)

lut = [1 1; 1 2; 2 1; 2 2];

for i = 1:length(contrastLevels)
    if contrastLevels(i)==0
        tp(i) = NaN;
        tt(i) = NaN;
    else
        tp(i) = lut(contrastLevels(i),1);
        tt(i) = lut(contrastLevels(i),2);
    end
end

