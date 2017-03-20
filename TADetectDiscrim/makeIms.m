function [madeBgIms, madeTargIms] = makeIms(display, backgroundIms, maskedIms)
%make contrast target images

madeBgIms = cell(size(backgroundIms)); %backgroundIms = phase
madeTargIms = cell(size(maskedIms)); %maskedIms = position x phase x condition

%make backgroundIms Textures
for iPhase=1:length(backgroundIms) % phase
    bgtex = Screen('MakeTexture',display.windowPtr,backgroundIms{iPhase}*255);
    madeBgIms{iPhase} = bgtex;
end

%make maskedIms Textures
for iPos=1:size(maskedIms,1) %position
    for iPhase=1:size(maskedIms,2) %phase
        for iCond=1:size(maskedIms,3) %condition
            tex = Screen('MakeTexture',display.windowPtr, maskedIms{iPos, iPhase, iCond}*255);
            madeTargIms{iPos, iPhase, iCond} = tex;
        end
    end
end