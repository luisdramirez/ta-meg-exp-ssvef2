function turnallwhite
%
% turnallwhite
%
% turns the background color of all open figures white

fh = findobj('Type','figure');

for f = fh 
    set(f,'color','w')
end

