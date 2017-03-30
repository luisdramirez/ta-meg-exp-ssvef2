function rd_supertitle2(supertitle)
%
% rd_supertitle(supertitle)
%
% places supertitle as the title of a current figure. useful for subplots.

set(gcf,'NextPlot','add');
axes;
h = title(supertitle);
set(gca,'Visible','off');
set(h,'Visible','on');
rd_raiseAxis(gca);