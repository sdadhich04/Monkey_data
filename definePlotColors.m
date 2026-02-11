
%A. Orsborn, updated 3-14-11
%defines a color scheme for plots (defined with matrix plotColors (colors x 3)
%sets these colors as the default color order for multi-line plots 
%also defines linestyle order to 'solid', 'dashed', 'dotted' so that lines can
%be distinguished if > plotColors lines are plotted on the same figure

plotColors = zeros(9,3);

plotColors(1,:) = [0 0 0];      %dark gray
plotColors(2,:) = [0 0.4 1];          %blue
plotColors(3,:) = [0.2 0.7 0.7];      %cyan
plotColors(4,:) = [0.7 0.7 0.7];       %gray
plotColors(5,:) = [0.5 .2 0.9];       %purple
plotColors(6,:) = [1 0.7 0.4];        %orange  pale orange [1 0.95 0.7]
plotColors(7,:) = [0.1 0.7 0.2];      %green
plotColors(8,:) = [0.95 0.25 0.25];  %red
plotColors(9,:) = [1 0.85 0.3];      %yellow/gold
plotColors(10,:) = [0.7 0.1 0.7];      %magenta
plotColors(11,:) = [0.95 0.55 0.9];    %pink


set(0, 'DefaultAxesColorOrder', plotColors)
set(0, 'DefaultAxesLineStyleOrder', '-|--')