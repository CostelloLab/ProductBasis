% plotBoolEvolution() plots the evolution of the genes in the variable bHistory

figure
imagesc(1-bHistory')
colormap(gray)
xlabel('time \rightarrow')

ax = gca;

ax.XTick = 1:size(bHistory, 1);
ax.YTick = 1:size(bHistory, 2);

hold on
for loopLine = 1:(size(bHistory, 1)-1)
    plot(loopLine+[0 0]+0.5, [0 size(bHistory, 2)]+0.5, ':')
end
for loopLine = 1:(size(bHistory, 2)-1)
    plot([0 size(bHistory, 1)]+0.5, loopLine+[0 0]+0.5, ':')
end

ax.YTickLabel = allVars;
%ax.YTickLabel = { 'SLP', 'wg', 'WG', 'en', 'EN', 'hh', 'HH', 'ptc', 'PTC', 'PH', 'SMO', 'ci', 'CI', 'CIA', 'CIR' };