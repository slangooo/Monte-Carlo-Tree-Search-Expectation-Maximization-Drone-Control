uload('simulation_results_all.mat')
figure
hold on
plot(results{1}{7},'B')
plot(results{1}{6},'R')
plot(results14{1}{7},'B.')
plot(results14{1}{6},'R.')



figure
hold on
plot(results14{1}{9},'B')
plot(results15{1}{9},'R')
plot(results16{1}{9},'C')
plot(results17{1}{9},'G')



plot(resultsOld2{4}{6},'G')

load('simulation_results_all.mat')
figure
hold on
idx = 8;
plot(results{3}{1}(idx,:),'B')
plot(results{3}{2}(idx,:),'R')
plot(resultsOld{3}{6},'G')

sinrThresholds = -5:22;

(results{1}{4}/500000).*(results{1}{8}.*15)/500