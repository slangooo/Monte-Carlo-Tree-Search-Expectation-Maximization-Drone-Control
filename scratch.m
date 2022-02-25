for i =1:length(statesHistory_K_05_N_5)-1
    all(statesHistory_K_05_N_5(i+1,:) == childHistory_K_05_N_5(i).sigma)
end


[length(childHistory_K_05_N_5),length(statesHistory_K_05_N_5(:,1))]


figure
plot(numberOfServedUsersPerIteration(2:end-1),'C')
hold on