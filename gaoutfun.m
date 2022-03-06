function [state,options,optchanged] = gaoutfun(options,state,~)
writematrix(state.Population,'CurrentPopulation.txt','Delimiter',' ');
writematrix(state.Score,'CurrentPopulationScore.txt','Delimiter',' ');
writematrix(state.Population,'TotalPopulation.txt','Delimiter',' ','WriteMode','append');
writematrix(state.Score,'TotalPopulationScore.txt','Delimiter',' ','WriteMode','append');
optchanged = false;
end