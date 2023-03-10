clear all
close all
clc

Nimp=52; %number of features to be mixed

fun1=@FeatureRoutine1d; %construct 1 mixed feature
fun2=@FeatureRoutine2d; %construct 2 mixed features

lb = 0;
ub = 1;

maxminu=3*60; % Max durantion time [minutes]

optionspa = optimoptions('gamultiobj','PopulationSize',120,...
          'ParetoFraction',0.7,'PlotFcn',@gaplotpareto,'MaxTime',maxminu*60);

[solution,ObjectiveValue] = gamultiobj(fun2,2*Nimp,...
                         [],[],[],[],[],[],optionspa);

save Pareto2DT35Batt52Pow solution ObjectiveValue