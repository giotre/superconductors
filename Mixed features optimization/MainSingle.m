clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Here we start from first tentative solutions (obtained with Genetic
%%%%%% algorithm) and minimize, according to kNN, the number of minimizing 
%%%%%% the number of negative neighbors around a positive sample in a given 
%%%%%% radius.
%%%%%% The minimization here is single objective, so we do not consider
%%%%%% anymore the variance of the distributions.
%%%%%% ATTENTION: any time one wants to change a threshold temperature, it
%%%%%% has to be changed also in the routine "FeatureRoutine2dSingolo".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nimp=52;

fun2Singolo=@FeatureRoutine2dSingolo;

lb = 0;

ub = 1;

maxminu=8*60; % Max durantion time [minutes]

optionspatt = optimoptions('patternsearch','Display','iter','PlotFcn','psplotbestf');

load OPTIMAL2D/Singolo2DT25kNN52FeatLinfineGS % Here we load a first tentative solution 
                                              % to be refined by
                                              % patternsearch
                                              % Some of the first tentative
                                              % solutions are in the folder
                                              % OPTIMAL2D
                                              

x=X;

Nx=max(size(x));

x0=zeros(1,Nx);

x1=x(1:Nx/2);
x2=x(Nx/2+1:Nx);

%%%%%Gram orthonormalization%%%%
x2=x2-((x1*x2')/(x1*x1'))*x1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1=x1/norm(x1);
x2=x2/norm(x2);

x0(1:Nx/2)=x1(:);
x0(Nx/2+1:Nx)=x2(:);

X = patternsearch(fun2Singolo,x0,[],[],[],[],[],[],optionspatt); % USO PATTERNSEARCH

save OPTIMAL2D/Singolo2DT25kNN52FeatLinfineGSpattNew.mat X