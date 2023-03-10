clear all
clc
close all

load LISTnor01 %if you want to construct the mixed feature with 30 SHAP features
%load LISTnor02 if you want to construct the mixed feature with 52 SHAP features

Nsample=max(size(TRAINDATAnor));

soglia=35;          % [k]

smo=-1;   % SMOTE if smo>0

mixing=1; % 1->power; 2->linear

utopian=0;

Nbins=100;

load Pareto1DT35Batt30Pow.mat

salva=1;

%%%%%%%%%%%% UTOPIA POINT FROM PARETO FRONT %%%%%%%%%%%%%%%%%%%
m1=min(ObjectiveValue(:,1));
mm1=max(ObjectiveValue(:,1));

m2=min(ObjectiveValue(:,2));
mm2=max(ObjectiveValue(:,2));

[m3,pic2]=min(ObjectiveValue(:,3));

if m3>0
ObjectiveValue(:,3)=ObjectiveValue(:,3)-abs(m3);
else
ObjectiveValue(:,3)=ObjectiveValue(:,3)+abs(m3);
end

m3=min(ObjectiveValue(:,3));
mm3=max(ObjectiveValue(:,3));

ObjectiveValue(:,1)=(ObjectiveValue(:,1)-m1)/(mm1-m1);
ObjectiveValue(:,2)=(ObjectiveValue(:,2)-m2)/(mm2-m2);
ObjectiveValue(:,3)=(ObjectiveValue(:,3)-m3)/(mm3-m3);

mx=10000;
pic=-1;

for i=1:max(size(ObjectiveValue))

    if norm([ObjectiveValue(i,1),ObjectiveValue(i,2),...
        ObjectiveValue(i,3)])<mx
    
    mx=norm([ObjectiveValue(i,1),ObjectiveValue(i,2),...
        ObjectiveValue(i,3)]);

    pic=i;
    end

end
%%%%%%%%%%%% Least POINT FROM PARETO FRONT %%%%%%%%%%%%%%%%%%%

if utopian==0 
    pic=pic2;
end

x=solution(pic,:)/norm(solution(pic,:));

Nx=max(size(x));

indexSI=find(Tc>=soglia);
indexNO=find(Tc<soglia);

Nsi=max(size(indexSI));
Nno=max(size(indexNO));

FEATSI=zeros(Nsi,Nx);
FEATNO=zeros(Nno,Nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Nx
    FEATSI(:,i)=TRAINDATAnor(indexSI,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Nx
    FEATNO(:,i)=TRAINDATAnor(indexNO,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if smo>0
%%%%%%%%%%%%%%%%%%%%%%%MINORITY CLASS OVERSAMPLING%%%%%%%%%%%%%%%%%%%%%%%%%
SINTET=[FEATSI Tc(indexSI)];

[SINTET,C,Xn,Cn] = smote(SINTET,smo,smo);
FEATSI=SINTET(:,1:Nx);
Nsi=max(size(FEATSI));
%%%%%%%%%%%%%%%%%%%%%%%MINORITY CLASS OVERSAMPLING%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%CONSTRUCTION FIRST FEATURE%%%%%%%%

if mixing==1

    ff1=ones(Nsample,1);
    f1si=ones(Nsi,1);
    f1no=ones(Nno,1);

else

    ff1=zeros(Nsample,1);
    f1si=zeros(Nsi,1);
    f1no=zeros(Nno,1);

end

%%% CONSTRUCTION FIRST FEATURE

for i=1:Nx

    if mixing==1

    ff1=ff1.*TRAINDATAnor(:,i).^x(i);
    f1si=f1si.*FEATSI(:,i).^x(i);
    f1no=f1no.*FEATNO(:,i).^x(i);
    
     else
  
    ff1=ff1+x(i)*TRAINDATAnor(:,i);
    f1si=f1si+x(i)*FEATSI(:,i);
    f1no=f1no+x(i)*FEATNO(:,i);
    end

end
%%%%%%%%%%%%%%%CONSTRUCTION FIRST FEATURE%%%%%%%%

b1min=min(ff1);
b1max=max(ff1);

figure(1)

%%%%%
f1sinor=(f1si-b1min)/(b1max-b1min);

f1nonor=(f1no-b1min)/(b1max-b1min);

hsi=histogram(f1sinor,0:1/Nbins:1,'Normalization','probability');
%%%%%

H1=hsi.Values;

W1=hsi.Values;
f1=hsi.BinEdges;
F1=f1(1:end-1);

delta=hsi.BinWidth;

hold on
hno=histogram(f1nonor,0:1/Nbins:1,'Normalization','probability');

legend('IN','OUT')

H2=hno.Values;
W2=hno.Values;
f2=hno.BinEdges;
F2=f2(1:end-1);

[N1,N2]=size(H1);

Nsi=max(size(FEATSI));
Nno=max(size(FEATNO));


%--------------------------------------------------------------------------
% Costruisci dati per i classificatori
%--------------------------------------------------------------------------

if smo>0
    TRAINDATAnor=[FEATSI;FEATNO];
    Tc=[SINTET(:,end);Tc(indexNO)];
end

Nsample=max(size(TRAINDATAnor));

if mixing ==1
    TTT(:,1)=ones(Nsample,1);
else
    TTT(:,1)=zeros(Nsample,1);
end

for i=1:Nx
    
if mixing==1
    TTT(:,1)=TTT(:,1).*TRAINDATAnor(:,i).^x(i);
else
    TTT(:,1)=TTT(:,1)+x(i)*TRAINDATAnor(:,i);
end    
   
TTT(:,i+1)=TRAINDATAnor(:,i);

end
%
%
TTT(:,Nx+2)=Tc>=soglia;

TTT(:,1)=(TTT(:,1)-b1min)/(b1max-b1min);

save mixing01D x b1min b1max

%%%%%%%%%%%

musi=F1*W1';
muno=F2*W2';

stdsi=0;
stdno=0;

for i=1:Nbins

    stdsi=stdsi+W1(i)*(F1(i)-musi)^2;
    stdno=stdno+W2(i)*(F2(i)-muno)^2;

end 

disp(sqrt(stdsi))
disp(sqrt(stdno))

%%%%%%%%%
if salva==1
fid = fopen('Power_Batt_1d_35_least_newcopy.txt', 'wt');
 
     fprintf(fid,'%s \n', 'Min-Max values');
     fprintf(fid,'%f %f \n', b1min, b1max);
     
     fprintf(fid,'%s \n', ' ');
     
     if mixing==1
     fprintf(fid,'%s \n', 'Weights-Power mixing');
     else
     fprintf(fid,'%s \n', 'Weights-Linear mixing');
     end
     
     fprintf(fid,'%12.12f \n',x);      

fclose(fid);
end