clear all
clc
close all

%load LISTnor01 if you want to construct the mixed feature with 30 SHAP features
load LISTnor02 %if you want to construct the mixed feature with 52 SHAP features


Nbins=35;

smo=-1;  % SMOTE parameter: If smo>o performs SMOTE

mixing=1; % 1->power ; 2->linear

Nsample=max(size(TRAINDATAnor));

soglia=35;          % [k]

utopian=0;  % Do you want to select the Utopian point (1) or the least distance (0)?

load Pareto2DT35Batt52Pow.mat

salva=0;

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
%%%%%%%%%%%% UTOPIA POINT FROM PARETO FRONT %%%%%%%%%%%%%%%%%%%

if utopian==0 
    pic=pic2;
end

x=solution(pic,:);

figure(1)
plot3(ObjectiveValue(:,1),ObjectiveValue(:,2),...
    ObjectiveValue(:,3),'m*')
grid on
box on

hold on
plot3(ObjectiveValue(pic,1),ObjectiveValue(pic,2),...
    ObjectiveValue(pic,3),'sg','MarkerSize',13,'LineWidth',3)
legend('Pareto Front','Selection')

figure(2)

Nx=max(size(x));

x1=x(1:Nx/2);
x2=x(Nx/2+1:Nx);

x1=x1/norm(x1);
x2=x2/norm(x2);

x(1:Nx/2)=x1(:);
x(Nx/2+1:Nx)=x2(:);

indexSI=find(Tc>=soglia);
indexNO=find(Tc<soglia);

Nsi=max(size(indexSI));
Nno=max(size(indexNO));

FEATSI=zeros(Nsi,Nx/2);
FEATNO=zeros(Nno,Nx/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Nx/2
    FEATSI(:,i)=TRAINDATAnor(indexSI,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Nx/2
    FEATNO(:,i)=TRAINDATAnor(indexNO,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if smo>0
SINTET=[FEATSI Tc(indexSI)];

%%%%%%%%%%%%%%%%%%%%%%%MINORITY CLASS OVERSAMPLING%%%%%%%%%%%%%%%%%%%%%%%%%
[SINTET,C,Xn,Cn] = smote(SINTET,smo,ceil(smo));
FEATSI=SINTET(:,1:Nx/2);
Nsi=max(size(FEATSI));
%%%%%%%%%%%%%%%%%%%%%%%MINORITY CLASS OVERSAMPLING%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%% costruzione prima feature mista

if mixing==1
    ff1=ones(Nsample,1);
    f1si=ones(Nsi,1);
    f1no=ones(Nno,1);
else
    ff1=zeros(Nsample,1);
    f1si=zeros(Nsi,1);
    f1no=zeros(Nno,1);
end


for i=1:Nx/2

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

%%% costruzione SECONDA feature mista

if mixing==1
    ff2=ones(Nsample,1);
    f2si=ones(Nsi,1);
    f2no=ones(Nno,1);
else
    ff2=zeros(Nsample,1);
    f2si=zeros(Nsi,1);
    f2no=zeros(Nno,1);
end

for i=1:Nx/2
    if mixing==1
        ff2=ff2.*TRAINDATAnor(:,i).^x(i+Nx/2);
        f2si=f2si.*FEATSI(:,i).^x(i+Nx/2);
        f2no=f2no.*FEATNO(:,i).^x(i+Nx/2);
    else
        ff2=ff2+x(i+Nx/2)*TRAINDATAnor(:,i);
        f2si=f2si+x(i+Nx/2)*FEATSI(:,i);
        f2no=f2no+x(i+Nx/2)*FEATNO(:,i);
    end
end

b1min=min(ff1);
b1max=max(ff1);

b2min=min(ff2);
b2max=max(ff2);

figure(2)

f1sinor=(f1si-b1min)/(b1max-b1min);
f2sinor=(f2si-b2min)/(b2max-b2min);

f1nonor=(f1no-b1min)/(b1max-b1min);
f2nonor=(f2no-b2min)/(b2max-b2min);

hsi=histogram2(f1sinor,f2sinor,0:1/Nbins:1,0:1/Nbins:1,...
    'Normalization','probability');

H1=hsi.Values;

yy1=0:1/(Nbins-1):1;
yy2=0:1/(Nbins-1):1;

[XX1,YY1] = meshgrid(yy1(1:end),yy2(1:end));

hold on

hno=histogram2(f1nonor,f2nonor,0:1/Nbins:1,0:1/Nbins:1,...
    'Normalization','probability','FaceColor','yellow');

H2=hno.Values;

zz1=0:1/(Nbins-1):1;
zz2=0:1/(Nbins-1):1;

[XX2,YY2] = meshgrid(zz1(1:end),zz2(1:end));

xlabel('First feature')
ylabel('Second feature')
legend('IN','OUT')

[F2,W2]=converhist(H2);
[F1,W1]=converhist(H1);

Nsi=max(size(FEATSI));
Nno=max(size(FEATNO));

if smo>0
% Costruisci dati per i classificatori
TRAINDATAnor=[FEATSI;FEATNO];
Tc=[SINTET(:,end);Tc(indexNO)];
Nsample=max(size(TRAINDATAnor));
end

if mixing==1
TTT(:,1)=ones(Nsample,1);
TTT(:,2)=ones(Nsample,1);
else
TTT(:,1)=zeros(Nsample,1);
TTT(:,2)=zeros(Nsample,1);
end

for i=1:Nx/2
    
if mixing==1
    TTT(:,1)=TTT(:,1).*TRAINDATAnor(:,i).^x(i);
    TTT(:,2)=TTT(:,2).*TRAINDATAnor(:,i).^x(i+Nx/2);
else
    TTT(:,1)=TTT(:,1)+x(i)*TRAINDATAnor(:,i);
    TTT(:,2)=TTT(:,2)+x(i+Nx/2)*TRAINDATAnor(:,i);
end
    TTT(:,i+2)=TRAINDATAnor(:,i);
end
%
TTT(:,Nx/2+3)=Tc>=soglia;


TTT(:,1)=(TTT(:,1)-b1min)/(b1max-b1min);
TTT(:,2)=(TTT(:,2)-b2min)/(b2max-b2min);

%save TESTCLASS/mixing02D x b1min b1max b2min b2max

%%%%%%%%%%%REMOVE ZEROS%%%%%%%%%%%%%
% for i=1:max(size(Tc))
% 
%     if Tc(i)==0
%     
%         TTT(i,1)=NaN;
%         TTT(i,2)=NaN;
%         Tc(i)=NaN;
%     end
% 
% end
%%%%%%%%%%REMOVE ZEROS%%%%%%%%%%%%%

figure
scatter(TTT(:,1),TTT(:,2),[],Tc,'filled')
colorbar

entr1=BivariateGauss(XX1,YY1,H1);
entr2=BivariateGauss(XX2,YY2,H2);

%%%%%%%%%
if salva==1
fid = fopen('Pow_Batt_2D_uto_mixed_52.txt', 'wt');
 
     fprintf(fid,'%s \n', 'Min-Max values (first)');
     fprintf(fid,'%f %f \n',b1min, b1max);
     

     fprintf(fid,'%s \n', ' ');
     
     if mixing==1
     fprintf(fid,'%s \n', 'Weights-Power mixing (first)');
     else
     fprintf(fid,'%s \n', 'Weights-Linear mixing (first)');
     end
     
     fprintf(fid,'%12.12f \n',x1);
    
     fprintf(fid,'%s \n', ' ');
    
     fprintf(fid,'%s \n', 'Min-Max values (second)');
     fprintf(fid,'%f %f \n',b2min, b2max);
     fprintf(fid,'%s \n', ' ');
     if mixing==1
     fprintf(fid,'%s \n', 'Weights-Power mixing (second)');
     else
     fprintf(fid,'%s \n', 'Weights-Linear mixing (second)');
     end
     
     fprintf(fid,'%12.12f \n',x2);

     fclose(fid);
end