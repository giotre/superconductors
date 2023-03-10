function difff=FeatureRoutine2d(x)

%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%
soglia=35;          % [k]

mixing=1;           % Mixing rule:        1->power; 2->linear

distancesel=2;      % Adopted distance:   1->EMD;   2-> Batt 

if distancesel==1
    Nbins=8;           % round(sqrt(Nsi));
elseif distancesel==2
    Nbins=15;
end
%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%

Nx=max(size(x));
% 
x1=x(1:Nx/2);
x2=x(Nx/2+1:Nx);

x1=x1/norm(x1);
x2=x2/norm(x2);

x(1:Nx/2)=x1(:);
x(Nx/2+1:Nx)=x2(:);

load LISTnor02

Nsample=max(size(TRAINDATAnor));

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

%%% costruction first mixed feature

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

%%% costruction second mixed feature

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

fig=figure; 
set(fig,'visible','off');

f1sinor=(f1si-b1min)/(b1max-b1min);
f2sinor=(f2si-b2min)/(b2max-b2min);

f1nonor=(f1no-b1min)/(b1max-b1min);
f2nonor=(f2no-b2min)/(b2max-b2min);

hsi=histogram2(f1sinor,f2sinor,0:1/Nbins:1,0:1/Nbins:1,...
    'Normalization','probability');

H1=hsi.Values;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yy1=0:1/(Nbins-1):1;
yy2=0:1/(Nbins-1):1;

[XX1,YY1] = meshgrid(yy1(1:end),yy2(1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close(fig);
%hold on
fig=figure; 
set(fig,'visible','off');
 
hno=histogram2(f1nonor,f2nonor,0:1/Nbins:1,0:1/Nbins:1,...
    'Normalization','probability');

H2=hno.Values;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zz1=0:1/(Nbins-1):1;
zz2=0:1/(Nbins-1):1;

[XX2,YY2] = meshgrid(zz1(1:end),zz2(1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close(fig);

[F2,W2]=converhist(H2);
[F1,W1]=converhist(H1);

if distancesel==1
   
    [~,fval]=emd(F1,F2,W1',W2',@gdf); % Earth mover distance
    difff(3)=-fval;

elseif distancesel==2

    [N1,N2]=size(H1);

    measure=zeros(N1,N2);

    for i=1:N1
   
        for j=1:N2
           measure(i,j)=sqrt(H1(i,j)*H2(i,j)); % Bhattacharyya distance
        end
   
    end

    difff(3)=sum(sum(measure));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

difff(1)=BivariateGauss(XX1,YY1,H1);
difff(2)=BivariateGauss(XX2,YY2,H2);
