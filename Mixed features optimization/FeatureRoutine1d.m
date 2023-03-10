function difff=FeatureRoutine1d(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soglia=35;         % [k]
mixing=2;          % Mixing rule:               1->power; 2->linear
distancesel=2;     % Adopted distance:          1->EMD; 2-> Batt; 3->Mutual info
Nbins=35;          % sqrt(Nno);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

difff=zeros(1,3);

x=x/norm(x);

load LISTnor01

Nsample=max(size(TRAINDATAnor));

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Nx
FEATNO(:,i)=TRAINDATAnor(indexNO,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mixing==1

    ff1=ones(Nsample,1);
    f1si=ones(Nsi,1);
    f1no=ones(Nno,1);

else
    
    ff1=zeros(Nsample,1);
    f1si=zeros(Nsi,1);
    f1no=zeros(Nno,1);

end

%%% costruzione prima feature mista

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

b1min=min(ff1);
b1max=max(ff1);

fig=figure; 
set(fig,'visible','off');

f1sinor=(f1si-b1min)/(b1max-b1min);
f1nonor=(f1no-b1min)/(b1max-b1min);

hsi=histogram(f1sinor,0:1/Nbins:1,'Normalization','probability');

delta=hsi.BinWidth;
W1=hsi.Values;
f1=hsi.BinEdges;
F1=f1(1:end-1)+delta/2;

close(fig);

fig=figure; 
set(fig,'visible','off');
 
hno=histogram(f1nonor,0:1/Nbins:1,'Normalization','probability');

W2=hno.Values;
f2=hno.BinEdges;
F2=f2(1:end-1)+delta/2;

close(fig);

if distancesel==1

    [~,fval] = emd(F1',F2',W1',W2',@gdf); % EMD
    difff(3)=-fval;

elseif distancesel==2

    N1=max(size(W1));

    measure=zeros(1,N1);

    for i=1:N1

        measure(i)=sqrt(W1(i)*W2(i)); % Bhattacharyya distance
    
    end

    difff(3)=sum(sum(measure));

    
else

    WW1=W1*Nbins;
    WW2=W2*Nbins;
    difff(3)=mutualinfo(WW1,WW2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

musi=F1*W1';
muno=F2*W2';

stdsi=0;
stdno=0;

for i=1:Nbins

stdsi=stdsi+W1(i)*(F1(i)-musi)^2;
stdno=stdno+W2(i)*(F2(i)-muno)^2;

end 

difff(1)=sqrt(stdsi);
difff(2)=sqrt(stdno);