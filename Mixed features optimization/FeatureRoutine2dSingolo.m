function difff=FeatureRoutine2dSingolo(x)

%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%
soglia=25;      % [K]  ATTENTION: Give here the right temperature!!!

mixing=2;       % Mixing rule: 1->power; 2->linear

r0=1e-2;        % radius for the kNN search. 
                % We suggest to avoid too coars parameter; in this work, we have used 1e-2.
%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%

Nx=max(size(x));
% 
x1=x(1:Nx/2);
x2=x(Nx/2+1:Nx);

%%%%%Gram orthonormalization%%%%
x2=x2-((x1*x2')/(x1*x1'))*x1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1=x1/norm(x1);
x2=x2/norm(x2);

x(1:Nx/2)=x1(:);
x(Nx/2+1:Nx)=x2(:);

load LISTnor02  % load the training set

                
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

f1sinor=(f1si-b1min)/(b1max-b1min);
f2sinor=(f2si-b2min)/(b2max-b2min);

f1nonor=(f1no-b1min)/(b1max-b1min);
f2nonor=(f2no-b2min)/(b2max-b2min);


Sii=[f1sinor f2sinor];
Noo=[f1nonor f2nonor];


fid = fopen('Pow_Batt_2D_knn_mixed_52.txt', 'wt');
 
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
    
difff=k2NN(Sii,Noo,r0);  