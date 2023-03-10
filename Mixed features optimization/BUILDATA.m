clear all
close all
clc

TRAINDATA = table2array(readtable('Train.xlsx'));
TESTDATA = table2array(readtable('Test.xlsx'));

testsoglia=35;

salva=1;

[Nsample,NF]=size(TRAINDATA);

%Nimp=30; % Importance = 70%
Nimp=52; % Importance = 90%

BOUND=zeros(Nimp,2);

for i=1:Nimp
    BOUND(i,1)=min(TRAINDATA(:,i));
    BOUND(i,2)=max(TRAINDATA(:,i));
end

TRAINDATAnor=zeros(Nsample,Nimp);
TESTDATAnor=zeros(max(size(TESTDATA)),Nimp);
Tc=zeros(1,Nsample);

%%%%%%%%%%FEATURE NORMALIZATION into [1-2]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for i=1:Nimp
    TRAINDATAnor(:,i)=(TRAINDATA(:,i)-BOUND(i,1))/...
        (BOUND(i,2)-BOUND(i,1))+1;
    TESTDATAnor(:,i)=(TESTDATA(:,i)-BOUND(i,1))/...
        (BOUND(i,2)-BOUND(i,1))+1;
end
%%%%%%%%%%FEATURE NORMALIZATION into [1-2]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

Tc=TRAINDATA(:,82);
Tct=TESTDATA(:,82);

%
if salva==1

    save LISTnor02.mat TRAINDATAnor Tc

    %save LISTTESTnor TESTDATAnor Tct

end
%
disp(100*sum(Tc>=testsoglia)/max(size(Tc)))
%
disp(100*sum(Tct>=testsoglia)/max(size(Tct)))


fid = fopen('PhysBound.txt', 'wt');
 
     fprintf(fid,'%s \n', 'Min             Max physical values of features');
     
     fprintf(fid,'%s \n', ' ');

     for i=1:Nimp
         fprintf(fid,'%12.12f %12.12f \n', BOUND(i,1), BOUND(i,2));
     end
   
   
fclose(fid);