% Bivariate Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all
%clc
%close all
%
% mu = [0 0];
% Sigma = [1 1; 1 1];
% 
% x1 = -5:0.2:5;
% x2 = -5:0.2:5;
% [X1,X2] = meshgrid(x1,x2);
% X = [X1(:) X2(:)];
% 
% y = mvnpdf(X,mu,Sigma);
% y = reshape(y,length(x2),length(x1));
% y=y/sum(sum(y));
% 
% surf(x1,x2,y)
% caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
% axis([-3 3 -3 3 0 0.4])
% xlabel('x1')
% ylabel('x2')
% zlabel('Probability Density')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function std=BivariateGauss(X1,X2,y)

mux=0;

for j=1:max(size(X1))

    mux=mux+X1(:,j)'*y(:,j);

end

muy=0;

for j=1:max(size(X2))

    muy=muy+X2(:,j)'*y(:,j);

end

% disp('mux:')
% disp(mux)
% 
% disp('muy:')
% disp(muy)

covxx=0;
covyy=0;
covxy=0;

for i=1:max(size(X1))

    for j=1:max(size(X2))

        covxx=covxx+y(i,j)*(X1(i,j)-mux)^2;
        covyy=covyy+y(i,j)*(X2(i,j)-muy)^2;
        covxy=covxy+y(i,j)*(X1(i,j)-mux)*(X2(i,j)-muy);
    
    end

end

Cov=[covxx,covxy;covxy,covyy];

%disp(Cov)

std=det(Cov);