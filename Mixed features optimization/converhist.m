function [F,W]=converhist(H1)

[N1,N2]=size(H1);

x=0:1/N1:1;
y=x;

cont=0;

for i=1:N1

    for j=1:N2

        cont=cont+1;

        F(cont,1)=x(i);
        
        F(cont,2)=y(j);
    
        W(cont)=H1(i,j);
    end

end