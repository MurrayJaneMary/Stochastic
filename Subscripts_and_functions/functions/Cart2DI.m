%Calculates Decs, Incs, and Rs from a series of XYZs
%Input A = [X Y Z]

function  B = Cart2DI(A)

X=(A(:,1));
Y=(A(:,2));
Z=(A(:,3));

R = sqrt(X.^2 + Y.^2 + Z.^2);

Xu = X./R; Yu = Y./R; Zu = Z./R;

N=size(A,1);

for i = 1:N;

    if X(i) >= 0
        if Y(i) >= 0
            Dec(i,:) = degrees(atan(Yu(i)./Xu(i) )); 
            %disp('ONE!!!')
        else
            Dec(i,:) = degrees(atan(Yu(i)./Xu(i) ))+360; 
            %disp('TWO!!!')
        end

    else

        Dec(i,:) = degrees(atan(Yu(i) ./Xu(i) ))+180; 
        %disp('THREE!!!')
   
    end

end

Inc = degrees(asin(Zu));
%Dec

B = [Dec Inc R];
    
