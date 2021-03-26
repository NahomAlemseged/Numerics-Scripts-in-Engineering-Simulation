% CODE BY NAHOM A. WORKU
% GAUSS JORDAN METHOD OF SOLVING SYSTEMS OF LINEAR EQUATIONS
A = [1 3 2; 4 2 5; 3 3 1];
B = [12;11;7];
n = length(A);
%% FORWARD ELIMINATION
for i=1:n-1
    for j=i+1:n
        m = A(j,i)/A(i,i);  %MULTIPLYING FACTOR m FOR MANIPULATION OF ROWS
      for k=1:n
          A(j,k) = A(j,k)-(m.*A(i,k));
      end
      B(j) = B(j) - m* B(i);
    end
end
%% BACKWARD ELIMINATION
C = [A B];
for i=n:-1:2
    for j = i-1:-1:1
        m = C(j,i)/C(i,i);
        C(j,:) = C(j,:) - m*C(i,:);
    end
end
for i=1:n
    X(i) = C(i,n+1)/C(i,i);
end
X