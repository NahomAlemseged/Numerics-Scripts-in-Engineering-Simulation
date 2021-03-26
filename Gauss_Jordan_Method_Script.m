%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE BY NAHOM A. WORKU
% GAUSS JORDAN METHOD OF SOLVING SYSTEMS OF LINEAR EQUATIONS
A = [1 3 2; 4 2 5; 3 3 1];  % COEFFICIENT MATRIX
B = [6;11;7];	% KNOWN VECTOR
n = length(A);
%% FORWARD ELIMINATION
for i=1:n-1
    for j=i+1:n
        m = A(j,i)/A(i,i);  %MULTIPLYING FACTOR m FOR MANIPULATION OF ROWS
      for k=1:n
          A(j,k) = A(j,k)-(m.*A(i,k));  %ROW ELIMINATION PROCESS FOR COEFFICINT MATRIX
      end
      B(j) = B(j) - m* B(i);	%ROW ELIMINATION PROCESS FOR COEFFICINT MATRIX
    end
end
%% BACKWARD SUBSTITUTION
C = [A B];
for i=n:-1:2
    for j = i-1:-1:1
        m = C(j,i)/C(i,i);
        C(j,:) = C(j,:) - m*C(i,:);
    end
end
for i=1:n
    X(i) = C(i,n+1)/C(i,i); 	% IDENTITY MATRIX
end
X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

