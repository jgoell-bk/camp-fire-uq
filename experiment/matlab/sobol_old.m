function [Si, STi] = sobol(fun, A, B)

% Computes first order Sobol indices using hybrid point Monte-Carlo:
%   - Si:  main effects Sobol indices
%   - STi: total effects Sobol indices
% Inputs:
%   - fun:  function to compute QoI
%   - A, B: samples matrices

M = size(A,1);
p = size(A,2);

% Construct the M x p x p matrix C: 
% 1. make p copies of B in the third dimension (can be done using 'repmat')
% 2. in the i^th copy of B, replace the i^th row by the values from A

%The comment says rows here but shouldn't it be columns? 
C = repmat(B,1,1,p); 
for i = 1:p
   C(:,i,i)=A(:,i); 
end

%calculate the model output for inputs A, B, and C
for i = 1:M
    
    % model outputs
    ya(i) = fun(A(i,:)); %The way I've written these they'll both be row vectors. 
    yb(i) = fun(B(i,:));
    
    for j = 1:p
        yc(i,j) = fun(C(i,:,j));
    end
end

% square of the mean 
f02 = (ya*transpose(yb))/M;

% first order sensitivity indices
for i = 1:p
    Si(i) = (mean(ya .* transpose(yc(:,i)))-f02)/((mean(ya .*ya))-f02);
end

% first order total sensitivity indices
for i = 1:p
    STi(i) = 1 - (mean(yb .* transpose(yc(:,i)))-f02)/((mean(ya .*ya))-f02);
end
