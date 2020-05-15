fName = "1.txt";

fid = fopen(fName,'r');

input = cell(0,1); % initialize our array that will hold the inputs
while ~feof(fid)
     line = fgetl(fid); % fetch next line
     input{end+1,1} = line; % add line to array
end

% Take input and feed to Horners
[count,~] = size(input);
if (feof(fid)) && count > 2
    hornersMethod(input)
end

fclose(fid);


function [alphas,alpha] = horners(a,n,x0)
    alphas = zeros(n+1,1);
    alphas(n+1,1) = a(n+1); % Array for storing coefficients being synthetically divided
    
    alpha = a(n+1); % initialize our alpha accumulator
    
    for i= n:-1:1
        alpha  = alpha*x0+a(i); % calculate the alpha for that position as per the algorithm
        alphas(i) = alpha; % store alphas in array so we can loop over it later
    end
end

function hornersMethod(input)
    in = sscanf(sprintf('%s ', input{:}), '%f'); % convert string in array to float
    [count,~] = size(in); % Size of the input, helper variable
    n = in(1,1); % First line is the degree n of the polynomial. Add 1 since matlab doesnt have zero indices.
    a = in(2:count-1); % Coefficients are located from second line until second to last line in file
    x0 = in(count); % Last line is the x0 value
    
    % We do the first iteration of horners and obtain the array we will be
    % synthetically dividing further, as well as our first P value.
    [alphaVals,P] = horners(a,n,x0);
    disp(sprintf("P(%.3f) = %.3f",x0,P))
    
    % Repeat synthetic division until we reach the final value using the
    % previous iteration's coefficients. 
    for i = 1:n
       [m,~] = size(alphaVals); % Size of coefficient array
       newVals = alphaVals(2:m,1); % Skip the first value as we do not need it for the current calculation
       
       % Print ith derivative of P
       [alphaVals,newP] = horners(newVals,n-i,x0);
       newP = newP * factorial(i); % Multiply our P based on the factorial of the derivative #
       disp(sprintf("P%d(%.3f) = %.3f",i,x0,newP))
    end

end
