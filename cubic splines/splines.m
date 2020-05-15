fName = "1.txt";
fid = fopen(fName,'r');

n = str2double(fgetl(fid)); % num pairs
x = zeros(n,1);
y = zeros(n,1);

for i = 1:n
    x(i) = str2double(fgetl(fid)); % x
    y(i) = str2double(fgetl(fid)); % y
end
% pairs = [x y]; 

c = solveForC(x,y,n);
d = solveForD(x,c,n);
b = solveForB(x,y,c,n);
plotSplines(x,y,c,d,b,n);

splineEquations(x,y,c,d,b,n)

fclose(fid);

function splineEquations(x,y,c,d,b,n)
    % formatString = "S(%d) = %f + %f(x-%f) + %f(x-%f)^2 + %d(x-%f)^3";
    formatString = "%f %f %f %f";
    
    % Print out the coefficients of each function
    for i = 1:n-1
       disp(sprintf(formatString,y(i),b(i),c(i),d(i))) 
    end
end
function plotSplines(x,y,c,d,b,n)
    % Plot each piecewise function
    for i=1:n-1
          fplot(@(k) y(i) + b(i)*(k-x(i))+c(i)*(k-x(i)).^2+d(i)*(k-x(i)).^3,[x(i),x(i+1)]);
          hold on 
    end
end
function c = solveForC(x,y,n)
    % find the Cs 
    h = zeros(n - 1,1);
    b = zeros(n,1);
    A = eye(n);
    
    % Find H values - Width of interval
    for i = 1:n-1
         h(i) = (x(i + 1) - x(i));
    end
    
    % Create our matrix to solve our vector of Cs according to the
    % algorithm
    for i = 2:n-1
        A(i,i-1) = h(i-1);
        A(i,i) = 2*(h(i-1)+h(i));
        A(i,i+1) = h(i);
        b(i) = 3 * (((y(i + 1) - y(i)) / h(i)) - ((y(i) - y(i - 1)) / h(i - 1)));
    end
    
    b = reshape(b, [n, 1]);
    
    c = A\b;
end
function d = solveForD(x,c,n)
    % find the Ds
    d = zeros(n - 1,1);
    % Follow the algorithm, last row of vector is 0
    for i = 1:n-1
        h = (x(i + 1) - x(i));
        d(i) = ((c(i + 1) - c(i)) / (3 * h));
    end
end
function b = solveForB(x,y,c,n)
  % find the Bs 
  b = zeros(n-1,1);
  % Follow the algorithm, last row of vector is 0
  for i=1:n-1
      h = (x(i+1) - x(i));
      b(i) = ((y(i + 1) - y(i)) / h) - h * (2 * c(i) + c(i + 1)) / 3;
  end
  
end