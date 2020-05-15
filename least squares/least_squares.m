fName = "1.txt";
fid = fopen(fName,'r');
n = str2double(fgetl(fid)); % num pairs
degree = str2double(fgetl(fid));  % degree
x = zeros(n,1);
y = zeros(n,1);
for i = 1:n
    x(i) = str2double(fgetl(fid)); % x value
    y(i) = str2double(fgetl(fid)); % y value
end
pairs = [x y];

leastSquares(n,x,y,degree)

fclose(fid);


function approx = leastSquares(n,x,y,degree)
  % Calculate our matrix of sums
  M = zeros(degree+1);
  M(1,1) = n;
  for i = 1:degree+1
      for j = 1:degree+1
          if (i == 1 && j == 1)
              M(i,j) = n;
          else
              M(i,j) = sum(x.^(i+j-2));
          end
      end
  end
  M;
  
  % Calculate vector of sums
  b = zeros(degree+1,1);
  for i = 0:degree
      if i == 0
          b(i+1) =  sum(y);
      else
          b(i+1) = sum(y.*(x.^i));
      end
  end
  b;
  
  % Solution vector with our coefficients
  a = inv(M)*b;
  
  xx = [0:0.1:n+10];
    
  % Generate the function for our approximation
  yy = a(1);
  for i = 2:degree+1
      yy = yy + a(i)*(xx.^(i-1));
  end
  yy;
  
  disp(a') % coefficients of the polynomial
  plot(xx,yy,x,y,'*'); % plot the splines
end