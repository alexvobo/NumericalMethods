fName = "1.txt";
fid = fopen(fName,'r');
fxn  = str2func("@(x)" + fgetl(fid)); % function to integrate
a = str2double(fgetl(fid)); % integrate from a
b = str2double(fgetl(fid)); % integrate to b
n = str2double(fgetl(fid)); % subintervals
simpsonsRule(fxn,a,b,n);
fclose(fid);

function approx = simpsonsRule(fxn,a,b,n)
    h = (b-a)/n; % to use for our increment operations
    odds = 0 ;
    evens = 0;
    
    for i = 1:n-1     % Go through intervals 1-(n-1). We compute interval 0 and interval n later.   
                
        xi = a+h*(i); % compute xi to plug into our fxn

        if mod(i,2) == 0
            evens = evens + fxn(xi); % sum when i is even
        else
            
            odds = odds+fxn(xi); % sum when i is odd
        end
    end
    
    approx = (h/3)*(fxn(a) + 2*evens + 4*odds + fxn(b)) % Formula for simpsons rule
end

