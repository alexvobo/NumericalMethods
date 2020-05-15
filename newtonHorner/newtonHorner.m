fName = "1.txt";

fid = fopen(fName,'r');

input = cell(0,1); % initialize our array that will hold the inputs
while ~feof(fid)
     line = fgetl(fid); % fetch next line
     input{end+1,1} = line; % add line to array
end

% Take input and feed to newtonHorners
[count,~] = size(input);
if (feof(fid)) && count > 2
    newtonHorners(input)
end

fclose(fid);

function root = newtonHorners(input)
    in = sscanf(sprintf('%s ', input{:}), '%f'); % convert string in array to float
    [count,~] = size(in);
    n = in(1,1); % First line is the degree n of the polynomial. Add 1 since matlab doesnt have zero indices.
    a = in(2:count-3); % Coefficients are located from second line until second to last line in file
    x0 = in(count-2); % initial value
    errorBound = in(count-1); % error tolerance 
    N = in(count); % max iterations
    
    currIterations = 1;

    % Ensure we have not hit the threshold for iterations
    while(currIterations <= N)
        alpha = a(n+1); % initialize alpha value
        beta = a(n+1); % initialize beta value
        for i=n:-1:1
            alpha = alpha*x0+a(i); % This will be our f(x)
            if i > 1
               beta = beta*x0+alpha; % This will be our f'(x)
            end
        end
        
        % Calculate x1 using newtons,
        % but f(x) is alpha and f'(x) is beta
        x1 = x0-(alpha/beta);
        e = (abs(x1-x0));
    
        % If the error of our approximation is smaller than the threshold,
        % stop and display solution
        if(e <= errorBound)
             % currIterations;
             % disp("Solution: " + x1);
              disp(x1)
            break;
        end
        % If the solution is not found our x0 for the next iteration
        % becomes x1
        x0 = x1;
        currIterations = currIterations+1;
    end
    % If we have not found the solution after N iterations stop
    if currIterations > N
        disp("No Solution");
    end
end