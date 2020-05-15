fName = "1.txt";

fid = fopen(fName,'r');

input = cell(0,1); % initialize our array that will hold the inputs
while ~feof(fid)
     line = fgetl(fid); % fetch next line
     input{end+1,1} = line; % add line to array
end

% Take input and feed to nevilles
[count,~] = size(input);
if (feof(fid)) && count > 2
    nevilles(input)
end

fclose(fid);


function sol = nevilles(input)
    in = sscanf(sprintf('%s ', input{:}), '%f'); % convert string in array to float
    [count,~] = size(in);
    n = in(1,1);
    a = in(2:count-1); % pairs
    x = in(count); % initial value
    n=n+1; % matlab does not like 0 indices
    
    points = zeros(n,2); % table of (x,y) values
    c = 1; % separate counter to create table of points
    for i=1:2:size(a)
        points(c,1) = a(i);
        points(c,2) = a(i+1);
        c=c+1;
    end
%     disp(points);

    P = eye(n).*points(:,2); % Solution matrix, set diagonal to y values
    for d=1:(n)
        for i=1:(n-d)
            j=i+d;
            xi = points(i,1);
            xj = points(j,1);
            P(i,j) = ((x-xi)*P(i+1,j)-(x-xj)*P(i,j-1))/(xj-xi);
        end
    end

    % Our solution is the last column in the first row
    disp("P("+x+")= "+P(1,n))
end

