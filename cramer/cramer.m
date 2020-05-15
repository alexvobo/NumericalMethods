fName = "1.txt";

fid = fopen(fName,'r');

input = cell(0,1); % initialize our array that will hold the inputs
while ~feof(fid)
     line = fgetl(fid); % fetch next line
     input{end+1,1} = line; % add line to array
end

% Take input and feed to Cramer
[count,~] = size(input);
if (feof(fid)) && count > 2
    Cramer(input);
end

fclose(fid);


function sol = Cramer(input)
    in = sscanf(sprintf('%s ', input{:}), '%f'); % convert string in array to float
    [count,~] = size(in); % size of input
    n = in(1,1); % number of columns in matrix
    
    aTemp = in(2:count-n); % create temp array of elements in matrix A before reshaping
    bTemp = in((count-n+1):count); % create temp array of elements in matrix B before reshaping
    [bTempSize,~] = size(bTemp); % counts elements in solution part
    a = reshape(aTemp,n,n)';
    b = reshape(bTemp,n,bTempSize/n);
    
    [rowB,colB] = size(b);
    a = [a,b];
    sign = 1; % Set sign variable to keep track of row swapping
    % For all columns below the main diagonal
    %det = 1; % Keep track of potential row swaps, will change det
    for j = 1:(n-1)
        %%%%% Find the row with the largest pivot for that column %%%%
        max = abs(a(j,j)); % prime max with first value on the diagonal we are under
        k = j;
        for r = (j+1):n % look through remaining rows under max to find a greater value
            if abs(a(r,j)) > max
                k = r;
                max = abs(a(r,j));
%                 disp(max);
                sign  = sign * -1; % Change sign because we swapped rows
            end
        end
        a = swapRow(a,j,k);
        
        %%%%% Check if theres no valid value left after row swap %%%%
        if (a(j,j) == 0)
            disp("Not Unique Solution")
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Do row operations  %%%%%%%%%%%%%%%%%%%
        for i = (j+1):n
                scalar = a(i,j)/max;
                scalarRow = multiplyRow(a,j,scalar);
                a = subtractRow(a,i,scalarRow);
        end
    end

    det = a(1,1); %initiate determinant with first element
    for i = 2:n
        det = det*a(i,i);
    end
    det = det*sign; % Need to account for sign changes based on row swapping
    disp("determinant A = " +  det)
    %%% Back sub, start from last row. %%%
    
    % Create column vector for solutions, size of B
    solutions = zeros(rowB,colB);
    
    % Size of Ax=b matrix in RREF
    [~,colA] = size(a);
    
    % Need to augment matrix we used Gaussian on, need access to right hand side
    for col = n+1 : colA 
        
        % Create a solution vector of size (n,1)
        xVars = zeros(n,1); 
        
        % Create an augmented matrix to perform backsub on
        aug = [a(1:n,1:n) a(1:n,col)];
        
        % Calc last row of augmented matrix by division (solves for Xn)
        xVars(n,1) = aug(n,n+1)/aug(n,n);
        
        % Perform algebra on the rest of the rows (solves for X1...Xn-1)
        for i=(n-1):-1:1
            sum = 0;
            for k = (i+1):n
                sum = sum + (aug(i,k).*xVars(k,1));
            end
            xVars(i,1) = (aug(i,n+1) - sum)/aug(i,i);
        end
        
        % Replace column of 0s in placeholder by the solution vector
        solutions(:,col-n) = xVars(:,1);
    end
    

    [rowSol,colSol] = size(solutions);
    
    % Print out determinants
    for col=1:colSol
         for row=1:rowSol
            disp("determinant A"+ row + " = " + solutions(row,col)*det);
         end
    end
    % Print out the solution vector(s)
    for col=1:colSol
         for row=1:rowSol
            disp("x"+ row + " = " + solutions(row,col));
         end
    end
   
end

function A = swapRow(A, row1, row2)
  % Swap row 1  with row 2, return a new matrix
  A([row1, row2], :) = A([row2, row1], :);
end
function A = multiplyRow(A,row1,scalar)
   % Multiplies entire row by a scalar, returns the row 
  A = A(row1,:).*scalar;
end
function A = subtractRow(A,row1,scalarRow)
  % Subtracts row we obtain in multiplyRow from a row1 we specify.
  A(row1,:) = A(row1,:) - scalarRow;
end