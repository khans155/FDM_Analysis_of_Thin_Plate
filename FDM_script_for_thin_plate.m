clear;
it=0;
maxit = 7;
temp_neo = zeros(5,5);
m_error = zeros(1,maxit);
while true
    tic;
    it = it+1;
    index = 1;    % Because Matlab indexs from 1, and equations are desribed from 0.
    nodes = 1+(4*2.^(it-1)); a = 0; b = 1; step = (b-a)/(nodes-1); 
    x = a:step:b;  y = a:step:b;  alpha = 1.2;
    dTx0 = 0; dTy0 = 0; Tx1 = y'.*100; Ty1 = x.*100;  %Boundary Conditions 
    neo = zeros(nodes,nodes)+25;  %creating matrix of nodes with WAG of 25
    neo(:,nodes) = Tx1; neo(nodes,:) = Ty1; %applying constant boundary conditions
    iter = 0;
    while true
       k = ((120-neo).^(1.52))./15; 
       iter = iter + 1;
       max_error = 0;
       
       %CALCULATING CONRNER NODE
       for n = 1
           temp=neo(0+index,0+index);   %STORING OLD VALUE 
           neo(0+index,0+index)=(2*k(0+index,1+index)*neo(0+index,1+index)+2*k(1+index,0+index)*neo(1+index,0+index))*1/(2*k(0+index,1+index)+2*k(1+index,0+index));
           neo(0+index,0+index)=alpha*neo(0+index,0+index)+(1-alpha)*temp;  %RELAXATION
           error = abs((temp-neo(0+index,0+index))/(neo(0+index,0+index)));
           if max_error<error
               max_error = error; 
           end
       end
       
       %CALCULATION OF LEFT NODES
       for j = 1:1:(nodes-2)
           temp=neo(j+index,0+index);  %STORING OLD VALUE 
           neo(j+index,0+index)=(2*k(j+index,1+index)*neo(j+index,1+index)+k(j+1+index,0+index)*neo(j+1+index,0+index)+k(j-1+index,0+index)*neo(j-1+index,0+index))*1/(2*k(j+index,1+index)+k(j+1+index,0+index)+k(j-1+index,0+index));
           neo(j+index,0+index)=alpha*neo(j+index,0+index)+(1-alpha)*temp;   %RELAXATION
           error = abs((temp-neo(j+index,0+index))/(neo(j+index,0+index)));
           if max_error<error
           max_error = error; 
           end
       end
       
       %CALCULATING LOWER NODES
       for i = 1:1:(nodes-2)
           temp=neo(0+index,i+index);  %STORING OLD VALUE 
           neo(0+index,i+index)=(k(0+index,i+1+index)*neo(0+index,i+1+index)+k(0+index,i-1+index)*neo(0+index,i-1+index)+2*k(1+index,i+index)*neo(1+index,i+index))*1/(k(0+index,i+1+index)+k(0+index,i-1+index)+2*k(1+index,i+index));
           neo(0+index,i+index)=alpha*neo(0+index,i+index)+(1-alpha)*temp;  %RELAXATION
           error = abs((temp-neo(0+index,i+index))/(neo(0+index,i+index)));
           if max_error<error
               max_error = error; 
           end
       end
       
       %CALCUATION OF INTERIOR NODES
       for i = 1:1:(nodes-2)
          for j = 1:1:(nodes-2)
              temp=neo(j+index,i+index);  %STORING OLD VALUE 
              neo(j+index,i+index)=(k(j+index,i+1+index)*neo(j+index,i+1+index)+k(j+index,i-1+index)*neo(j+index,i-1+index)+k(j+1+index,i+index)*neo(j+1+index,i+index)+k(j-1+index,i+index)*neo(j-1+index,i+index))*1/(k(j+index,i+1+index)+k(j+index,i-1+index)+k(j+1+index,i+index)+k(j-1+index,i+index));
              neo(j+index,i+index)=alpha*neo(j+index,i+index)+(1-alpha)*temp;   %RELAXATION
              error = abs((temp-neo(j+index,i+index))/(neo(j+index,i+index)));
              if max_error<error
                  max_error = error; 
              end
          end
       end
       
       if iter == 10000
           break
       end
    end
    
    m_error(it) = max_error;
    it
    toc
    if it == maxit
        break
    end
end
contour(x,y,neo,200)
title('Contour of T(x,y)')
xlabel('x (m)')
ylabel('y (m)')