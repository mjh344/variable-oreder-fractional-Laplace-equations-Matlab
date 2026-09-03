
function A = coef_matrix(N, alpha_fun)

    if nargin < 2
        %alpha_fun = @(x) 0.9 - 0.1*x.^2;
        alpha_fun = @(x) 0.9 - 0.1*sin(abs(pi*x)); 
        %alpha_fun = @(x) 0.9 - 0.3*abs(x).^3;
        %alpha_fun = @(x) 0.9 - 0.1*abs(x).^3;


       
    end
    a = -1;
    b = 1;
    x = graded_mesh(a, b, N, alpha_fun);   
    h = x(2:end) - x(1:end-1);
    alp = alpha_fun(x);
    d = 1./alp./(1-alp);

    
    aii = (x(2:end-1)-x(1:end-2)).^(1-alp(2:end-1))./h(1:end-1) + ...
          (x(3:end)-x(2:end-1)).^(1-alp(2:end-1))./h(2:end);
    aii1 = -(h(3:end)+h(2:end-1)).*(x(3:end-1)-x(2:end-2)).^(1-alp(2:end-2))./h(2:end-1)./h(3:end) + ...
           (x(4:end)-x(2:end-2)).^(1-alp(2:end-2))./h(3:end);
    aii_1 = (x(3:end-1)-x(1:end-3)).^(1-alp(3:end-1))./h(1:end-2) - ...
            (h(1:end-2)+h(2:end-1)).*(x(3:end-1)-x(2:end-2)).^(1-alp(3:end-1))./h(1:end-2)./h(2:end-1);

    temp_r = cell(N-1,1);
    temp_l = cell(N-1,1);
    ai = cell(N-1,1);
    ak = cell(N-1,1);

    for i = 1:N-3
        temp_r{i} = (x(i+2:end-2)-x(i+1)).^(1-alp(i+1))./h(i+2:end-1) - ...
                    (h(i+3:end)+h(i+2:end-1)).*(x(i+3:end-1)-x(i+1)).^(1-alp(i+1))./(h(i+2:end-1).*h(i+3:end)) + ...
                    (x(i+4:end)-x(i+1)).^(1-alp(i+1))./h(i+3:end);
        ai{i} = [zeros(1,i+1), temp_r{i}];
    end

    for k = 3:N-1
        temp_l{k} = (x(k+1)-x(1:k-2)).^(1-alp(k+1))./h(1:k-2) - ...
                    (h(2:k-1)+h(1:k-2)).*(x(k+1)-x(2:k-1)).^(1-alp(k+1))./(h(1:k-2).*h(2:k-1)) + ...
                    (x(k+1)-x(3:k)).^(1-alp(k+1))./h(2:k-1);
        ak{k} = [temp_l{k}, zeros(1,N-k+1)];
    end

    A1 = diag(aii) + diag(aii_1,-1) + diag(aii1,1) + ...
         vertcat(ai{1:N-3}, zeros(2,N-1)) + vertcat(zeros(2,N-1), ak{3:N-1});
    A1 = d(2:end-1)'.*A1;

    A = [zeros(N-1,1), A1, zeros(N-1,1)];
    A = [[1, zeros(1,N)]; A; [zeros(1,N), 1]];
end
