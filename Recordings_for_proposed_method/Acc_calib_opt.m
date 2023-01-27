function [M,B] = Acc_calib_opt(data)
%   Inputs:
%   data          (N,3)    columns ax, ay and az acceleration average
%                          values obtain in m/sec^2
%                          N represent N different static positions.
%
%   Outputs:
%   M        dims: [3,3]     Misalignment and scale Matrix
%   B        dims: [3,1]     Bias Vector

% Configurable variables

lambda = 1;      % Damping Gain - Start with 1
kl = 0.01;       % Damping paremeter - has to be less than 1. Changing this will affect rate of convergence ; Recommend to use k1 between 0.01 - 0.05
tol = 1e-9;      % Convergence criterion threshold
Rold = 100000;   % Better to leave this No. big.
itr = 2000;       % No. Of iterations. If your solutions don't converge then try increasing this. Typically it should converge within 20 iterations


% Initial Guess values of M and B.  Change this only if you need to

Mxx0 = 1;
Mxy0 = 0;
Mxz0 = 0;
Myy0 = 1;
Myz0 = 0;
Mzz0 = 1;

Bx0 = 0;
By0 = 0;
Bz0 = 0;

%   Actual Algorithm
g = 9.79;
V = data/g;
[r, c] = size(V);

if r < 9
    disp('Need atleast 9 Measurements for the calibration procedure!')
    return
end

if c ~= 3
    disp('Not enough columns in the data')
    return
end


% f is the error function given by Ax^2+Ay^2+Az^2 - g, where g = 1
f = @(Vx, Vy, Vz, Mxx,Mxy,Mxz,Myy,Myz,Mzz,Bx,By,Bz) (Mxx*(Bx - Vx) + Mxy*(By - Vy) + Mxz*(Bz - Vz))^2 + (Mxy*(Bx - Vx) + Myy*(By - Vy) + Myz*(Bz - Vz))^2 + (Mxz*(Bx - Vx) + Myz*(By - Vy) + Mzz*(Bz - Vz))^2-1;

% Functions f1 to f9 are the elements of the Jacobian vector (partial
% derivatives of the error function with respect to the gain and bias
% components)
f1 = @(Vx, Vy, Vz, Mxx,Mxy,Mxz,Myy,Bx,By,Bz) (Bx - Vx)*(Mxx*(Bx - Vx) + Mxy*(By - Vy) + Mxz*(Bz - Vz));
f2 = @(Vx, Vy, Vz, Mxx,Mxy,Mxz,Myy,Myz,Bx,By,Bz) 2*(By - Vy)*(Mxx*(Bx - Vx) + Mxy*(By - Vy) + Mxz*(Bz - Vz)) + 2*(Bx - Vx)*(Mxy*(Bx - Vx) + Myy*(By - Vy) + Myz*(Bz - Vz));
f3 = @(Vx,Vy,Vz,Mxx,Mxy,Mxz,Myz,Mzz,Bx,By,Bz) 2*(Bx - Vx)*(Mxz*(Bx - Vx) + Myz*(By - Vy) + Mzz*(Bz - Vz)) + 2*(Bz - Vz)*(Mxx*(Bx - Vx) + Mxy*(By - Vy) + Mxz*(Bz - Vz));
f4 = @(Vx,Vy,Vz,Mxy,Myy,Myz,Bx,By,Bz) 2*(By - Vy)*(Mxy*(Bx - Vx) + Myy*(By - Vy) + Myz*(Bz - Vz));
f5=  @(Vx,Vy,Vz,Mxy,Mxz,Myy,Myz,Mzz,Bx,By,Bz) 2*(By - Vy)*(Mxz*(Bx - Vx) + Myz*(By - Vy) + Mzz*(Bz - Vz)) + 2*(Bz - Vz)*(Mxy*(Bx - Vx) + Myy*(By - Vy) + Myz*(Bz - Vz));
f6=  @(Vx,Vy,Vz,Mxz,Myz,Mzz,Bx,By,Bz) 2*(Bz - Vz)*(Mxz*(Bx - Vx) + Myz*(By - Vy) + Mzz*(Bz - Vz));                                    
f7 = @(Vx,Vy,Vz,Mxx,Mxy,Mxz,Myy,Myz,Mzz,Bx,By,Bz) 2*Mxx*(Mxx*(Bx - Vx) + Mxy*(By - Vy) + Mxz*(Bz - Vz)) + 2*Mxy*(Mxy*(Bx - Vx) + Myy*(By - Vy) + Myz*(Bz - Vz)) + 2*Mxz*(Mxz*(Bx - Vx) + Myz*(By - Vy) + Mzz*(Bz - Vz));
f8 = @(Vx,Vy,Vz,Mxx,Mxy,Mxz,Myy,Myz,Mzz,Bx,By,Bz) 2*Mxy*(Mxx*(Bx - Vx) + Mxy*(By - Vy) + Mxz*(Bz - Vz)) + 2*Myy*(Mxy*(Bx - Vx) + Myy*(By - Vy) + Myz*(Bz - Vz)) + 2*Myz*(Mxz*(Bx - Vx) + Myz*(By - Vy) + Mzz*(Bz - Vz));
f9 = @(Vx,Vy,Vz,Mxx,Mxy,Mxz,Myy,Myz,Mzz,Bx,By,Bz) 2*Mxz*(Mxx*(Bx - Vx) + Mxy*(By - Vy) + Mxz*(Bz - Vz)) + 2*Myz*(Mxy*(Bx - Vx) + Myy*(By - Vy) + Myz*(Bz - Vz)) + 2*Mzz*(Mxz*(Bx - Vx) + Myz*(By - Vy) + Mzz*(Bz - Vz));


Vx = V(:,1);
Vy = V(:,2);
Vz = V(:,3);


m = length(V);
R = zeros(m, 1);
J = zeros(m, 2);
v = [Mxx0, Mxy0, Mxz0, Myy0, Myz0, Mzz0, Bx0, By0, Bz0]';

for n=0:itr % iterate
    % Calculate the Jacobian at every iteration
    for i=1:length(Vx)
        R(i)    =  f (Vx(i),Vy(i),Vz(i),Mxx0,Mxy0,Mxz0,Myy0,Myz0,Mzz0,Bx0,By0,Bz0);
        J(i, 1) =  f1(Vx(i),Vy(i),Vz(i),Mxx0,Mxy0,Mxz0,Myy0,Bx0,By0,Bz0);
        J(i, 2) =  f2(Vx(i),Vy(i),Vz(i),Mxx0,Mxy0,Mxz0,Myy0,Myz0,Bx0,By0,Bz0);
        J(i, 3) =  f3(Vx(i),Vy(i),Vz(i),Mxx0,Mxy0,Mxz0,Myz0,Mzz0,Bx0,By0,Bz0);
        J(i, 4) =  f4(Vx(i),Vy(i),Vz(i),Mxy0,Myy0,Myz0,Bx0,By0,Bz0);
        J(i, 5) =  f5(Vx(i),Vy(i),Vz(i),Mxy0,Mxz0,Myy0,Myz0,Mzz0,Bx0,By0,Bz0);
        J(i, 6) =  f6(Vx(i),Vy(i),Vz(i),Mxz0,Myz0,Mzz0,Bx0,By0,Bz0);
        J(i, 7) =  f7(Vx(i),Vy(i),Vz(i),Mxx0,Mxy0,Mxz0,Myy0,Myz0,Mzz0,Bx0,By0,Bz0);
        J(i, 8) =  f8(Vx(i),Vy(i),Vz(i),Mxx0,Mxy0,Mxz0,Myy0,Myz0,Mzz0,Bx0,By0,Bz0);
        J(i, 9) =  f9(Vx(i),Vy(i),Vz(i),Mxx0,Mxy0,Mxz0,Myy0,Myz0,Mzz0,Bx0,By0,Bz0);
    end

    Rnew = norm(R);

    H = inv(J'*J); % Hessian matrix
    D = (J'*R)';

    v = v - lambda*(D*H)';
    % This is to make sure that the error is decereasing with every
    % iteration
    if (Rnew <= Rold)
        lambda = lambda-kl*lambda;
    else
        lambda = kl*lambda;
    end

    % Iterations are stopped when the following convergence criteria is
    % satisfied
    if  (n>1)
        disp(sprintf('%d', abs(max(2*(v-vold)/(v+vold)))));
        if (abs(max(2*(v-vold)/(v+vold))) <= tol)
            disp('Convergence achieved');
            break;
        end
    end

    Mxx0 = v(1);
    Mxy0 = v(2);
    Mxz0 = v(3);
    Myy0 = v(4);
    Myz0 = v(5);
    Mzz0 = v(6);
    Bx0 = v(7);
    By0 = v(8);
    Bz0 = v(9);
    vold = v;
    Rold = Rnew;
end

% Save Outputs

M = [v(1) v(2) v(3); v(2) v(4) v(5); v(3) v(5) v(6)];
B = [v(7);v(8);v(9)];
save M.mat;
save B.mat;
