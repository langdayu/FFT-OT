function [dxx,dxy,dyy,x,y] = FDM(phi,M,N) %Finite differential methods
h_x = 2 / M;
h_y = 2 / N;
x = zeros(M+2,N+2);
y = zeros(M+2,N+2);
Phi = zeros(M+2,N+2); %ghost cells
dxx = zeros(M,N);
dxy = zeros(M,N);
dyy = zeros(M,N);
for i = 2:M+1
    for j = 2:N+1
        x(i,j) = -1 + h_x/2 + (i-2) * h_x;
        y(i,j) = -1 + h_y/2 + (j-2) * h_y;
        Phi(i,j) = phi(i-1,j-1);
        Phi(1,j) = Phi(2,j);
        Phi(M+2,j) = Phi(M+1,j);
    end
    Phi(i,1) = Phi(i,2);
    Phi(i,N+2) = Phi(i,N+1);
end

for i = 2:M+1
    for j = 2:N+1
        dxx(i-1,j-1) = (Phi(i+1,j) + Phi(i-1,j) -2 * Phi(i,j));
        dyy(i-1,j-1) = (Phi(i,j+1) + Phi(i,j-1) -2 * Phi(i,j));
        dxy(i-1,j-1) = (Phi(i+1,j+1) + Phi(i-1,j-1) -2 * Phi(i,j));
    end
end
end

