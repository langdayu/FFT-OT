%possion equation solver: Delta phi = rhs ;
% input : phi_n
% output: U,i.e. phi_(n+1) 
function U = Possolver(phi,src) 
[M,N] = size(src);
% h_x = 2 / M;
% h_y = 2 / N;
% Rhs = zeros(M,N);
[dxx,dxy,dyy,x,y] = FDM(phi,M,N);
u1 = 0.1;
u2 = 0.1;
sigma1 = 0.3;
sigma2 = 0.3;
rou = 0.1;
X = x(2:M+1,2:N+1);
Y = y(2:M+1,2:N+1);
% f: 2d normal distribution density function
f = 1/(2*pi*sigma1*sigma2*sqrt(1-rou*rou)).*exp(-1/(2*(1-rou^2)).*((X-u1).*(X-u1)/(sigma1*sigma1)-2*rou*(X-u1).*(Y-u2)/(sigma1*sigma2)+(Y-u2).*(Y-u2)/(sigma2*sigma2)));
f_dens = f*M*N/sum(sum(f));
% surf(X,Y,f);
% shading interp;
% colorbar;
for i = 1:M
    for j = 1:N
    rhs = (dxx(i,j)+1)^2+(dyy(i,j)+1)^2+2*dxy(i,j)^2+f_dens(i,j);
    Rhs(i,j) = sqrt(rhs)-2;
    end
end

F = dct2(Rhs);          %DCT变换
for i = 1:M
    for j = 1:N
         den = 8*(cos(i*pi/M)-1) + 8*(cos(j*pi/N)-1);
         u(i,j) = F(i,j) / den;
    end
end
u(1,1) = 0;
U = idct2(u);

 