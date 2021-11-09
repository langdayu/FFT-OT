% main function
clear
phi = zeros(600);
src = phi;
epsilon = 1e-4;
Stopcriterion = 1;
n = 0;
while Stopcriterion > epsilon
    n=n+1;
    new_phi = Possolver(phi,src);
    Stopcriterion = norm(phi - new_phi);
    phi = new_phi;
end 
T = eye(600,600)+gradient(phi);

