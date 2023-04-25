using LinearAlgebra

E = 21e10 # N/m^2
I = 3055e-8 # m^4
mu = 60 # kg/m
k = 3e6 # N/m^2
F = 100 #kN
l = 0.01 #s



# Define the given coefficients
P0 = 1.03377065e7
P1 = 4.8744172e4
P2 = 4.55205728e2
P3 = 1.13950676
P4 = 4.83245447e-3
P5 = 3.11054119e-6

Q1 = 4.90648418e-3
Q2 = 2.71147593e-5
Q3 = 5.42912842e-8
Q4 = 6.94958499e-12

# Define the linear system
A = [Q1 Q2 Q3 Q4 0; 1 Q1 Q2 Q3 Q4; 0 1 Q1 Q2 Q3; 0 0 1 Q1 Q2; 0 0 0 1 Q1]
B = [1; 0; 0; 0; 0]
r = [P0; P1; P2; P3; P4]

# Solve the linear system for ẑ
z = A \ (r - B * P5)
println(z)
# Compute the value of K̃(Ω) for a given frequency Ω
function K̃(Ω)
    iΩ = im * Ω
    return P0 + iΩ * (P1 + iΩ * (P2 + iΩ * (P3 + iΩ * (P4 + iΩ * P5 / (1 + iΩ * Q1 + iΩ^2 * Q2 + iΩ^3 * Q3 + iΩ^4 * Q4)))))
end



