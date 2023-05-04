using LinearAlgebra

E = 21e10 # N/m^2
I = 3055e-8 # m^4
mu = 60 # kg/m
k = 3e6 # N/m^2
F = 100 #kN
l = 0.01 #s

#iΩ = im * Ω

#  Koeffizienten einer gebrochenrationalen Approximation der Steifigkeit
P0 = 1.03377065e7
P1 = 4.8744172e4
P2 = 4.55205728e2
P3 = 1.13950676
P4 = 4.83245447e-3
P5 = 3.11054119e-6

Q0 = 1
Q1 = 4.90648418e-3
Q2 = 2.71147593e-5
Q3 = 5.42912842e-8
Q4 = 6.94958499e-12

###############################################################################
#Aufgabe a)


#Abspaltung

s1_0 = P5
s0_0 = P4
r3_0 = P3
r2_0 = P2
r1_0 = P1
r0_0 = P0


s1_1 = Q4
s0_1 = Q3
r2_1 = Q2
r1_1 = Q1
r0_1 = Q0


s1_2 = r3_0
s0_2 = r2_0
r1_2 = r1_0
r0_2 = r0_0


s1_3 = r2_1
s0_3 = r1_1
r0_3 = r0_1


s0_4 = r0_2/r0_3
s1_4 = r1_2/r0_3


f0 = P0/Q0


A = [s1_0 0 0 0 0;
     0 -s1_1 0 0 0;
     0 0 s1_2 0 0;
     0 0 0 -s1_3 0;
     0 0 0 0 s1_4]

B = [s0_0 1 0 0 0;
     1 -s0_1 -1 0 0;
     0 -1 s0_2 1 0;
     0 0 1 -s0_3 -1;
     0 0 0 -1 s0_4]



r = [f0;0;0;0;0]

println(A)
println(B)
println(r)
#Ausgabe in Latex

###############################################################################
#Aufgabe b)

###############################################################################
#Aufgabe c)



# Define the linear system


# Solve the linear system for ẑ
#z = A \ (r - B * P5)
# println(z)
# Compute the value of K̃(Ω) for a given frequency Ω

# function K̃(Ω)
#     iΩ = im * Ω
#     return P0 + iΩ * (P1 + iΩ * (P2 + iΩ * (P3 + iΩ * (P4 + iΩ * P5 / (1 + iΩ * Q1 + iΩ^2 * Q2 + iΩ^3 * Q3 + iΩ^4 * Q4)))))
# end
F = [100; 0; 0; 0; 0]
t = 0.01
z0 = zeros(5)
z2 = zeros(5)

# z1 = (2*F*t-(0.01/2*B-A)*z0)\(A+0.01/2*B)
# println(z1)

function f(z0,A,B,t)
     return (A+0.01/2*B)^-1 * (2*F*t-(0.01/2*B-A)*z0)
end

for i = 1:100
     z1 = f(z0,A,B,t)
     global z2 = z2 + z1
     global z0 = z1
end

println(z2)
#test
