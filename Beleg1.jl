using LinearAlgebra
using QuadGK


#  Koeffizienten einer gebrochenrationalen Approximation der Steifigkeit
P0 = 10337706.5
P1 = 48744.172
P2 = 455.205728
P3 = 1.13950676
P4 = 0.00483245447
P5 = 3.11054119e-6

Q0 = 1
Q1 = 0.00490648418
Q2 = 2.71147593e-5
Q3 = 5.42912842e-8
Q4 = 6.94958499e-12

###############################################################################



#Abspaltung

s1_0 = P5/Q4
s0_0 = (P4-s1_0*Q3)/Q4
r3_0 = P3-s0_0*Q3-s1_0*Q2
r2_0 = P2-s0_0*Q2-s1_0*Q1
r1_0 = P1-s0_0*Q1-s1_0
r0_0 = P0-s0_0


s1_1 = Q4/r3_0
s0_1 = (Q3-s1_1*r2_0)/r3_0
r2_1 = Q2-s0_1*r2_0-s1_1*r1_0
r1_1 = Q1-s0_1*r1_0-s1_1*r0_0
r0_1 = Q0-s0_1*r0_0


s1_2 = r3_0/r2_1
s0_2 = (r2_0-s1_2*r1_1)/r2_1
r1_2 = r1_0-s0_2*r1_1-s1_2*r0_1
r0_2 = r0_0-s0_2*r0_1


s1_3 = r2_1/r1_2
s0_3 = (r1_1-s1_3*r0_2)/r1_2
r0_3 = r0_1-s0_3*r0_2


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



###############################################################################


#Zeitschrittlänge
t = 0.000001 #s

#Anzahl der Iteratiosschritte
j_ges = Int(1/t+1)
j_Schritt = Int(1/t)
#Befüllen der Anfangsarrays
z0 = [0;0;0;0;0]
z2 = [0;0;0;0;0]


# Erstellen eines Arrays über die Zeitschritte, für in die Integration
j = zeros(j_ges, 1)


for i in 1:j_ges
    j[i, 1] = (i - 1) * t
end



# Funktionen von F(t)

function F1(k)
     return [1;0;0;0;0]*100000/0.01*k
end

F2 = [1;0;0;0;0]*100000

function F3(k)
     return [1;0;0;0;0]*(100000 - 100000/0.01*k)
end



for i =1:j_ges

     if i == 1
          z1 = (A+0.001/2*B)^-1 * (F1(j[i])+F1(j[i+1]))/2*t
          global z2 = z2 + z1
          global z0 = z1

     end

     if i>1 && i<=0.1*j_Schritt
          z1 = inv((A+t/2*B)) * ((F1(j[i])+F1(j[i+1]))/2*t-(t/2*B-A)*z0)
          global z2 = z2 + z1
          global z0 = z1

     end

     if i>0.1*j_Schritt && i<=0.2*j_Schritt
          z1 = inv((A+t/2*B)) * (F2*t-(t/2*B-A)*z0)
          global z2 = z2 + z1
          global z0 = z1


     end

     if i>0.2*j_Schritt && i<=0.3*j_Schritt
          z1 = inv((A+t/2*B)) * ((F3(j[i])+F3(j[i+1]))/2*t-(t/2*B-A)*z0)
          global z2 = z2 + z1
          global z0 = z1

     end

     if i>0.3*j_Schritt
          z1 = inv(A+t/2*B) * -(t/2*B-A)*z0
          global z2 = z2 + z1
          global z0 = z1
     end

end

println("Lösung für z1")
println(z0)

