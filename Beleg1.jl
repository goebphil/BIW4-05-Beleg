
using Plots

# Koeffizienten einer gebrochenrationalen Approximation der Steifigkeit
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

# Abspaltung
s1_0 = P5 / Q4
s0_0 = (P4 - s1_0 * Q3) / Q4
r3_0 = P3 - s0_0 * Q3 - s1_0 * Q2
r2_0 = P2 - s0_0 * Q2 - s1_0 * Q1
r1_0 = P1 - s0_0 * Q1 - s1_0
r0_0 = P0 - s0_0

s1_1 = Q4 / r3_0
s0_1 = (Q3 - s1_1 * r2_0) / r3_0
r2_1 = Q2 - s0_1 * r2_0 - s1_1 * r1_0
r1_1 = Q1 - s0_1 * r1_0 - s1_1 * r0_0
r0_1 = Q0 - s0_1 * r0_0

s1_2 = r3_0 / r2_1
s0_2 = (r2_0 - s1_2 * r1_1) / r2_1
r1_2 = r1_0 - s0_2 * r1_1 - s1_2 * r0_1
r0_2 = r0_0 - s0_2 * r0_1

s1_3 = r2_1 / r1_2
s0_3 = (r1_1 - s1_3 * r0_2) / r1_2
r0_3 = r0_1 - s0_3 * r0_2

s0_4 = r0_2 / r0_3
s1_4 = r1_2 / r0_3

f0 = 0

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

r = [0; 0; 0; 0; 0]

println(A)
println(B)
println(r)

# Zeitschrittlänge wird hier im Skript mit t anstatt delta t beschrieben
delta_t = 0.00001 # s
t = 0
t_end = 1
k = 1

# Array für die Darstellung des Plots
z3 = Float64[]
z0 = [0;0;0;0;0]

function Belastung(t)
    if t < 0 && t < 0.01
        return t * [100000; 0; 0; 0; 0]
    elseif t <= 0.01 && t < 0.02
        return [100000; 0; 0; 0; 0]
    elseif t <= 0.02 && t < 0.03
        return -t * [100000; 0; 0; 0; 0]
    else
        return [0;0;0;0;0]
    end
end

while true
    z1 = inv(A + delta_t/2*B) * ((Belastung(t) + Belastung(t+delta_t))/2*delta_t - (delta_t/2*B-A)*z0)


    global z0 = z1

    push!(z3, z1[1])
    global k += 1
    global t += delta_t
    if t >= t_end
        break
    end
end

println("Lösung für z1")
println(z0)

plot(z3, xlabel="Sekunden", ylabel="Werte für z1")
savefig("plot.png")
