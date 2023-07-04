


r0 = 1.5
w = 2*pi*40

G0 = 80
roh0 = 1500
v0 = 0.25

G1 = 0.5*G0
roh1 = roh0
v1 = 0.3

G2 = 0.2*G0
roh2 = 0.89*roh0
v2 = 1/3

D = 0.05

#println(4*G0*r0/(1-v0))
#println(8*G0*r0^3/(1-v0)/3)

cs1 = sqrt(G1/roh1)
println("cs1: ",cs1)

a0 = w*r0/cs1/1000
println(a0)


function Sv()
    4*G0*r0/(1-v0)*(kv+i*a0*cv)
end
