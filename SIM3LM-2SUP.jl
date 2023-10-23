#André Luiz Pianca - Estagiario - Grupo MArÉ - 29/09/2023
#
#


using LinearAlgebra
using Plots
using ControlSystems
include(raw"Material Properties\Properties.jl")



#mini dac

g = 1/1000 #g p/ kg
m_dac = 51.273g
c_dac = Average(4.2,40,Cu_SH)*0.95+Average(4.2,40,Be_SH)*0.05#properties.jl
#c_dac =Average(4.2,40,Cu_SH)*.95+Average(4.2,40,Be_SH)*0.5 
LM1 = m_dac*c_dac

m_s1= 0.332
c_s1= Average(40,300,SS304_SH)#properties.jl
LM2 = m_s1*c_s1
LM3 = LM2





q1 = 0#2.7
q2 = 10e-3
q3 = 0#45
Ta = 4.2
Tb = 40



R1 = 1 #3 casos, para cada alternativa de contato com o cryo
R12 = 1
R23 = 1
R3  = 1

E = [LM1  0    0
     0    LM2  0
     0    0    LM3]

ustat = [q1;q2;q3;Ta;Tb]
#sys = ss(A,B,C,D)
C = I
D = zeros(3,5)

function get_stationary_temperature(R1, R12, R23, R3)


    K = [1/R1+1/R12  -1/R12          0
        -1/R12     +1/R12+1/R23    -1/R23
         0       -1/R23        +1/R23+1/R3]

    L = [1  1   0   1/R1    0
         0  0   0   0       0
         0  0   1   0       1/R3]


    A = -inv(E) * K
    B =  inv(E) * L
 
    xstat = -inv(A)*B*ustat
    Tstat = C*xstat+D*ustat
return Tstat
end

test1 = get_stationary_temperature(1,200,200,1)


R1 = exp10.(LinRange(0,3,3))#dividir em casos
R12 = LinRange(1,200,15)
R23 = LinRange(1,200,30)
R3 = exp10.(LinRange(0,3,3))

function Iterate_stat_temp(R1,R12,R23,R3)




    tLM1 = []
    tLM2 = []
    tLM3 = []





    for i1 in eachindex(R1)
        for i2 in eachindex(R3)
            for i3 in eachindex(R12)
                for i4 in eachindex(R23)
                    t = get_stationary_temperature(R1[i1],R12[i3],R23[i4],R3[i2])
                    push!(tLM1, t[1]) 
                    push!(tLM2, t[2])
                    push!(tLM3, t[3]) 
                end
            end
        end
    end




    tLM1=reshape(tLM1,(length(R23),length(R12),3,3))#(L,C,caser3,caser1)
    tLM2=reshape(tLM2,(length(R23),length(R12),3,3))
    tLM3=reshape(tLM3,(length(R23),length(R12),3,3))
    return tLM1,tLM2,tLM3
end


tLM1,tLM2,tLM3 = Iterate_stat_temp(R1,R12,R23,R3)

axes(tLM1)

tLM1[:,:,1,3]
#cada linha é R23 e cada coluna é R12

for i in 1:1:3
    for j in 1:1:3
        c=contourf(R12,R23,tLM1[:,:,j,i])
        title!("Contour plot for i=$i and j=$j")
        display(c)
    end
end


R1
R3