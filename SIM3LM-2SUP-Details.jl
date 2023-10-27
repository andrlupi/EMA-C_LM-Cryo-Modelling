#André Luiz Pianca - Estagiario - Grupo MArÉ - 23/10/2023
#
#

using .Threads
using ControlSystems
using LinearAlgebra
using Plots
include(raw"Material Properties\Properties.jl")



const q1::Float64 = -1
const q2::Float64 = 10e-3
const q3::Float64 = -45
const Ta::Float64 = 4.2
const Tb::Float64 = 40
const size::Int32 = 1000

#definir cases p/ LM1 e LM5
#
#case LM1 = LM5

const M1::Float64 = 134.603e-3

c1(T1) = Cu_SH(T1)

const M5::Float64 = 1*M1

c5(T5) = Cu_SH(T5)

const A1::Float64 = 1030.73e-2 #area de contato da face da braid (cm^2)

const r1::Float64 =  50 #cm^2*K/W           
#resistencia de contato da braid com o criostato ( Cu-Cu)

const R1::Float64  = r1/A1

#LM2 é uma mini DAC de aço inox
const r12c::Float64 = 20 
#resistencia de contato entre braid e dac (CU-SS)

const re2::Float64 = 11.25e-3
#mini dac external radius

const ri2::Float64 = 0.35e-3
#dac internal radius

const l2::Float64 = 11.205e-3
#minidac lenght 


const A12::Float64 = (3/2)*π*re2*l2

const R12c::Float64 = r12c*A12



#const TC2::Float64 = Average(4.2,40,SS304_TC)



R2A(T)::Float64 = log(re2/ri2)/((3/2)*π*SS304_TC(T)*l2)

R2B(T)::Float64 = log(re2/ri2)/((3/2)*π*SS304_TC(T)*l2)

#R23 é paramento variavel
#R34 é paramentro variavel
#R45c depende do material mas para essa abordagem vou considerar (braid)CU-SS(sup2)
#const TCcu::Float64= Average(4,40,Cu_TC1)
L12 = 60e-3
TC_m = Cu_TC1(111.4)
A12aux = 0.314 * L12 / TC_m
#aqui, considerando o resultado de resistencia do artigo e a condutividade do cobre nesta temperatura
# eu calculei uma area equivalente para aproximar a braid como um cilindro solido A X L
R12(T)::Float64 = L12/(A12aux*Cu_TC1(T))

R45(T)::Float64 = R12(T)

const R45c::Float64 = R12c

const R5::Float64 = R1


const R23::Array = LinRange(1,1000,size)

const R34::Array = LinRange(1,1000,size)

const M2::Float64 = 51.273e-3 #check minidac massa

c2(T)::Float64 = SS304_SH(T)

const LM3::Array = LinRange(1,200,size)

const LM4::Array = LinRange(1,200,size)



g1 = 1/R1
g12(T) = 1/(R12(T)+R12c+R2A(T))
g23(T) = 1 ./(R2B(T) .+ R23)
g34 = 1 ./R34
g45(T) = 1/(R45(T)+R45c)
g5 = R5





function K_matrix(i,j,T_in)
    K = [g1+g12(T_in[1])   -g12(T_in[2])                     0                        0                      0
         -g12(T_in[1])      g12(T_in[2]).+g23(T_in[2])[i]   -g23(T_in[3])[i]          0                      0 
         0                 -g23(T_in[2])[i]                  g23(T_in[3])[i].+g34[j] -g34[j]                 0
         0                  0                               -g34[j]                   g34[j].+g45(T_in[4])  -g45(T_in[5])
         0                  0                                0                       -g45(T_in[4])           g45(T_in[5])+g5]
    return K
end

function E_matrix(i,j,T_in)
    E = [M1*c1(T_in[1])   0               0      0       0
         0                M2*c2(T_in[2])  0      0       0
         0                0               LM3[i] 0       0
         0                0               0      LM4[j]  0
         0                0               0      0       M5*c5(T_in[5])]  
    return E
end

const L::Array = [1    0   0   g1  0
                  0    1   0   0   0
                  0    0   0   0   0
                  0    0   0   0   0
                  0    0   1   0   g5]



const u::Array = [q1;q2;q3;Ta;Tb]


const D::Array = zeros(5,5)





const T_in::Array = [4.2;4.2;25;40;40]


function stat_T1(Ki,Kj,T_in)  
    A = -inv(E_matrix(size÷2,size÷2,T_in))*K_matrix(Ki,Kj,T_in)
    B =  inv(E_matrix(size÷2,size÷2,T_in))*L
    xstat = -inv(A)*B*u
    T_out1 = I*xstat+D*u
    return T_out1
end






#function steady_state1(size,T_in)
#    T1::Array = zeros(Float64, size, size)
#    T2::Array = zeros(Float64, size, size)
#    T3::Array = zeros(Float64, size, size)
#    T4::Array = zeros(Float64, size, size)
#    T5::Array = zeros(Float64, size, size)
#    T_Place = []
#    
#
#    @threads for i in 1:size
#        @threads for j in 1:size
#            T_Place = stat_T1(i,j,T_in)
#
#
#            T1[i,j] = T_Place[1]
#            T2[i,j] = T_Place[2]
#            T3[i,j] = T_Place[3]
#            T4[i,j] = T_Place[4]
#            T5[i,j] = T_Place[5]
#
#
#            
#        end
#    end
#    return T1, T2, T3, T4, T5
#end

function get_ssT(Ki,Kj,T_in)
    T_out2 = inv(K_matrix(Ki,Kj,T_in))*L*u 
    return T_out2
end
function iter_resistA(size,T_in,iters)
    k=0
    T_prev = T_in
    T_out = zeros(5, size, size)
    δ = [0.1;0.1;0.1;0.1;0.1]
    @threads for i in 1:size
        @threads for j in 1:size
            while true
                T_out[:,i,j]= get_ssT(i,j,T_prev)
                if norm(T_out[:,i,j]-T_prev) <= norm(δ) || k==iters
                    break
                else
                    T_prev = T_out[:,i,j]
                    k=k+1
                end
            end 
        end
    end
    return T_out
end





T_B0 = iter_resistA(size,T_in,0)
contourf(T_B0[1,:,:],title="LM1 - 1x",xlabel="R34",ylabel="R23",colorbar=true)
contourf(T_B0[2,:,:],title="LM2 - 1x",xlabel="R34",ylabel="R23",colorbar=true)
contourf(T_B0[3,:,:],title="LM3 - 1x",xlabel="R34",ylabel="R23",colorbar=true)
contourf(T_B0[4,:,:],title="LM4 - 1x",xlabel="R34",ylabel="R23",colorbar=true)
contourf(T_B0[5,:,:],title="LM5 - 1x",xlabel="R34",ylabel="R23",colorbar=true)



min_val, min_index = findmin(T_B0[2,:,:])
println("The smallest value is $min_val at index $min_index")


#T11 , T21, T31, T41, T51 = steady_state1(size,T_in)
#
#T1B , T2B, T3B, T4B, T5B = steady_state1(size,T_in)÷2

T_B0[:,size,size]
T_B0[1,size,size] = 4.2
#aqui eu acho a matriz K e a E com a temperatura que eu obtive anteriormente
A = -inv(E_matrix(500,500,T_B0[:,size,size]))*K_matrix(size,size,T_B0[:,size,size])
B =  inv(E_matrix(500,500,T_B0[:,size,size]))*L
u0::Array=[q1;q2;q3;Ta;Tb]
X1::Array=-inv(A)*B*u0#nao estou usando este X1 como estado inicial, apenas para verificar o estado que as matrizes que eu encontrei levam
X0::Array=T_B0[:,size,size]


stime = 1*60*60
t = 0:0.1:stime
ut = [fill(q1,1,length(t));fill(q2,1,length(t));fill(q3,1,length(t));fill(Ta,1,length(t));fill(Tb,1,length(t))]
 

sys = ss(A,B,I,D)
ys, ts, xs, us = lsim(sys,ut,t,X0)

transient = plot(ts,ys',label=["LM1"    "LM2"   "LM3"   "LM4"   "LM5"], title="Transient Lumped Mass Temperature")
display(transient)

plot(SS304_TC,0,300)