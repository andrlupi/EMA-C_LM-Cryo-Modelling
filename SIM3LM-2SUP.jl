#André Luiz Pianca - Estagiario - Grupo MArÉ - 29/09/2023
#
#

using Trapz
using DelimitedFiles
using LinearAlgebra
include(raw"Material Properties\Properties.jl")


R1 = 1 #3 casos, para cada alternativa de contato com o cryo
g1 = 1/R1

m_dac = 0.439
c_dac = #properties.jl
DAC = m_dac*c_dac

m_s1= 0.332
c_s1= Average(40,300,SS304_SH)#properties.jl
S1 = m_s1*c_s1
S2 = S1

R12 = 
R23 = 
R3  = 








g12 = 1/R12
g23 = 1/R23
g3 = 1/R3

E = [DAC  0   0
      0   S1  0
      0   0   S2]
K = [g1+g12  -g12          0
    -g12     +g12+g23    -g23
     0       -g23        +g23+g3]