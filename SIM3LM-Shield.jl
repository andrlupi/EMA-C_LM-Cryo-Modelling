#André Luiz Pianca - Estagiario - Grupo MArÉ - 29/09/2023
#
#considerando as dimenções maximas do sistema em 100mm de altura e ⌀150mm

import PhysicalConstants.CODATA2018:σ
using Trapz
using DelimitedFiles
using LinearAlgebra


mm = 1e-3
cu_data = readdlm(raw"Material Properties\CU-TC1_TC2+TC3_TC4_TC5_SH_EC.csv", ',')
be_data = readdlm(raw"Material Properties\BE-SH_LE1_LE2_LE3.csv", ',')
cube_data = 
Cu_c = cu_data[2:10,2:8]
Be_c = be_data[2:10,2:4]

#funções // vou salvar em arquivos separados para invocar elas de forma organizada
hrad(T) =σ*ϵ*(T^2+4.2)*(T+4.2)
Cu_TC1(T) = 10^((dot(Cu_c[1:2:9,1],T.^(0:0.5:2)))/(1+dot(Cu_c[2:2:9,1],T.^(0.5:0.5:2))))
Cu_SH(T) = 10^(dot(Cu_c[1:9,6],log10(T).^(0:1:8)))
Be_SH(T) = 10^(dot(Be_c[1:9,2],log10(T).^(0:1:8)))

#lumped quantities vector
M = zeros(3,1)
C = zeros(3,1)

# MODELANDO O ESCUDO
H = 100mm                    # da altura maxima
D = 150mm                   # do diametro maximo
R = D/2                     #raio
Lado = 2*pi*R*H             #area do lado do escudo
Tampa = pi*R^2              #area superior do escudo
Espessura = 3*mm            #estimando
A =Lado+2*Tampa              #area
V = A*Espessura             #volume
ρ = 8.9e+3                   #densidade
M[3] = ρ*V                     #massa
ϵ = 0.04                        #emissividade do cobre polido





tv = range(40,300,length=1000)
hrad_m = trapz(tv,[hrad(T) for T=tv])/(300-40)   
Cu_TC_m = trapz(tv,[Cu_TC1(T) for T=tv])/(300-40)
Cu_S̄H̄ = trapz(tv,[Cu_SH(T) for T=tv])/(300-40)


tv2 = range(4.2,300,length=1000)
Be_SH_m =  trapz(tv2,[Be_SH(T) for T=tv2])/(300-4.2)
Cu_TC_m2 = trapz(tv2,[Cu_TC1(T) for T=tv2])/(300-4.2)



LM = Cu_SH_m*M[3]                 #cap termica media

Biot_shield = hrad_m/Cu_TC_m 

g13 = hrad_m * A
g35 = Espessura/(Cu_TC_m*Tampa)

#DAC

M[1] = 0.5 #massa da DAC


CuBe_SH_m = 0.95*Cu_TC_m2 +0.05*Be_SH_m   #Kopp's law p/ 95% cobre 5% berilio 

DAC = M[1]*CuBe_SH_m   # LM da dac

#numero de biot da DAC

#resistencia interna da DAC

