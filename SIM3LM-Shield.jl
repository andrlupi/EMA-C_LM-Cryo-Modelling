#André Luiz Pianca - Estagiario - Grupo MArÉ - 29/09/2023
#
#considerando as dimenções maximas do sistema em 100mm de altura e ⌀150mm

import PhysicalConstants.CODATA2018:σ
using Trapz
using DelimitedFiles
using LinearAlgebra
include(raw"Material Properties\Properties.jl")


mm = 1e-3

#funções // vou salvar em arquivos separados para invocar elas de forma organizada posteriormente




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



#integrando os valores medios das quantidades
#medias entre 40K e 300K(suporte, e escudos)
tv = range(40,300,length=1000)
hrad_m = trapz(tv,[hrad(T) for T=tv])/(300-40)   
Cu_TC_m = trapz(tv,[Cu_TC1(T) for T=tv])/(300-40)
Cu_SH_m = trapz(tv,[Cu_SH(T) for T=tv])/(300-40)
SS304_TC_m = trapz(tv,[SS304_TC(T) for T=tv])/(300-40)




#medias entre 4.2K e 300K (amostra)
tv2 = range(4.2,300,length=1000)
Be_SH_m =  trapz(tv2,[Be_SH(T) for T=tv2])/(300-4.2)
Cu_SH_m2 = trapz(tv2,[Cu_SH(T) for T=tv2])/(300-4.2)
CuBe_TC_m = trapz(tv2,[CuBe_TC(T) for T=tv2])/(300-4.2)


SHIELD= Cu_SH_m*M[3]                 #cap termica media do escudo


Biot_shield = hrad_m/Cu_TC_m 

g13 = hrad_m * A
g35 = Espessura/(Cu_TC_m*Tampa)

#DAC

M[1] = 0.5 #massa da DAC


CuBe_SH_m = 0.95*Cu_SH_m2 +0.05*Be_SH_m   #Kopp's law p/ 95% cobre 5% berilio 

DAC = M[1]*CuBe_SH_m   # LM da dac

#numero de biot da DAC

#resistencia interna da DAC
#simplificando que a DAC é um cilindro furado
r1 = 0.35mm

r2 = 27mm

DAC_L = 37mm

Rint_DAC = log(r2/r1)/(2π*CuBe_TC_m*DAC_L )

#caso radial - radiação

Bi_DAC_rad = (Rint_DAC)*(A*hrad_m)

#caso condução, simplificando que o suporte é um anel
r3 = 36mm

SUP_L = 25mm

R_Sup = log(r3/r2)/(2π*SS304_TC_m*SUP_L)

Bi_DAC_cond = (Rint_DAC)/(R_Sup)

Bi_max = 0.1
