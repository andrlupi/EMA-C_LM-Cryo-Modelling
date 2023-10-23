using DelimitedFiles
using Trapz

cu_data = readdlm(raw"Material Properties\CU-TC1_TC2+TC3_TC4_TC5_SH_EC.csv", ',')
be_data = readdlm(raw"Material Properties\BE-SH_LE1_LE2_LE3.csv", ',')
cube_data =  readdlm(raw"Material Properties\CuBe-TC_LE.csv", ',')
SS304_data = readdlm(raw"Material Properties\SS304-TC_SH_YM_YM_LE.csv", ',')
Cu_c = cu_data[2:10,2:8]
Be_c = be_data[2:10,2:5]
CuBe_c = cube_data[2:10,2:3]
SS304_c = SS304_data[2:10,2:6]




n = 0:1:8
hrad(T) =σ*ϵ*(T^2+4.2)*(T+4.2)
Cu_TC1(T) = 10^((dot(Cu_c[1:2:9,1],T.^(0:0.5:2)))/(1+dot(Cu_c[2:2:9,1],T.^(0.5:0.5:2))))
Cu_SH(T) = 10^(dot(Cu_c[1:9,6],log10(T).^n))
Be_SH(T) = 10^(dot(Be_c[1:9,1],log10(T).^n))
CuBe_TC(T) = 10^(dot(CuBe_c[1:9,1],log10(T).^n))
SS304_TC(T) = 10^(dot(SS304_c[1:9,1],log10(T).^n))
SS304_SH(T) = 10^(dot(SS304_c[1:9,2],log10(T).^n))




function Average(Tmin,Tmax,f::Function)
    tv = range(Tmin,Tmax,length=1000)
    Average= trapz(tv,[f(T) for T=tv])/(Tmax-Tmin)
end