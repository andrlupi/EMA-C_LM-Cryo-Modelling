"""
Módulo Materials
================

Modelagem analítica das propriedades termofísicas e mecânicas de materiais criogênicos,
integrando dados empíricos do NIST e cálculo exato de fluxo condutivo via Transformada de Kirchhoff.

CONCEITOS FÍSICOS CENTRAIS:
1. CONDUTIVIDADE TÉRMICA k(T):
   - Em metais puros (ex: Cobre OFHC), elétrons dominam a condução de calor. A baixas temperaturas,
     a competição entre espalhamento de elétrons em impurezas e em fônons cria um pico violento
     de condutividade térmica por volta de 15 K - 25 K (~1000 W/m·K), caindo para ~300 W/m·K a 4.2 K.
   - Em ligas desordenadas (ex: Aço Inoxidável SS304), o espalhamento estático em defeitos de rede
     é tão intenso que k(T) é ordens de grandeza menor (~0.29 W/m·K a 4.2 K), subindo monotonicamente
     com a temperatura.

2. CALOR ESPECÍFICO cp(T):
   - Segundo a teoria de Debye, para T << θ_Debye, os modos vibracionais da rede congelam,
     resultando em cp ∝ T³ (além do termo linear eletrônico γT para metais).
   - Isso significa que a 4 K os materiais têm capacidade térmica quase nula; qualquer
     aporte minúsculo de calor (ex: feixe síncrotron de 10 mW) causa variações rápidas de temperatura.

3. TRANSFORMADA DE KIRCHHOFF E INTEGRAL DE CONDUTIVIDADE:
   - A lei de Fourier unidimensional é:
         q = -k(T) * (dT/dx)
   - Discretizar com a condutividade na temperatura média k((T1+T2)/2) gera erros graves porque
     k(T) não é linear no intervalo.
   - A solução exata integrada entre as temperaturas T1 e T2 ao longo de uma barra é:
         Q = (A / L) * ∫_{T1}^{T2} k(T) dT
   - Definindo o potencial térmico Θ(T) = ∫_0^T k(T') dT', o fluxo de calor torna-se estritamente linear:
         Q = (A / L) * [Θ(T2) - Θ(T1)]
"""
module Materials

using LinearAlgebra
include("NIST_Data.jl")
using .NIST_Data

export AbstractMaterial,
       CopperOFHC,
       StainlessSteel304,
       BerylliumCopper,
       Beryllium,
       TitaniumTi6Al4V,
       thermal_conductivity,
       specific_heat,
       density,
       youngs_modulus,
       thermal_expansion,
       thermal_conductivity_integral

"""
    AbstractMaterial

Super-tipo abstrato para qualquer material criogênico no sistema.
Subtipos implementam:
- `thermal_conductivity(mat, T)`: k(T) em W/(m·K)
- `specific_heat(mat, T)`: cp(T) em J/(kg·K)
- `density(mat)`: ρ em kg/m³
- `youngs_modulus(mat, T)`: E(T) em Pa (opcional para componentes estruturais)
- `thermal_expansion(mat, T)`: ΔL/L_293 (adimensional)
"""
abstract type AbstractMaterial end

# =============================================================================
# QUADRATURA DE GAUSS-LEGENDRE DE 7 PONTOS
# =============================================================================
# POR QUE USAR GAUSS-LEGENDRE EM VEZ DE TRAPÉZIOS OU SIMPSON?
# 1. Eficiência: Uma quadratura de Gauss de n pontos integra de forma EXATA qualquer
#    polinômio de grau até 2n - 1. Com n = 7, ela é exata para polinômios de grau até 13!
# 2. Desempenho: Exige exatamente 7 avaliações de k(T), enquanto métodos trapezoidais
#    precisariam de centenas de passos para obter erro < 10⁻⁵.
# 3. Diferenciabilidade e Suavidade: Por ser uma soma ponderada analítica fixa,
#    ela preserva a suavidade C^∞ das funções NIST, ideal para solvers de Newton-Raphson.
#
# Nós (raízes dos polinômios de Legendre P_7(x)) no intervalo [-1, 1]:
const GL_X = [
    -0.9491079123427585,
    -0.7415311855993944,
    -0.4058451513773972,
     0.0,
     0.4058451513773972,
     0.7415311855993944,
     0.9491079123427585
]

# Pesos correspondentes de quadratura:
const GL_W = [
    0.1294849661688697,
    0.2797053914892766,
    0.3818300505051189,
    0.4179591836734694,
    0.3818300505051189,
    0.2797053914892766,
    0.1294849661688697
]

"""
    thermal_conductivity_integral(mat, T1, T2)

Calcula a Integral de Condutividade Térmica entre T1 e T2 [W/m]:
    Θ(T1, T2) = ∫_{T1}^{T2} k(T) dT

Faz a mudança de variável linear de [T1, T2] para [-1, 1]:
    T(x) = (T1 + T2)/2 + ((T2 - T1)/2) * x
    dT = ((T2 - T1)/2) * dx
"""
function thermal_conductivity_integral(mat::AbstractMaterial, T1::Real, T2::Real)
    if T1 == T2
        return zero(promote_type(typeof(T1), typeof(T2)))
    end
    
    half_span = (T2 - T1) / 2
    midpoint  = (T1 + T2) / 2
    integral  = zero(promote_type(typeof(T1), typeof(T2)))
    
    # Avaliação nos 7 pontos de Gauss
    for i in 1:7
        T_eval = midpoint + half_span * GL_X[i]
        integral += GL_W[i] * thermal_conductivity(mat, T_eval)
    end
    
    return half_span * integral
end

# -----------------------------------------------------------------------------
# Cache Global para os dados NIST (carregados uma única vez na memória)
# -----------------------------------------------------------------------------
const NIST_REF = Ref{Union{Nothing, NISTMaterialData}}(nothing)

function get_nist()
    if NIST_REF[] === nothing
        NIST_REF[] = load_nist_data()
    end
    return NIST_REF[]
end

# =============================================================================
# 1. COBRE ELETROLÍTICO OFHC (Oxygen-Free High Conductivity)
# =============================================================================
"""
    CopperOFHC(rrr=100; density=8960.0)

Cobre de alta pureza e livre de oxigênio (ASTM C10100/C10200).
- `rrr`: Residual Resistance Ratio (razão entre resistividade elétrica a 273 K e 4.2 K).
  Tipicamente RRR = 50 a 100 para cordoalhas e condutores criogênicos comerciais.
- Utilizado no sistema: Cordoalhas flexíveis (thermal braids) de extração térmica da DAC.
"""
struct CopperOFHC <: AbstractMaterial
    rrr::Int
    density_val::Float64
    function CopperOFHC(rrr::Int = 100; density = 8960.0)
        new(rrr, density)
    end
end

density(m::CopperOFHC) = m.density_val

"""
Condutividade Térmica do Cobre OFHC pelo ajuste fracionário racional do NIST:
    log10(k) = (a + c*T^0.5 + e*T + g*T^1.5 + i*T^2) / (1 + b*T^0.5 + d*T + f*T^1.5 + h*T^2)
"""
function thermal_conductivity(m::CopperOFHC, T::Real)
    nist = get_nist()
    # TC1 da tabela NIST corresponde a cobre com RRR ~ 100
    num_p = [1.0, T^0.5, T, T^1.5, T^2]
    den_p = [T^0.5, T, T^1.5, T^2]
    
    # Produto escalar com os coeficientes do numerador e denominador
    num = dot(nist.cu_c[1:2:9, 1], num_p)
    den = 1.0 + dot(nist.cu_c[2:2:9, 1], den_p)
    return 10.0^(num / den)
end

"""
Calor Específico do Cobre OFHC por polinômio logarítmico NIST:
    log10(cp) = ∑_{n=0}^8 a_n * (log10(T))^n
"""
function specific_heat(m::CopperOFHC, T::Real)
    nist = get_nist()
    logT = log10(max(T, 1.0)) # salvaguarda contra T <= 0
    p = [logT^i for i in 0:8]
    val = dot(nist.cu_c[1:9, 6], p)
    return 10.0^val
end

youngs_modulus(m::CopperOFHC, T::Real) = 128e9 # ~128 GPa em temperaturas criogênicas

# =============================================================================
# 2. AÇO INOXIDÁVEL AUSTENÍTICO AISI 304 (SS304)
# =============================================================================
"""
    StainlessSteel304(; density=7900.0)

Aço inoxidável austenítico com 18% Cr e 8% Ni.
- Características criogênicas: Mantém tenacidade sem transição dúctil-frágil acentuada;
  possui baixíssima condutividade térmica a 4 K (0.29 W/m·K) e alto módulo de rigidez (~210 GPa).
- Utilizado no sistema: Suportes estruturais tubulares/hastes de sustentação da DAC.
"""
struct StainlessSteel304 <: AbstractMaterial
    density_val::Float64
    function StainlessSteel304(; density = 7900.0)
        new(density)
    end
end

density(m::StainlessSteel304) = m.density_val

"""
Condutividade Térmica do Inox 304 por expansão logarítmica NIST (coluna TC):
    log10(k) = ∑_{n=0}^8 a_n * (log10(T))^n
"""
function thermal_conductivity(m::StainlessSteel304, T::Real)
    nist = get_nist()
    logT = log10(max(T, 1.0))
    p = [logT^i for i in 0:8]
    val = dot(nist.ss304_c[1:9, 1], p)
    return 10.0^val
end

"""
Calor Específico do Inox 304 (coluna SH):
    log10(cp) = ∑_{n=0}^7 a_n * (log10(T))^n
"""
function specific_heat(m::StainlessSteel304, T::Real)
    nist = get_nist()
    logT = log10(max(T, 1.0))
    p = [logT^i for i in 0:7]
    val = dot(nist.ss304_c[1:8, 2], p)
    return 10.0^val
end

"""
Módulo de Young E(T) em Pascals [N/m²] do Inox 304 (coluna YM1):
    E(T) = a + b*T + c*T^2 + d*T^3 + e*T^4  [em GPa, convertido para Pa (* 1e9)]
"""
function youngs_modulus(m::StainlessSteel304, T::Real)
    nist = get_nist()
    p = [1.0, T, T^2, T^3, T^4]
    E_GPa = dot(nist.ss304_c[1:5, 3], p)
    return E_GPa * 1e9
end

"""
Contração Térmica Integrada ΔL/L_293 do Inox 304 (coluna LE):
Retorna a variação fracionária de comprimento em relação a 293 K.
"""
function thermal_expansion(m::StainlessSteel304, T::Real)
    nist = get_nist()
    p = [1.0, T, T^2, T^3, T^4]
    return dot(nist.ss304_c[1:5, 5], p) * 1e-5
end

# =============================================================================
# 3. COBRE-BERÍLIO (CuBe - Alloy 25 / C17200)
# =============================================================================
"""
    BerylliumCopper(; density=8250.0)

Liga de cobre com ~1.8% a 2% de Berílio endurecida por precipitação.
- Características: Une altíssima resistência mecânica (tensão de escoamento > 1 GPa),
  ausência de magnetismo (crítico para linhas de luz síncrotron) e condutividade térmica
  intermediária entre o cobre puro e o inox.
- Utilizado no sistema: Corpo da célula de bigorna de diamante (Mini-DAC).
"""
struct BerylliumCopper <: AbstractMaterial
    density_val::Float64
    function BerylliumCopper(; density = 8250.0)
        new(density)
    end
end

density(m::BerylliumCopper) = m.density_val

"""
Condutividade Térmica do CuBe por expansão logarítmica NIST:
    log10(k) = ∑_{n=0}^7 a_n * (log10(T))^n
"""
function thermal_conductivity(m::BerylliumCopper, T::Real)
    nist = get_nist()
    logT = log10(max(T, 1.0))
    p = [logT^i for i in 0:7]
    val = dot(nist.cube_c[1:8, 1], p)
    return 10.0^val
end

"""
Calor Específico do CuBe pela Lei de Kopp-Neumann:
Em ligas metálicas sólidas, a capacidade calorífica molar é a média ponderada
das frações molares dos componentes:
    cp(CuBe) ≈ 0.95 * cp(Cu) + 0.05 * cp(Be)
"""
function specific_heat(m::BerylliumCopper, T::Real)
    cu = CopperOFHC()
    be = Beryllium()
    return 0.95 * specific_heat(cu, T) + 0.05 * specific_heat(be, T)
end

youngs_modulus(m::BerylliumCopper, T::Real) = 131e9 # ~131 GPa em frio criogênico

# =============================================================================
# 4. BERÍLIO PURO (Be)
# =============================================================================
"""
    Beryllium(; density=1850.0)

Metal de baixo número atômico (Z=4) e baixíssima densidade.
- Características: Excelente transparência a raios X; utilizado em janelas ópticas
  e nos cálculos de capacidade calorífica da liga CuBe via Kopp-Neumann.
"""
struct Beryllium <: AbstractMaterial
    density_val::Float64
    function Beryllium(; density = 1850.0)
        new(density)
    end
end

density(m::Beryllium) = m.density_val

# =============================================================================
# 5. TITÂNIO GRAU 5 (Ti-6Al-4V)
# =============================================================================
"""
    TitaniumTi6Al4V(; density=4430.0)

Liga aeroespacial de titânio (Ti-6%Al-4%V).
- Características criogênicas: É uma das alternativas mais famosas em criogenia de precisão.
  Possui condutividade térmica ainda menor que o Inox 304 (~0.12 W/m·K a 4.2 K)
  com excelente resistência mecânica e densidade quase metade da do aço (4430 kg/m³).
- Utilizado no sistema: Alternativa estrutural avançada para isolamento dos suportes.
"""
struct TitaniumTi6Al4V <: AbstractMaterial
    density_val::Float64
    function TitaniumTi6Al4V(; density = 4430.0)
        new(density)
    end
end

density(m::TitaniumTi6Al4V) = m.density_val

"""
Condutividade Térmica do Ti-6Al-4V por ajuste empírico NIST:
    log10(k) = -1.944 + 1.258*log10(T) - 0.177*(log10(T))^2 + 0.015*(log10(T))^3
"""
function thermal_conductivity(m::TitaniumTi6Al4V, T::Real)
    logT = log10(max(T, 1.0))
    val = -1.944 + 1.258*logT - 0.177*(logT^2) + 0.015*(logT^3)
    return 10.0^val
end

"""
Calor Específico do Ti-6Al-4V:
    log10(cp) = -2.18 + 2.52*log10(T) - 0.49*(log10(T))^2
"""
function specific_heat(m::TitaniumTi6Al4V, T::Real)
    logT = log10(max(T, 1.0))
    val = -2.18 + 2.52*logT - 0.49*(logT^2)
    return 10.0^val
end

youngs_modulus(m::TitaniumTi6Al4V, T::Real) = 118e9 # ~118 GPa em criogenia

end # module Materials
