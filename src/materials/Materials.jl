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
using ForwardDiff
include("NIST_Data.jl")
using .NIST_Data

# Desempacotamento escalar recursivo para números Dual do ForwardDiff de qualquer profundidade
_scalar_val(x::Real) = Float64(x)
_scalar_val(x::ForwardDiff.Dual) = _scalar_val(ForwardDiff.value(x))

# Constantes físicas de decaimento sub-Kelvin linear (Terceira Lei da Termodinâmica)
const CU_TC_1 = 136.8184232866769 # k(1.0 K) exato do ajuste NIST
const SS304_CP_1 = 2.065985684821771 / 4.0 # cp(1.0 K) = 0.5164964212054427 J/(kg·K)
const BE_CP_1 = 0.0003318870374978916 / 5.0 # cp(1.0 K) = 6.637740749957832e-5 J/(kg·K)

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
# QUADRATURA DE GAUSS-LEGENDRE DE 7 PONTOS (ESTÁTICA, ZERO ALOCAÇÕES)
# =============================================================================
const GL_X = (
    -0.9491079123427585,
    -0.7415311855993944,
    -0.4058451513773972,
     0.0,
     0.4058451513773972,
     0.7415311855993944,
     0.9491079123427585,
)

const GL_W = (
    0.1294849661688697,
    0.2797053914892766,
    0.3818300505051189,
    0.4179591836734694,
    0.3818300505051189,
    0.2797053914892766,
    0.1294849661688697,
)

"""
    thermal_conductivity_integral(mat, T1, T2)

Calcula a Integral de Condutividade Térmica entre T1 e T2 [W/m]:
    Θ(T1, T2) = ∫_{T1}^{T2} k(T) dT

Avaliação estática por quadratura de Gauss-Legendre de 7 pontos com zero alocações na heap.
"""
function thermal_conductivity_integral(mat::AbstractMaterial, T1::Real, T2::Real)
    if T1 == T2
        return zero(promote_type(typeof(T1), typeof(T2), Float64))
    end
    
    half_span = (T2 - T1) * 0.5
    midpoint  = (T1 + T2) * 0.5
    
    # Avaliação desenrolada nos 7 nós de Gauss para evitar alocação de iteradores ou vetores
    s = GL_W[1] * thermal_conductivity(mat, midpoint + half_span * GL_X[1]) +
        GL_W[2] * thermal_conductivity(mat, midpoint + half_span * GL_X[2]) +
        GL_W[3] * thermal_conductivity(mat, midpoint + half_span * GL_X[3]) +
        GL_W[4] * thermal_conductivity(mat, midpoint + half_span * GL_X[4]) +
        GL_W[5] * thermal_conductivity(mat, midpoint + half_span * GL_X[5]) +
        GL_W[6] * thermal_conductivity(mat, midpoint + half_span * GL_X[6]) +
        GL_W[7] * thermal_conductivity(mat, midpoint + half_span * GL_X[7])
    
    return half_span * s
end

# -----------------------------------------------------------------------------
# Cache Global para os dados NIST (carregados uma única vez na memória)
# -----------------------------------------------------------------------------
const NIST_REF = Ref{Union{Nothing, NISTMaterialData}}(nothing)

function get_nist()::NISTMaterialData
    if NIST_REF[] === nothing
        NIST_REF[] = load_nist_data()
    end
    return (NIST_REF[])::NISTMaterialData
end

function __init__()
    get_nist()
end

# =============================================================================
# 1. COBRE ELETROLÍTICO OFHC (Oxygen-Free High Conductivity)
# =============================================================================
"""
    CopperOFHC(rrr=100; density=8960.0)

Cobre de alta pureza e livre de oxigênio (ASTM C10100/C10200).
- `rrr`: Residual Resistance Ratio (razão entre resistividade elétrica a 273 K e 4.2 K).
- Utilizado no sistema: Cordoalhas flexíveis (thermal braids) de extração térmica da DAC.
"""
struct CopperOFHC <: AbstractMaterial
    rrr::Int
    density_val::Float64
    function CopperOFHC(rrr::Int = 100; density::Real = 8960.0)
        new(rrr, Float64(density))
    end
end

density(m::CopperOFHC) = m.density_val

"""
Condutividade Térmica do Cobre OFHC pelo ajuste fracionário racional do NIST:
    log10(k) = (a + c*T^0.5 + e*T + g*T^1.5 + i*T^2) / (1 + b*T^0.5 + d*T + f*T^1.5 + h*T^2)
"""
function thermal_conductivity(m::CopperOFHC, T::Real)
    if T > 350.0
        T_num = _scalar_val(T)
        @warn "Temperature $T_num K exceeds NIST validity range (T <= 350 K)" maxlog=1
        T_eval = typeof(T)(350.0)
    elseif T < 1.0
        # Terceira Lei da Termodinâmica: k -> 0 quando T -> 0 com derivada não-nula para diferenciação automática
        return CU_TC_1 * ifelse(T >= zero(T), T, zero(T))
    else
        T_eval = T
    end
    
    nist = get_nist()
    u = sqrt(T_eval)
    num = evalpoly(u, nist.cu_tc_num)
    den = evalpoly(u, nist.cu_tc_den)
    return 10.0^(num / den)
end

"""
Calor Específico do Cobre OFHC por polinômio logarítmico NIST:
    log10(cp) = ∑_{n=0}^8 a_n * (log10(T))^n
"""
function specific_heat(m::CopperOFHC, T::Real)
    nist = get_nist()
    if T > 350.0
        T_num = _scalar_val(T)
        @warn "Temperature $T_num K exceeds NIST validity range (T <= 350 K)" maxlog=1
        T_eval = typeof(T)(350.0)
    elseif T < 1.0
        # Em regime sub-Kelvin, calor específico eletrônico linear cp = γ*T domina com derivada finita
        cp_1 = 10.0^nist.cu_sh[1]
        return cp_1 * ifelse(T >= zero(T), T, zero(T))
    else
        T_eval = T
    end
    
    logT = log10(T_eval)
    val = evalpoly(logT, nist.cu_sh)
    return 10.0^val
end

youngs_modulus(m::CopperOFHC, T::Real) = 128e9 # ~128 GPa em temperaturas criogênicas

function thermal_expansion(m::CopperOFHC, T::Real)
    T_eval = clamp(T, 0.0, 350.0)
    nist = get_nist()
    return evalpoly(T_eval, nist.cu_le) * 1e-5
end

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
    function StainlessSteel304(; density::Real = 7900.0)
        new(Float64(density))
    end
end

density(m::StainlessSteel304) = m.density_val

"""
Condutividade Térmica do Inox 304 por expansão logarítmica NIST (coluna TC):
    log10(k) = ∑_{n=0}^8 a_n * (log10(T))^n
"""
function thermal_conductivity(m::StainlessSteel304, T::Real)
    nist = get_nist()
    if T > 350.0
        T_num = _scalar_val(T)
        @warn "Temperature $T_num K exceeds NIST validity range (T <= 350 K)" maxlog=1
        T_eval = typeof(T)(350.0)
    elseif T < 1.0
        k_1 = 10.0^nist.ss304_tc[1]
        return k_1 * ifelse(T >= zero(T), T, zero(T))
    else
        T_eval = T
    end
    
    logT = log10(T_eval)
    val = evalpoly(logT, nist.ss304_tc)
    return 10.0^val
end

"""
Calor Específico do Inox 304 (coluna SH):
    log10(cp) = ∑_{n=0}^7 a_n * (log10(T))^n
"""
function specific_heat(m::StainlessSteel304, T::Real)
    nist = get_nist()
    if T > 350.0
        T_num = _scalar_val(T)
        @warn "Temperature $T_num K exceeds NIST validity range (T <= 350 K)" maxlog=1
        T_eval = typeof(T)(350.0)
    elseif T < 4.0
        # NIST SH para SS304 é válido para T >= 4.0 K. Abaixo de 4.0 K, o polinômio diverge artificialmente para 1e22.
        # Adota-se decaimento físico linear cp(T) = SS304_CP_1 * T (~0.5165 T J/kg·K).
        # Para T <= 1.0 K, garante cp(1.0) = SS304_CP_1 > 0 e Terceira Lei com derivada positiva.
        return SS304_CP_1 * ifelse(T >= zero(T), T, zero(T))
    else
        T_eval = T
    end
    
    logT = log10(T_eval)
    val = evalpoly(logT, nist.ss304_sh)
    return 10.0^val
end

"""
Módulo de Young E(T) em Pascals [N/m²] do Inox 304 (coluna YM2: válida de 4 K a 300 K):
    E(T) = a + b*T + c*T^2 + d*T^3 + e*T^4  [em GPa, convertido para Pa (* 1e9)]
"""
function youngs_modulus(m::StainlessSteel304, T::Real)
    T_eval = clamp(T, 0.0, 350.0)
    nist = get_nist()
    E_GPa = evalpoly(T_eval, nist.ss304_ym)
    return E_GPa * 1e9
end

"""
Contração Térmica Integrada ΔL/L_293 do Inox 304 (coluna LE):
Retorna a variação fracionária de comprimento em relação a 293 K.
"""
function thermal_expansion(m::StainlessSteel304, T::Real)
    T_eval = clamp(T, 0.0, 350.0)
    nist = get_nist()
    return evalpoly(T_eval, nist.ss304_le) * 1e-5
end

# =============================================================================
# 3. BERÍLIO PURO (Be)
# =============================================================================
"""
    Beryllium(; density=1850.0)

Metal de baixo número atômico (Z=4) e baixíssima densidade.
- Características: Excelente transparência a raios X; utilizado em janelas ópticas
  e nos cálculos de capacidade calorífica da liga CuBe via Kopp-Neumann.
"""
struct Beryllium <: AbstractMaterial
    density_val::Float64
    function Beryllium(; density::Real = 1850.0)
        new(Float64(density))
    end
end

density(m::Beryllium) = m.density_val

const BE_TC_COEFFS = (-0.45, 2.05, -0.62, 0.051)

"""
Condutividade Térmica do Berílio policristalino:
    log10(k) = ∑ a_n * (log10(T))^n
"""
function thermal_conductivity(m::Beryllium, T::Real)
    if T > 350.0
        T_num = _scalar_val(T)
        @warn "Temperature $T_num K exceeds NIST validity range (T <= 350 K)" maxlog=1
        T_eval = typeof(T)(350.0)
    elseif T < 1.0
        k_1 = 10.0^BE_TC_COEFFS[1]
        return k_1 * ifelse(T >= zero(T), T, zero(T))
    else
        T_eval = T
    end
    logT = log10(T_eval)
    val = evalpoly(logT, BE_TC_COEFFS)
    return 10.0^val
end

"""
Calor Específico do Berílio puro por expansão logarítmica NIST (coluna SH):
    log10(cp) = ∑_{n=0}^8 a_n * (log10(T))^n
"""
function specific_heat(m::Beryllium, T::Real)
    nist = get_nist()
    if T > 350.0
        T_num = _scalar_val(T)
        @warn "Temperature $T_num K exceeds NIST validity range (T <= 350 K)" maxlog=1
        T_eval = typeof(T)(350.0)
    elseif T < 5.0
        # NIST SH para Be puro é válido para T >= 5.0 K. Abaixo de 5.0 K, o polinômio submerge em underflow (10^-526 a 1 K).
        # Adota-se decaimento físico linear cp(T) = BE_CP_1 * T (~6.638e-5 T J/kg·K).
        # Para T <= 1.0 K, garante cp(1.0) = BE_CP_1 > 0 e Terceira Lei com derivada positiva.
        return BE_CP_1 * ifelse(T >= zero(T), T, zero(T))
    else
        T_eval = T
    end
    
    logT = log10(T_eval)
    val = evalpoly(logT, nist.be_sh)
    return 10.0^val
end

youngs_modulus(m::Beryllium, T::Real) = 287e9 # ~287 GPa em criogenia

function thermal_expansion(m::Beryllium, T::Real)
    T_eval = clamp(T, 0.0, 350.0)
    nist = get_nist()
    return evalpoly(T_eval, nist.be_le) * 1e-5
end

# =============================================================================
# 4. COBRE-BERÍLIO (CuBe - Alloy 25 / C17200)
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
    function BerylliumCopper(; density::Real = 8250.0)
        new(Float64(density))
    end
end

density(m::BerylliumCopper) = m.density_val

const DEFAULT_CU = CopperOFHC()
const DEFAULT_BE = Beryllium()

"""
Condutividade Térmica do CuBe por expansão logarítmica NIST:
    log10(k) = ∑_{n=0}^7 a_n * (log10(T))^n
"""
function thermal_conductivity(m::BerylliumCopper, T::Real)
    nist = get_nist()
    if T > 350.0
        T_num = _scalar_val(T)
        @warn "Temperature $T_num K exceeds NIST validity range (T <= 350 K)" maxlog=1
        T_eval = typeof(T)(350.0)
    elseif T < 1.0
        k_1 = 10.0^nist.cube_tc[1]
        return k_1 * ifelse(T >= zero(T), T, zero(T))
    else
        T_eval = T
    end
    
    logT = log10(T_eval)
    val = evalpoly(logT, nist.cube_tc)
    return 10.0^val
end

"""
Calor Específico do CuBe pela Lei de Kopp-Neumann:
Em ligas metálicas sólidas, a capacidade calorífica molar é a média ponderada
das frações molares dos componentes:
    cp(CuBe) ≈ 0.95 * cp(Cu) + 0.05 * cp(Be)
"""
function specific_heat(m::BerylliumCopper, T::Real)
    return 0.95 * specific_heat(DEFAULT_CU, T) + 0.05 * specific_heat(DEFAULT_BE, T)
end

youngs_modulus(m::BerylliumCopper, T::Real) = 131e9 # ~131 GPa em frio criogênico

function thermal_expansion(m::BerylliumCopper, T::Real)
    T_eval = clamp(T, 0.0, 350.0)
    nist = get_nist()
    return evalpoly(T_eval, nist.cube_le) * 1e-5
end

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
    function TitaniumTi6Al4V(; density::Real = 4430.0)
        new(Float64(density))
    end
end

density(m::TitaniumTi6Al4V) = m.density_val

const TI_TC_COEFFS = (-1.944, 1.258, -0.177, 0.015)
const TI_SH_COEFFS = (-2.18, 2.52, -0.49)
const TI_LE_COEFFS = (-173.9439, -0.20462, 0.0087853, -3.37692e-5, 4.46548e-8)

"""
Condutividade Térmica do Ti-6Al-4V por ajuste empírico NIST:
    log10(k) = -1.944 + 1.258*log10(T) - 0.177*(log10(T))^2 + 0.015*(log10(T))^3
"""
function thermal_conductivity(m::TitaniumTi6Al4V, T::Real)
    if T > 350.0
        T_num = _scalar_val(T)
        @warn "Temperature $T_num K exceeds NIST validity range (T <= 350 K)" maxlog=1
        T_eval = typeof(T)(350.0)
    elseif T < 1.0
        k_1 = 10.0^TI_TC_COEFFS[1]
        return k_1 * ifelse(T >= zero(T), T, zero(T))
    else
        T_eval = T
    end
    
    logT = log10(T_eval)
    val = evalpoly(logT, TI_TC_COEFFS)
    return 10.0^val
end

"""
Calor Específico do Ti-6Al-4V:
    log10(cp) = -2.18 + 2.52*log10(T) - 0.49*(log10(T))^2
"""
function specific_heat(m::TitaniumTi6Al4V, T::Real)
    if T > 350.0
        T_num = _scalar_val(T)
        @warn "Temperature $T_num K exceeds NIST validity range (T <= 350 K)" maxlog=1
        T_eval = typeof(T)(350.0)
    elseif T < 1.0
        cp_1 = 10.0^TI_SH_COEFFS[1]
        return cp_1 * ifelse(T >= zero(T), T, zero(T))
    else
        T_eval = T
    end
    
    logT = log10(T_eval)
    val = evalpoly(logT, TI_SH_COEFFS)
    return 10.0^val
end

youngs_modulus(m::TitaniumTi6Al4V, T::Real) = 118e9 # ~118 GPa em criogenia

function thermal_expansion(m::TitaniumTi6Al4V, T::Real)
    T_eval = clamp(T, 0.0, 350.0)
    return evalpoly(T_eval, TI_LE_COEFFS) * 1e-5
end

end # module Materials
