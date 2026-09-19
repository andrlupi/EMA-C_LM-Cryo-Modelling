"""
Módulo ColdInterfaces
=====================

Modelagem física detalhada das duas principais alternativas de interface térmica entre
o dedo frio do criocooler e a amostra/DAC:

1. CORDOALHA FLEXÍVEL COM JUNTA APARAFUSADA E FOLHA DE ÍNDIO:
   - Mecanismo: Condução por contato mecânico em vácuo ultra-alto.
   - Física: Deformação plástica do Índio sob torque dos parafusos preenchendo as microasperezas.
   - Parâmetros: Torque de aperto T_torque, diâmetro do parafuso, número de parafusos,
     área de contato aparente e presença/ausência de folha de Índio.

2. GÁS DE TROCA DE HÉLIO-4 (HELIUM EXCHANGE GAS):
   - Mecanismo: Condução térmica gasosa em folga (gap) milimétrica sem nenhum contato mecânico.
   - Vantagem crítica: Desacoplamento vibracional total (zero transmissão de vibração da bomba criogênica).
   - Regimes de Rarefação (Número de Knudsen Kn = λ / d):
     * Regime Molecular Livre (Kn > 10, vácuo alto / baixa pressão): Moléculas colidem apenas
       com as paredes; fluxo de calor q ∝ p_gas (Kennard law).
     * Regime Contínuo (Kn < 0.01, alta pressão): Condução clássica de Fourier q = k_He * A * ΔT / d.
     * Regime de Transição (0.01 ≤ Kn ≤ 10): Fórmula de Sherman-Lees interpolando suavemente os regimes.
"""
module ColdInterfaces

using ..Materials
using ..Network

export IndiumBoltedJoint,
       HeliumExchangeGasGap,
       IndiumContactLink,
       ExchangeGasLink,
       clamping_pressure,
       contact_conductance_indium,
       gas_mean_free_path,
       knudsen_number,
       gas_gap_conductance

# =============================================================================
# 1. INTERFACE MECÂNICA: CORDOALHA COM FOLHA DE ÍNDIO (BOLTED JOINT)
# =============================================================================

"""
    IndiumBoltedJoint

Representa a interface mecânica aparafusada com ou sem folha de Índio.

PARÂMETROS:
- `contact_area`: Área aparente de contato das sapatas da cordoalha [m²] (ex: 5 a 15 cm²)
- `num_bolts`: Quantidade de parafusos de aperto (ex: 2 a 4 parafusos)
- `bolt_diameter`: Diâmetro nominal da rosca dos parafusos [m] (ex: M3 = 3e-3 m, M4 = 4e-3 m)
- `torque`: Torque de aperto aplicado em cada parafuso [N·m] (ex: 0.8 N·m para M3, 2.0 N·m para M4)
- `torque_factor`: Fator de atrito da rosca K_t (adimensional, ~0.20 para aço inox lubrificado/padrão)
- `has_indium`: Se `true`, utiliza folha de Índio puro (0.1 mm a 0.5 mm); se `false`, contato metal-metal seco
"""
struct IndiumBoltedJoint
    contact_area::Float64
    num_bolts::Int
    bolt_diameter::Float64
    torque::Float64
    torque_factor::Float64
    has_indium::Bool
    
    function IndiumBoltedJoint(contact_area::Real;
                               num_bolts::Int = 2,
                               bolt_diameter::Real = 3e-3,
                               torque::Real = 0.8,
                               torque_factor::Real = 0.2,
                               has_indium::Bool = true)
        new(Float64(contact_area), num_bolts, Float64(bolt_diameter),
            Float64(torque), Float64(torque_factor), has_indium)
    end
end

"""
    clamping_pressure(joint::IndiumBoltedJoint)

Calcula a pressão média de aperto aparente [Pa] gerada pelos parafusos:
    F_parafuso = Torque / (K_t * d_parafuso)
    Pressão = (N_parafusos * F_parafuso) / Area_contato
"""
function clamping_pressure(joint::IndiumBoltedJoint)
    F_per_bolt = joint.torque / (joint.torque_factor * joint.bolt_diameter) # Força de tração [N]
    total_force = joint.num_bolts * F_per_bolt
    return total_force / joint.contact_area # Pressão em Pascals [N/m²]
end

"""
    contact_conductance_indium(joint::IndiumBoltedJoint, T::Real)

Calcula a condutância térmica de contato h_c [W/(m²·K)] a partir da pressão e temperatura.

FÍSICA DA CORRELAÇÃO (Baseada em dados experimentais do NIST e Ekin, 2006):
- Com folha de Índio:
  O Índio escoa plasticamente acima de P ~ 3 a 5 MPa. A condutância de contato a baixas
  temperaturas (2 K a 30 K) segue uma relação do tipo:
      h_c(P, T) = C_in * (P / P_ref)^0.7 * (T / 4.2)^1.2
  onde C_in ≈ 1800 W/(m²·K) a 4.2 K sob P_ref = 5 MPa.
- Sem folha de Índio (contato seco Cu-Cu):
  As asperezas microscópicas não se conformam plasticamente. A condutância é cerca de 20 a 30 vezes menor:
      h_c,seco(P, T) = C_seco * (P / P_ref)^0.6 * (T / 4.2)^1.0
  com C_seco ≈ 80 W/(m²·K).
"""
function contact_conductance_indium(joint::IndiumBoltedJoint, T::Real)
    P = clamping_pressure(joint)
    P_ref = 5e6 # Pressão de referência: 5 MPa
    T_norm = max(T, 1.0) / 4.2
    P_norm = max(P, 1e5) / P_ref
    
    if joint.has_indium
        # Junta com folha de Índio: excelente conformação plástica
        C_in = 1800.0 # W/(m²·K) a 4.2 K e 5 MPa
        return C_in * (P_norm^0.7) * (T_norm^1.2)
    else
        # Contato seco: asperezas rígidas com vácuo intermediário
        C_dry = 80.0  # W/(m²·K) a 4.2 K e 5 MPa
        return C_dry * (P_norm^0.6) * (T_norm^1.0)
    end
end

"""
    IndiumContactLink(name, node_a, node_b, joint)

Elo térmico dinâmico que calcula a resistência de contato R_c = 1 / (h_c * A) em cada temperatura.
"""
struct IndiumContactLink <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    joint::IndiumBoltedJoint
end

function Network.heat_flow(link::IndiumContactLink, Ta::Real, Tb::Real)
    T_mean = (Ta + Tb) / 2
    hc = contact_conductance_indium(link.joint, T_mean)
    conductance = hc * link.joint.contact_area # Condutância total [W/K]
    return conductance * (Ta - Tb)
end

# =============================================================================
# 2. INTERFACE FLUIDODINÂMICA: GÁS DE TROCA DE HÉLIO-4 (EXCHANGE GAS GAP)
# =============================================================================

# Constante de Boltzmann: k_B = 1.380649 × 10⁻²³ J/K
const K_BOLTZMANN = 1.380649e-23
# Diâmetro cinético molecular do átomo de Hélio: d_He ≈ 2.18 × 10⁻¹⁰ m
const D_HELIUM = 2.18e-10
# Massa de um átomo de Hélio-4: m_He = 4.0026 u ≈ 6.6465 × 10⁻²⁷ kg
const MASS_HE4 = 6.6464764e-27

"""
    HeliumExchangeGasGap

Representa a folga de gás de troca de Hélio entre a superfície fria e a amostra/DAC.

PARÂMETROS:
- `gap_distance`: Distância da folga entre superfícies d [m] (ex: 1 mm a 5 mm)
- `area`: Área da superfície de troca de calor [m²]
- `pressure`: Pressão do gás Hélio na câmara [Pa] (1 mbar = 100 Pa)
- `accommodation_factor`: Coeficiente de acomodação térmica α_eff (0.4 a 0.7 para metais polidos em He)
"""
struct HeliumExchangeGasGap
    gap_distance::Float64
    area::Float64
    pressure::Float64
    accommodation_factor::Float64
    
    function HeliumExchangeGasGap(gap_distance::Real, area::Real;
                                  pressure::Real = 10.0, # 0.1 mbar = 10 Pa padrão
                                  accommodation_factor::Real = 0.5)
        new(Float64(gap_distance), Float64(area), Float64(pressure), Float64(accommodation_factor))
    end
end

"""
    gas_mean_free_path(pressure, T)

Calcula o Livre Caminho Médio (λ) dos átomos de Hélio [m]:
    λ = (k_B * T) / (√2 * π * d_mol² * p)
"""
function gas_mean_free_path(pressure::Real, T::Real)
    p_safe = max(pressure, 1e-8) # evita divisão por zero em vácuo absoluto
    return (K_BOLTZMANN * T) / (sqrt(2) * π * (D_HELIUM^2) * p_safe)
end

"""
    knudsen_number(gap::HeliumExchangeGasGap, T::Real)

Calcula o Número de Knudsen:
    Kn = λ / d
Determina se o gás se comporta como um contínuo clássico (Kn < 0.01) ou gás rarefeito/molecular (Kn > 10).
"""
function knudsen_number(gap::HeliumExchangeGasGap, T::Real)
    lambda = gas_mean_free_path(gap.pressure, T)
    return lambda / gap.gap_distance
end

"""
    gas_gap_conductance(gap::HeliumExchangeGasGap, T::Real)

Calcula o coeficiente de transferência de calor por condução no gás h_gas [W/(m²·K)]
utilizando o modelo unificado de Sherman-Lees (interpolação suave contínuo-molecular).

FÍSICA DOS REGIMES:
1. Regime Contínuo (Kn << 1):
   h_cont = k_He(T) / d
   onde k_He(T) ≈ 2.78e-3 * T^0.7 [W/(m·K)] é a condutividade térmica do gás Hélio.
   (Independe da pressão, pois densidade e livre caminho médio se cancelam mutuamente).

2. Regime Molecular Livre (Kn >> 1):
   As moléculas viajam diretamente de uma parede à outra sem colisões intermoleculares.
   Pela lei de Kennard:
       h_mol = α_eff * ((γ + 1)/(γ - 1)) * √(R_gas / (8π * M * T)) * p
   Para o Hélio (monoatômico, γ = 5/3):
       h_mol ≈ 2.15 * α_eff * (p / √T)

3. Interpolação de Sherman-Lees:
       h_gas = h_cont / (1 + h_cont / h_mol)
   Comporta-se de forma assintoticamente exata em baixas pressões (h ∝ p)
   e satura no limite contínuo em pressões elevadas.
"""
function gas_gap_conductance(gap::HeliumExchangeGasGap, T::Real)
    T_safe = max(T, 1.0)
    p = gap.pressure
    d = gap.gap_distance
    α = gap.accommodation_factor
    
    # 1. Condutância no limite contínuo (Fourier clássico)
    k_He = 2.78e-3 * (T_safe^0.7) # Condutividade térmica do gás He [W/(m·K)]
    h_cont = k_He / d
    
    # 2. Condutância no limite molecular livre (Kennard)
    # Constante para He monoatômico: √(k_B / (2π * m_He))
    molecular_factor = sqrt(K_BOLTZMANN / (2π * MASS_HE4 * T_safe))
    h_mol = α * 2.5 * p * molecular_factor # [W/(m²·K)]
    
    # 3. Interpolação harmônica de Sherman-Lees
    return (h_cont * h_mol) / (h_cont + h_mol)
end

"""
    ExchangeGasLink(name, node_a, node_b, gap)

Elo térmico representando a troca de calor através do gás de troca de Hélio.
"""
struct ExchangeGasLink <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    gap::HeliumExchangeGasGap
end

function Network.heat_flow(link::ExchangeGasLink, Ta::Real, Tb::Real)
    T_mean = (Ta + Tb) / 2
    h_gas = gas_gap_conductance(link.gap, T_mean)
    conductance = h_gas * link.gap.area # Condutância líquida [W/K]
    return conductance * (Ta - Tb)
end

end # module ColdInterfaces
