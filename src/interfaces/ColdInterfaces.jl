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

FÍSICA DA CORRELAÇÃO:
- Em temperaturas criogênicas (T <= 25 K), segue a correlação de Ekin / NIST:
      h_c(P, T) = C_in * (P / P_ref)^0.7 * (T / 4.2)^1.2
- Em temperaturas intermediárias e ambiente (T > 25 K até 300 K), a resistência interfacial
  é regularizada fisicamente com saturação assintótica suave (C1 contínua), evitando a
  divergência não-física de T^1.2 que geraria centenas de milhares de W/(m²·K) a 300 K.
- Para T < 1.0 K (regime sub-Kelvin), decai linearmente respeitando a Terceira Lei (h_c -> 0 quando T -> 0).
"""
function contact_conductance_indium(joint::IndiumBoltedJoint, T::Real)
    P = clamping_pressure(joint)
    P_ref = 5e6 # Pressão de referência: 5 MPa
    P_norm = max(P, 1e5) / P_ref
    
    # Salvaguarda sub-Kelvin: decaimento linear proporcional a T com sensibilidade estrita
    if T < 1.0
        T_safe = ifelse(T >= zero(T), T, zero(T))
        hc_1 = contact_conductance_indium(joint, 1.0)
        return hc_1 * T_safe
    end
    
    # Ponto de transição do regime criogênico (Ekin) para saturação plástica macroscópica
    T_trans = 25.0
    
    if joint.has_indium
        C_in = 1800.0 # W/(m²·K) a 4.2 K e 5 MPa
        h_base = C_in * (P_norm^0.7)
        if T <= T_trans
            return h_base * ((T / 4.2)^1.2)
        else
            # Saturação suave C1 contínua para T > 25 K (evita explosão em 300 K)
            h_trans = h_base * ((T_trans / 4.2)^1.2)
            dh_dT = 1.2 * h_trans / T_trans
            h_max = 1.6 * h_trans
            delta_h = h_max - h_trans
            k_decay = dh_dT / delta_h
            return h_trans + delta_h * (1.0 - exp(-k_decay * (T - T_trans)))
        end
    else
        C_dry = 80.0  # W/(m²·K) a 4.2 K e 5 MPa
        h_base = C_dry * (P_norm^0.6)
        if T <= T_trans
            return h_base * (T / 4.2)
        else
            h_trans = h_base * (T_trans / 4.2)
            dh_dT = h_trans / T_trans
            h_max = 1.6 * h_trans
            delta_h = h_max - h_trans
            k_decay = dh_dT / delta_h
            return h_trans + delta_h * (1.0 - exp(-k_decay * (T - T_trans)))
        end
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
    T_mean = (Ta + Tb) * 0.5
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
    T_safe = ifelse(T >= zero(T), T, zero(T))
    return (K_BOLTZMANN * T_safe) / (sqrt(2) * π * (D_HELIUM^2) * p_safe)
end

"""
    knudsen_number(gap::HeliumExchangeGasGap, T::Real)

Calcula o Número de Knudsen:
    Kn = λ / d
"""
function knudsen_number(gap::HeliumExchangeGasGap, T::Real)
    lambda = gas_mean_free_path(gap.pressure, T)
    return lambda / gap.gap_distance
end

"""
    gas_gap_conductance(gap::HeliumExchangeGasGap, T::Real)

Calcula o coeficiente de transferência de calor por condução no gás h_gas [W/(m²·K)]
utilizando o modelo unificado de Sherman-Lees (interpolação suave contínuo-molecular).
"""
function gas_gap_conductance(gap::HeliumExchangeGasGap, T::Real)
    if T < 1.0
        h_1 = gas_gap_conductance(gap, 1.0)
        return h_1 * ifelse(T >= zero(T), T, zero(T))
    end
    
    p = gap.pressure
    d = gap.gap_distance
    α = gap.accommodation_factor
    
    # 1. Condutância no limite contínuo (Fourier clássico)
    k_He = 2.78e-3 * (T^0.7) # Condutividade térmica do gás He [W/(m·K)]
    h_cont = k_He / d
    
    # 2. Condutância no limite molecular livre (Kennard)
    molecular_factor = sqrt(K_BOLTZMANN / (2π * MASS_HE4 * T))
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
    T_mean = (Ta + Tb) * 0.5
    h_gas = gas_gap_conductance(link.gap, T_mean)
    conductance = h_gas * link.gap.area # Condutância líquida [W/K]
    return conductance * (Ta - Tb)
end

end # module ColdInterfaces
