"""
Módulo Network
==============

Abstração orientada a grafos para redes térmicas e acoplamento termo-mecânico criogênico.

TOPOLOGIA DE UMA REDE TÉRMICA DE PARÂMETROS CONCENTRADOS (LUMPED NETWORK):
O sistema físico é discretizado em dois elementos fundamentais:
1. NÓS TÉRMICOS (`ThermalNode`):
   - Representam os corpos ou subsistemas isotérmicos (Mini-DAC, cordoalha, suportes).
   - Possuem massa M, material, capacidade calorífica C(T) = M * cp(T) e aportes de calor externos.
   - Podem ser 'Livres' (temperatura calculada dinamicamente pelo balanço térmico) ou
     'Fixos' (condições de contorno de Dirichlet, como o dedo frio do criostato a 4.2 K
     ou o primeiro estágio do escudo a 40 K, tratados como reservatórios térmicos infinitos).

2. ELOS TÉRMICOS (`AbstractThermalLink`):
   - Conectam um par de nós (nó A e nó B) e determinam a taxa de transferência de calor Q_{A -> B}.
   - Convenção de sinal: Q > 0 significa calor fluindo de A para B (ou seja, quando Ta > Tb).
   - Tipos físicos de elos:
     * `ConductionLink`: Condução unidimensional por sólidos usando a integral exata de Kirchhoff.
     * `ContactLink`: Resistência de constrição/interface em juntas parafusadas com folha de Índio.
     * `RadiationLink`: Radiação térmica de corpo cinzento via lei de Stefan-Boltzmann (σ * ϵ_eff * A * (Ta⁴ - Tb⁴)).
     * `GenericResistorLink`: Resistência arbitrária R [K/W] para varreduras paramétricas e estudos de trade-off.

BALANÇO DE ENERGIA (PRIMEIRA LEI DA TERMODINÂMICA):
Para cada nó livre i, a equação que governa sua temperatura é:
    M_i * cp_i(T_i) * (dT_i/dt) = Q_externo_i + ∑_{j conectada} Q_{j -> i}(T)
No regime estacionário (steady-state), dT/dt = 0, logo o resíduo líquido de energia deve ser nulo:
    Res_i = Q_externo_i + ∑ Q_{j -> i} = 0
"""
module Network

using LinearAlgebra
using ..Materials

export ThermalNode,
       AbstractThermalLink,
       ConductionLink,
       ContactLink,
       RadiationLink,
       GenericResistorLink,
       ThermalSystem,
       heat_flow,
       stiffness_axial,
       energy_balance_residuals,
       energy_balance_residuals!

"""
    ThermalNode{M<:AbstractMaterial}

Representa um nó da rede térmica (massa concentrada ou reservatório de fronteira).
Parametrizado no tipo do material para evitar instabilidade de tipo.

CAMPOS:
- `name`: Identificador legível (ex: "Mini DAC", "Dedo Frio Criostato 4.2K")
- `mass`: Massa em kg (usar `Inf` para reservatórios térmicos de fronteira)
- `material`: Instância concreta de `M<:AbstractMaterial` (define cp(T), k(T), ρ)
- `heat_load`: Potência térmica externa aplicada diretamente em Watts (ex: feixe síncrotron q_beam = 10 mW)
- `is_fixed`: Se `true`, a temperatura deste nó é fixa (fronteira de Dirichlet)
- `fixed_temp`: Temperatura fixa em Kelvin se `is_fixed == true`
"""
mutable struct ThermalNode{M<:AbstractMaterial}
    name::String
    mass::Float64
    material::M
    heat_load::Float64
    is_fixed::Bool
    fixed_temp::Float64
    
    function ThermalNode(name::String, mass::Real, material::M;
                         heat_load::Real = 0.0, is_fixed::Bool = false, fixed_temp::Real = 0.0) where {M<:AbstractMaterial}
        new{M}(name, Float64(mass), material, Float64(heat_load), is_fixed, Float64(fixed_temp))
    end
end

"""
Super-tipo abstrato para qualquer elemento condutor de calor entre dois nós.
"""
abstract type AbstractThermalLink end

# Constante universal de Stefan-Boltzmann: σ = 5.670374419 × 10⁻⁸ W / (m² · K⁴)
const SIGMA_SB = 5.670374419e-8

# =============================================================================
# 1. ELO DE CONDUÇÃO SÓLIDA 1D (COM RIGIDEZ MECÂNICA ACOPLADA)
# =============================================================================
"""
    ConductionLink{M<:AbstractMaterial}(name, node_a, node_b, area, length, material)

Elo de condução sólida unidimensional (hastes de suporte, cordoalhas térmicas, tubos).
Parametrizado no tipo do material para garantir despacho estático e zero boxing.

CO-DESIGN TERMO-MECÂNICO:
Este elo unifica as duas físicas concorrentes do projeto da nanoestação EMA:
1. Térmica: O calor é calculado pela integral exata de Kirchhoff:
       Q_{a -> b} = (A / L) * ∫_{Tb}^{Ta} k(T) dT
2. Mecânica: A rigidez elástica axial [N/m] na temperatura média é dada por:
       k_axial = E(T_médio) * A / L
"""
struct ConductionLink{M<:AbstractMaterial} <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    area::Float64      # Área da seção transversal [m²]
    length::Float64    # Comprimento condutivo [m]
    material::M
    
    function ConductionLink(name::String, node_a_idx::Int, node_b_idx::Int,
                            area::Real, length::Real, material::M) where {M<:AbstractMaterial}
        new{M}(name, node_a_idx, node_b_idx, Float64(area), Float64(length), material)
    end
end

function heat_flow(link::ConductionLink{M}, Ta::Real, Tb::Real) where {M<:AbstractMaterial}
    geom_factor = link.area / link.length
    return geom_factor * thermal_conductivity_integral(link.material, Tb, Ta)
end

function stiffness_axial(link::ConductionLink{M}, T_mean::Real) where {M<:AbstractMaterial}
    E = youngs_modulus(link.material, T_mean)
    return E * link.area / link.length # [N/m]
end

# =============================================================================
# 2. ELO DE RESISTÊNCIA DE CONTATO MECÂNICO
# =============================================================================
"""
    ContactLink(name, node_a, node_b, resistance)

Resistência térmica de interface entre dois componentes prensados ou parafusados.
"""
struct ContactLink <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    resistance::Float64 # Resistência de contato R_c [K/W]
    
    function ContactLink(name::String, node_a_idx::Int, node_b_idx::Int, resistance::Real)
        new(name, node_a_idx, node_b_idx, Float64(resistance))
    end
end

function heat_flow(link::ContactLink, Ta::Real, Tb::Real)
    return (Ta - Tb) / link.resistance
end

# =============================================================================
# 3. ELO DE RADIAÇÃO TÉRMICA (STEFAN-BOLTZMANN)
# =============================================================================
"""
    RadiationLink(name, node_a, node_b, area, emissivity_eff)

Transferência de calor por radiação entre superfícies cinzentas no vácuo:
    Q_{a -> b} = σ * ϵ_eff * A * (Ta⁴ - Tb⁴)
"""
struct RadiationLink <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    area::Float64           # Área radiativa de troca [m²]
    emissivity_eff::Float64 # Emissividade mútua efetiva (fator de forma e acabamento superficial)
    
    function RadiationLink(name::String, node_a_idx::Int, node_b_idx::Int, area::Real, emissivity_eff::Real)
        new(name, node_a_idx, node_b_idx, Float64(area), Float64(emissivity_eff))
    end
end

function heat_flow(link::RadiationLink, Ta::Real, Tb::Real)
    return SIGMA_SB * link.emissivity_eff * link.area * (Ta^4 - Tb^4)
end

# =============================================================================
# 4. ELO RESISTIVO GENÉRICO (PARA VARREDURAS PARAMÉTRICAS)
# =============================================================================
"""
    GenericResistorLink(name, node_a, node_b, resistance)

Elo linear genérico Q = (Ta - Tb)/R.
"""
struct GenericResistorLink <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    resistance::Float64
    
    function GenericResistorLink(name::String, node_a_idx::Int, node_b_idx::Int, resistance::Real)
        new(name, node_a_idx, node_b_idx, Float64(resistance))
    end
end

function heat_flow(link::GenericResistorLink, Ta::Real, Tb::Real)
    return (Ta - Tb) / link.resistance
end

# =============================================================================
# 5. MONTAGEM DO SISTEMA E BALANÇO DE ENERGIA
# =============================================================================
"""
    ThermalSystem(nodes, links)

Agrega o grafo completo de nós térmicos e elos de ligação do instrumento.
"""
struct ThermalSystem
    nodes::Vector{ThermalNode}
    links::Vector{AbstractThermalLink}
    
    function ThermalSystem(nodes::AbstractVector{<:ThermalNode}, links::AbstractVector = AbstractThermalLink[])
        new(Vector{ThermalNode}(nodes), Vector{AbstractThermalLink}(links))
    end
end

"""
    energy_balance_residuals!(residuals, sys, T_full)

Versão in-place pré-alocada do balanço de energia nos nós térmicos.
"""
function energy_balance_residuals!(residuals::AbstractVector, sys::ThermalSystem, T_full::AbstractVector)
    n = length(sys.nodes)
    for i in 1:n
        residuals[i] = sys.nodes[i].heat_load
    end
    
    for link in sys.links
        idx_a = link.node_a_idx
        idx_b = link.node_b_idx
        Ta = T_full[idx_a]
        Tb = T_full[idx_b]
        
        Q = heat_flow(link, Ta, Tb)
        residuals[idx_a] -= Q
        residuals[idx_b] += Q
    end
    
    return residuals
end

"""
    energy_balance_residuals(sys, T_full)

Calcula o vetor de resíduos de conservação de energia para todas as massas do sistema.
"""
function energy_balance_residuals(sys::ThermalSystem, T_full::AbstractVector)
    residuals = zeros(eltype(T_full), length(sys.nodes))
    return energy_balance_residuals!(residuals, sys, T_full)
end

end # module Network
