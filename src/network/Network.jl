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
       energy_balance_residuals

"""
    ThermalNode

Representa um nó da rede térmica (massa concentrada ou reservatório de fronteira).

CAMPOS:
- `name`: Identificador legível (ex: "Mini DAC", "Dedo Frio Criostato 4.2K")
- `mass`: Massa em kg (usar `Inf` para reservatórios térmicos de fronteira)
- `material`: Instância de `AbstractMaterial` (define cp(T), k(T), ρ)
- `heat_load`: Potência térmica externa aplicada diretamente em Watts (ex: feixe síncrotron q_beam = 10 mW)
- `is_fixed`: Se `true`, a temperatura deste nó é fixa (fronteira de Dirichlet)
- `fixed_temp`: Temperatura fixa em Kelvin se `is_fixed == true`
"""
mutable struct ThermalNode
    name::String
    mass::Float64
    material::AbstractMaterial
    heat_load::Float64
    is_fixed::Bool
    fixed_temp::Float64
    
    function ThermalNode(name::String, mass::Real, material::AbstractMaterial;
                         heat_load::Real = 0.0, is_fixed::Bool = false, fixed_temp::Real = 0.0)
        new(name, Float64(mass), material, Float64(heat_load), is_fixed, Float64(fixed_temp))
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
    ConductionLink(name, node_a, node_b, area, length, material)

Elo de condução sólida unidimensional (hastes de suporte, cordoalhas térmicas, tubos).

CO-DESIGN TERMO-MECÂNICO:
Este elo unifica as duas físicas concorrentes do projeto da nanoestação EMA:
1. Térmica: O calor é calculado pela integral exata de Kirchhoff:
       Q_{a -> b} = (A / L) * ∫_{Tb}^{Ta} k(T) dT
2. Mecânica: A rigidez elástica axial [N/m] na temperatura média é dada por:
       k_axial = E(T_médio) * A / L
   Isso permite que ao alterar a área A ou o comprimento L, tanto a condutância
   térmica quanto a frequência natural da estrutura sejam atualizadas em conjunto.
"""
struct ConductionLink <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    area::Float64      # Área da seção transversal [m²]
    length::Float64    # Comprimento condutivo [m]
    material::AbstractMaterial
end

function heat_flow(link::ConductionLink, Ta::Real, Tb::Real)
    # Integral de Tb até Ta: se Ta > Tb, a integral é positiva (calor flui de a para b)
    geom_factor = link.area / link.length
    return geom_factor * thermal_conductivity_integral(link.material, Tb, Ta)
end

function stiffness_axial(link::ConductionLink, T_mean::Real)
    E = youngs_modulus(link.material, T_mean)
    return E * link.area / link.length # [N/m]
end

# =============================================================================
# 2. ELO DE RESISTÊNCIA DE CONTATO MECÂNICO
# =============================================================================
"""
    ContactLink(name, node_a, node_b, resistance)

Resistência térmica de interface entre dois componentes prensados ou parafusados.

FÍSICA DO CONTATO CRIOGÊNICO:
No vácuo de uma linha de luz síncrotron, não há convecção de ar entre as peças.
O contato físico ocorre apenas no topo das micro-asperezas superficiais microscópicas
(área real de contato < 1% da área geométrica aparente).
- Em contato metal-metal seco a 4 K, a resistência de constrição é enorme (dezenas de K/W).
- Para contornar isso, utiliza-se uma folha fina de Índio puro (Indium foil, 0.1 mm a 0.5 mm)
  entre as faces. O índio tem baixíssima tensão de escoamento mesmo a frio, deformando-se
  plasticamente sob aperto e preenchendo as cavidades microscópicas, reduzindo a resistência R_c.
"""
struct ContactLink <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    resistance::Float64 # Resistência de contato R_c [K/W]
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

ONDE SE APLICA NA EMA:
Entre o escudo térmico (~40 K) e a Mini-DAC (~4 K), ou entre a câmara externa de vácuo (300 K)
e o escudo térmico. Devido à dependência de T⁴, a 4 K a radiação costuma ser ordens de grandeza
menor do que a condução, mas entre 300 K e 40 K ela representa uma carga de calor crítica.
"""
struct RadiationLink <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    area::Float64           # Área radiativa de troca [m²]
    emissivity_eff::Float64 # Emissividade mútua efetiva (fator de forma e acabamento superficial)
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
Utilizado para varrer faixas de resistência térmica desconhecidas durante a fase conceitual.
"""
struct GenericResistorLink <: AbstractThermalLink
    name::String
    node_a_idx::Int
    node_b_idx::Int
    resistance::Float64
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
end

"""
    energy_balance_residuals(sys, T_full)

Calcula o vetor de resíduos de conservação de energia para todas as massas do sistema.

EQUAÇÃO DO BALANÇO:
Para cada nó i:
    Resíduo_i = Q_externo_i + ∑_{elo conectado} Q_entrando - ∑_{elo conectado} Q_saindo

Fisicamente:
- Se Resíduo_i > 0: Há excesso líquido de potência entrando no nó i (temperatura tenderia a subir).
- Se Resíduo_i < 0: Há déficit líquido de potência (o calor que sai é maior do que o que entra).
- Se Resíduo_i == 0: O nó está em equilíbrio térmico estacionário exato.
"""
function energy_balance_residuals(sys::ThermalSystem, T_full::AbstractVector)
    n = length(sys.nodes)
    residuals = zeros(eltype(T_full), n)
    
    # 1. Aportes externos diretos (ex: feixe síncrotron na amostra)
    for i in 1:n
        residuals[i] += sys.nodes[i].heat_load
    end
    
    # 2. Somatório das trocas térmicas em cada elo do grafo
    for link in sys.links
        idx_a = link.node_a_idx
        idx_b = link.node_b_idx
        Ta = T_full[idx_a]
        Tb = T_full[idx_b]
        
        # Fluxo líquido que parte de A em direção a B
        Q = heat_flow(link, Ta, Tb)
        
        residuals[idx_a] -= Q  # Energia sai do nó A
        residuals[idx_b] += Q  # Energia entra no nó B
    end
    
    return residuals
end

end # module Network
