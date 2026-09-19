"""
Estudo Comparativo das Interfaces do Dedo Frio: Braid com Índio vs Gás de Hélio
==============================================================================

Este script gera uma figura vetorial com CairoMakie detalhando:
1. Painel A: A mecânica de contato por aperto de parafusos com folha de Índio vs contato seco.
   Mostra o ganho de > 20x na condutância de contato e a redução drástica na temperatura da DAC.
2. Painel B: A física da folga de gás de troca de Hélio (He-4 exchange gas).
   Mapeia a transição entre o regime molecular livre (baixa pressão) e o contínuo (alta pressão),
   demonstrando o funcionamento como chave térmica ajustável sem acoplamento mecânico.
"""

using CairoMakie
using LinearAlgebra

# Carrega o framework
include(joinpath(@__DIR__, "..", "src", "CryoThermal.jl"))
using .CryoThermal

println("Iniciando estudo detalhado das interfaces frias...")

# Materiais
cu   = CopperOFHC()
ss   = StainlessSteel304()
cube = BerylliumCopper()

area_contact = 10e-4 # 10 cm² de área de sapata
q_beam       = 10e-3 # 10 mW

# =============================================================================
# SIMULAÇÃO 1: VARREDURA DE TORQUE / PRESSÃO COM E SEM ÍNDIO
# =============================================================================
torques = range(0.1, 2.0, length=50) # Torque de 0.1 N·m a 2.0 N·m (2 parafusos M3)

P_pressures = Float64[]
Rc_indium   = Float64[]
Rc_dry      = Float64[]
T_dac_in    = Float64[]
T_dac_dry   = Float64[]

for tau in torques
    joint_in  = IndiumBoltedJoint(area_contact; num_bolts=2, bolt_diameter=3e-3, torque=tau, has_indium=true)
    joint_dry = IndiumBoltedJoint(area_contact; num_bolts=2, bolt_diameter=3e-3, torque=tau, has_indium=false)
    
    P = clamping_pressure(joint_in) / 1e6 # em MPa
    push!(P_pressures, P)
    
    hc_in  = contact_conductance_indium(joint_in, 4.2)
    hc_dry = contact_conductance_indium(joint_dry, 4.2)
    
    push!(Rc_indium, 1.0 / (hc_in * area_contact))
    push!(Rc_dry,    1.0 / (hc_dry * area_contact))
    
    # Resolver a temperatura de equilíbrio da DAC para cada caso
    for (has_in, t_dest) in [(true, T_dac_in), (false, T_dac_dry)]
        j = IndiumBoltedJoint(area_contact; num_bolts=2, bolt_diameter=3e-3, torque=tau, has_indium=has_in)
        
        n1 = ThermalNode("ColdHead", Inf, cu; is_fixed=true, fixed_temp=4.2)
        n2 = ThermalNode("Braid", 0.134, cu)
        n3 = ThermalNode("DAC", 0.0513, cube; heat_load=q_beam)
        n4 = ThermalNode("Support", 0.332, ss)
        n5 = ThermalNode("Shield", Inf, ss; is_fixed=true, fixed_temp=40.0)
        
        l1 = IndiumContactLink("BraidContact", 1, 2, j)
        l2 = ConductionLink("BraidCond", 2, 3, 12e-6, 60e-3, cu)
        l3 = ConductionLink("SupportStrut", 3, 4, 15e-6, 25e-3, ss)
        l4 = ConductionLink("SupportShield", 4, 5, 25e-6, 20e-3, ss)
        l5 = RadiationLink("Rad", 5, 3, 10e-4, 0.05)
        
        sys = ThermalSystem([n1, n2, n3, n4, n5], [l1, l2, l3, l4, l5])
        res = solve_steady_state(sys; tol=1e-6)
        push!(t_dest, res.temperatures[3])
    end
end

# =============================================================================
# SIMULAÇÃO 2: VARREDURA DE PRESSÃO DO GÁS DE TROCA DE HÉLIO
# =============================================================================
pressures_mbar = exp10.(range(-5, 2, length=60)) # 10⁻⁵ mbar até 100 mbar
h_gas_vals     = Float64[]
Kn_vals        = Float64[]
T_dac_gas      = Float64[]

gap_dist = 2e-3 # folga de 2 mm

for p_mbar in pressures_mbar
    p_pa = p_mbar * 100.0 # 1 mbar = 100 Pa
    gap = HeliumExchangeGasGap(gap_dist, area_contact; pressure=p_pa)
    
    h = gas_gap_conductance(gap, 4.2)
    Kn = knudsen_number(gap, 4.2)
    
    push!(h_gas_vals, h)
    push!(Kn_vals, Kn)
    
    # Resolve sistema substituindo a cordoalha pelo gás de troca
    n1 = ThermalNode("ColdHead", Inf, cu; is_fixed=true, fixed_temp=4.2)
    n3 = ThermalNode("DAC", 0.0513, cube; heat_load=q_beam)
    n4 = ThermalNode("Support", 0.332, ss)
    n5 = ThermalNode("Shield", Inf, ss; is_fixed=true, fixed_temp=40.0)
    
    # Nós no sys_gas:
    # 1: ColdHead (Fixo 4.2 K)
    # 2: DAC (Livre)
    # 3: Support (Livre)
    # 4: Shield (Fixo 40.0 K)
    l_gas = ExchangeGasLink("HeliumGap", 1, 2, gap) # Troca direta dedo frio (1) -> DAC (2)
    l3 = ConductionLink("SupportStrut", 2, 3, 15e-6, 25e-3, ss)  # DAC (2) -> Support (3)
    l4 = ConductionLink("SupportShield", 3, 4, 25e-6, 20e-3, ss) # Support (3) -> Shield (4)
    l5 = RadiationLink("Rad", 4, 2, 10e-4, 0.05)                # Shield (4) -> DAC (2)
    
    sys_gas = ThermalSystem([n1, n3, n4, n5], [l_gas, l3, l4, l5])
    res_gas = solve_steady_state(sys_gas; tol=1e-6)
    push!(T_dac_gas, res_gas.temperatures[2]) # nó 2 é a DAC
end

println("Simulações concluídas. Renderizando painel duplo com CairoMakie...")

# =============================================================================
# GERAÇÃO DA FIGURA DUPLA COM CAIROMAKIE
# =============================================================================
fig = Figure(size = (1100, 750), fontsize = 13)

# -----------------------------------------------------------------------------
# PAINEL 1: JUNTA APARAFUSADA (COM vs SEM ÍNDIO)
# -----------------------------------------------------------------------------
ax1 = Axis(fig[1, 1],
    title = "A) Interface Mecânica: Cordoalha com Folha de Índio vs Contato Seco",
    xlabel = "Pressão Média de Aperto dos Parafusos [MPa]",
    ylabel = "Temperatura da Mini-DAC T_DAC [K]",
    limits = ((0.3, 6.8), (4.0, 15.0))
)

lines!(ax1, P_pressures, T_dac_dry, color = :red, linewidth = 2.5,
       label = "Contato Seco Cu-Cu (Sem Índio)")
lines!(ax1, P_pressures, T_dac_in, color = :blue, linewidth = 2.5,
       label = "Contato com Folha de Índio (t = 0.2 mm)")

hlines!(ax1, [6.0], color = :darkred, linestyle = :dash, linewidth = 1.5)
text!(ax1, 0.5, 6.3, text = "Meta de Operação Criogênica (T_DAC ≤ 6 K)", color = :darkred, fontsize = 11)

# Zona de escoamento plástico do Índio (P > 3 MPa)
vspan!(ax1, 3.0, 6.8, color = (:blue, 0.08))
text!(ax1, 3.2, 13.5, text = "Zona de Escoamento Plástico do Índio\n(P ≥ 3 MPa: vedação microscópica)", color = :navy, fontsize = 11)

axislegend(ax1, position = :rt, framevisible = true, bgcolor = (:white, 0.9))

# -----------------------------------------------------------------------------
# PAINEL 2: GÁS DE TROCA DE HÉLIO (TRANSIÇÃO DE KNUDSEN)
# -----------------------------------------------------------------------------
ax2 = Axis(fig[2, 1],
    title = "B) Interface Térmica sem Contato Mecânico: Gás de Troca de Hélio (Gap = 2 mm)",
    xlabel = "Pressão do Gás Hélio [mbar]",
    ylabel = "Temperatura da Mini-DAC T_DAC [K]",
    xscale = log10,
    limits = ((1e-5, 1e2), (4.0, 35.0))
)

lines!(ax2, pressures_mbar, T_dac_gas, color = :purple, linewidth = 2.5,
       label = "T_DAC vs Pressão de Hélio")

# Delimitação dos Regimes Físicos de Knudsen
vlines!(ax2, [1.5e-3], color = :gray, linestyle = :dot, linewidth = 1.5)
vlines!(ax2, [1.5],    color = :gray, linestyle = :dot, linewidth = 1.5)

text!(ax2, 1.5e-5, 26.0, text = "REGIME MOLECULAR LIVRE\n(Kn > 10, Alto Vácuo)\nq ∝ p_gas\n[Chave Térmica Aberta]", color = :darkgray, fontsize = 10)
text!(ax2, 5e-2, 26.0, text = "TRANSIÇÃO\n(Sherman-Lees)", color = :darkgray, fontsize = 10)
text!(ax2, 4.0, 26.0, text = "REGIME CONTÍNUO\n(Kn < 0.01)\nq = k_He·A·ΔT/d\n[Chave Fechada]", color = :darkgray, fontsize = 10)

hlines!(ax2, [6.0], color = :darkred, linestyle = :dash, linewidth = 1.5)
text!(ax2, 2.0, 7.5, text = "T_DAC ≤ 6 K (Requer p > 2 mbar)", color = :darkred, fontsize = 11)

axislegend(ax2, position = :rt, framevisible = true, bgcolor = (:white, 0.9))

# Salva a figura
figures_dir = joinpath(@__DIR__, "..", "figures")
mkpath(figures_dir)
output_path = joinpath(figures_dir, "cold_interfaces_study.png")

save(output_path, fig, px_per_unit = 2)
println("✓ Gráfico salvo com sucesso em: ", output_path)
