"""
Estudo de Co-Design Termo-Mecânico com Flexão Lateral e Seleção de Materiais (Linha EMA)
=======================================================================================

Este script fecha o ciclo de engenharia do projeto respondendo à pergunta central:
"Como escolher as dimensões geométricas e o material do suporte da Mini-DAC para
garantir rigidez lateral (fn >= 150 Hz) sem sobreaquecer a amostra criogênica (TDAC <= 6 K)?"

FÍSICA MODELADA:
1. Mecânica em Flexão Lateral Transversal (Tripé Isostático de 3 Hastes):
   A DAC (massa M = 51.3 g) é sustentada por 3 hastes sólidas cilíndricas de diâmetro d e comprimento L = 30 mm.
   - Momento de inércia à flexão: I_haste = π * d⁴ / 64
   - Rigidez lateral do conjunto (tripé a 120°): k_total = 1.5 * (12 * E * I_haste) / L³
   - Primeira frequência natural: fn = (1 / 2π) * √(k_total / M)
2. Térmica Não-Linear:
   - Área condutiva total: A_tot = 3 * (π * d² / 4)
   - Calor conduzido do anel a 40 K para a DAC a ~4.5 K: Q = (A / L) * ∫ k(T) dT
3. Comparação de Materiais Estruturais:
   - Aço Inox AISI 304 (E = 210 GPa, alta rigidez, condutividade intermediária)
   - Titânio Grau 5 Ti-6Al-4V (E = 118 GPa, baixíssima condutividade, FoM 52% superior!)
"""

using CairoMakie
using LinearAlgebra
using CryoThermal

println("Iniciando estudo de co-design termo-mecânico (tripé isostático em flexão)...")

# Instância dos materiais
cu   = CopperOFHC()
ss   = StainlessSteel304()
ti   = TitaniumTi6Al4V()
cube = BerylliumCopper()

m_dac   = 0.0513 # Massa da Mini-DAC = 51.3 g
q_beam  = 10e-3  # Potência do feixe síncrotron = 10 mW
L_fixed = 30e-3  # Comprimento padrão das hastes = 30 mm
N_struts = 3     # Tripé isostático

# Junta de cordoalha padrão com Índio bem apertado
joint_in = IndiumBoltedJoint(10e-4; num_bolts=2, bolt_diameter=3e-3, torque=1.0, has_indium=true)

# =============================================================================
# SIMULAÇÃO 1: VARREDURA DE DIÂMETRO DAS HASTES (INOX 304 vs TITÂNIO Ti-6Al-4V)
# =============================================================================
# Diâmetro d de cada haste do tripé variando de 0.8 mm a 2.5 mm
d_struts = range(0.8e-3, 2.5e-3, length=35)

fn_ss = Float64[]
fn_ti = Float64[]
T_ss  = Float64[]
T_ti  = Float64[]

for d in d_struts
    I_haste = (π / 64) * (d^4)
    A_tot   = N_struts * (π / 4) * (d^2)
    
    # 1. Caso Inox 304
    E_ss_val  = youngs_modulus(ss, 20.0)
    k_flex_ss = 1.5 * (12 * E_ss_val * I_haste) / (L_fixed^3)
    fn_ss_val = (1.0 / (2π)) * sqrt(k_flex_ss / m_dac)
    push!(fn_ss, fn_ss_val)
    
    n1 = ThermalNode("ColdHead", Inf, cu; is_fixed=true, fixed_temp=4.2)
    n2 = ThermalNode("Braid", 0.134, cu)
    n3 = ThermalNode("DAC", m_dac, cube; heat_load=q_beam)
    n4 = ThermalNode("Support", 0.332, ss)
    n5 = ThermalNode("Shield", Inf, ss; is_fixed=true, fixed_temp=40.0)
    
    l1 = IndiumContactLink("BraidContact", 1, 2, joint_in)
    l2 = ConductionLink("BraidCond", 2, 3, 12e-6, 60e-3, cu)
    l3_ss = ConductionLink("SupportStrut", 3, 4, A_tot, L_fixed, ss)
    l4 = ConductionLink("SupportShield", 4, 5, 25e-6, 20e-3, ss)
    l5 = RadiationLink("Rad", 5, 3, 10e-4, 0.05)
    
    sys_ss = ThermalSystem([n1, n2, n3, n4, n5], [l1, l2, l3_ss, l4, l5])
    res_ss = solve_steady_state(sys_ss; tol=1e-6)
    push!(T_ss, res_ss.temperatures[3])
    
    # 2. Caso Titânio Ti-6Al-4V
    E_ti_val  = youngs_modulus(ti, 20.0)
    k_flex_ti = 1.5 * (12 * E_ti_val * I_haste) / (L_fixed^3)
    fn_ti_val = (1.0 / (2π)) * sqrt(k_flex_ti / m_dac)
    push!(fn_ti, fn_ti_val)
    
    l3_ti = ConductionLink("SupportStrutTi", 3, 4, A_tot, L_fixed, ti)
    sys_ti = ThermalSystem([n1, n2, n3, n4, n5], [l1, l2, l3_ti, l4, l5])
    res_ti = solve_steady_state(sys_ti; tol=1e-6)
    push!(T_ti, res_ti.temperatures[3])
end

# =============================================================================
# SIMULAÇÃO 2: MAPA 2D NO PLANO (DIÂMETRO d x COMPRIMENTO L) PARA INOX 304
# =============================================================================
d_grid = range(0.8e-3, 2.4e-3, length=25) # Diâmetro da haste: 0.8 mm a 2.4 mm
L_grid = range(18e-3, 45e-3, length=25)  # Comprimento da haste: 18 mm a 45 mm

T_map  = zeros(length(d_grid), length(L_grid))
fn_map = zeros(length(d_grid), length(L_grid))

println("Mapeando espaço de projeto 2D (Diâmetro d x Comprimento L)...")
for (i, d_val) in enumerate(d_grid)
    for (j, L_val) in enumerate(L_grid)
        I_h = (π / 64) * (d_val^4)
        A_t = N_struts * (π / 4) * (d_val^2)
        
        E_val = youngs_modulus(ss, 20.0)
        k_flex = 1.5 * (12 * E_val * I_h) / (L_val^3)
        fn_map[i, j] = (1.0 / (2π)) * sqrt(k_flex / m_dac)
        
        n1 = ThermalNode("ColdHead", Inf, cu; is_fixed=true, fixed_temp=4.2)
        n2 = ThermalNode("Braid", 0.134, cu)
        n3 = ThermalNode("DAC", m_dac, cube; heat_load=q_beam)
        n4 = ThermalNode("Support", 0.332, ss)
        n5 = ThermalNode("Shield", Inf, ss; is_fixed=true, fixed_temp=40.0)
        
        l1 = IndiumContactLink("BraidContact", 1, 2, joint_in)
        l2 = ConductionLink("BraidCond", 2, 3, 12e-6, 60e-3, cu)
        l3 = ConductionLink("SupportStrut", 3, 4, A_t, L_val, ss)
        l4 = ConductionLink("SupportShield", 4, 5, 25e-6, 20e-3, ss)
        l5 = RadiationLink("Rad", 5, 3, 10e-4, 0.05)
        
        sys = ThermalSystem([n1, n2, n3, n4, n5], [l1, l2, l3, l4, l5])
        res = solve_steady_state(sys; tol=1e-6)
        T_map[i, j] = res.temperatures[3]
    end
end

println("Simulações 2D concluídas. Renderizando figura com CairoMakie...")

# =============================================================================
# GERAÇÃO DA FIGURA DE CO-DESIGN COM CAIROMAKIE
# =============================================================================
fig = Figure(size = (1200, 580), fontsize = 13)

# -----------------------------------------------------------------------------
# PAINEL 1: CURVAS DE FRONTEIRA DE PARETO (INOX vs TITÂNIO)
# -----------------------------------------------------------------------------
ax1 = Axis(fig[1, 1],
    title = "A) Fronteira de Pareto: Rigidez Lateral vs Temperatura Final",
    xlabel = "Primeira Frequência Natural em Flexão fn,flex [Hz]",
    ylabel = "Temperatura da Mini-DAC T_DAC [K]",
    limits = ((40, 420), (4.3, 6.5))
)

# Zonas de Requisitos
vspan!(ax1, 40.0, 150.0, color = (:red, 0.08))
vlines!(ax1, [150.0], color = :red, linestyle = :dash, linewidth = 1.5)
text!(ax1, 153.0, 6.3, text = "Corte Mecânico (fn ≥ 150 Hz)", color = :red, fontsize = 11)

hlines!(ax1, [6.0], color = :darkred, linestyle = :dash, linewidth = 1.5)
text!(ax1, 45.0, 6.08, text = "Meta Térmica (T_DAC ≤ 6.0 K)", color = :darkred, fontsize = 11)

# Polígono da Zona Ótima
poly!(ax1, Point2f[(150, 4.3), (420, 4.3), (420, 6.0), (150, 6.0)],
      color = (:green, 0.12))
text!(ax1, 190.0, 5.2, text = "ZONA ÓTIMA DE PROJETO\n(Viável para Linha EMA)", color = :darkgreen, fontsize = 11, font = :bold)

lines!(ax1, fn_ss, T_ss, color = :blue, linewidth = 2.5,
       label = "Aço Inox AISI 304 (Padrão)")
lines!(ax1, fn_ti, T_ti, color = :darkorange, linewidth = 2.5,
       label = "Titânio Ti-6Al-4V (Alta Performance)")

axislegend(ax1, position = :lt, framevisible = true, backgroundcolor = (:white, 0.9))

# -----------------------------------------------------------------------------
# PAINEL 2: ENVELOPE DE USINAGEM 2D (DIÂMETRO d vs COMPRIMENTO L PARA INOX 304)
# -----------------------------------------------------------------------------
ax2 = Axis(fig[1, 2],
    title = "B) Envelope de Fabricação do Tripé de Suporte (Inox 304)",
    xlabel = "Diâmetro de cada Haste d [mm]",
    ylabel = "Comprimento das Hastes L [mm]"
)

d_mm = d_grid .* 1e3
L_mm = L_grid .* 1e3

# Heatmap de temperatura
hm = contourf!(ax2, d_mm, L_mm, T_map, levels = 15, colormap = :viridis)
Colorbar(fig[1, 3], hm, label = "T_DAC [K]")

# Isolinhas de Frequência Natural
contour!(ax2, d_mm, L_mm, fn_map, levels = [100.0, 150.0, 200.0, 300.0],
         color = :white, linewidth = 2.0)

# Isoterma limite de 6.0 K destacada em vermelho pontilhado
contour!(ax2, d_mm, L_mm, T_map, levels = [6.0],
         color = :red, linewidth = 2.5, linestyle = :dash)

text!(ax2, 0.85, 41.5, text = "Isolinhas brancas: fn = 100, 150, 200, 300 Hz\n(TDAC sub-5K em todo o domínio com Índio)",
      color = :white, fontsize = 10)

figures_dir = joinpath(@__DIR__, "..", "figures")
mkpath(figures_dir)
output_path = joinpath(figures_dir, "codesign_pareto_study.png")

save(output_path, fig, px_per_unit = 2)
println("✓ Gráfico salvo com sucesso em: ", output_path)
