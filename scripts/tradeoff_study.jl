"""
Estudo de Compromisso Termo-Mecânico (Co-Design) da Nanoestação EMA
==================================================================

OBJETIVO:
Mapear o trade-off entre a estabilidade mecânica (frequência natural fn) e o desempenho
térmico criogênico (temperatura final da Mini-DAC TDAC), variando:
1. A geometria do suporte mecânico de Inox 304 (área transversal A e rigidez axial k_axial);
2. A qualidade da interface térmica da cordoalha (resistência de contato R_contato,
   representando montagem com folha de Índio vs contato seco).

REQUISITOS DA LINHA EMA (SÍRIUS/LNLS):
- Mecânico: fn >= 150 Hz (evitar ressonância com vibrações do criocooler e piso do sincrotron)
- Térmico:  TDAC <= 6.0 K (manter a amostra em regime criogênico sob feixe de raios X de 10 mW)
"""

using CairoMakie
using LinearAlgebra

# Carrega o framework CryoThermal
include(joinpath(@__DIR__, "..", "src", "CryoThermal.jl"))
using .CryoThermal

println("Iniciando varredura paramétrica termo-mecânica...")

# Materiais
cu   = CopperOFHC()
ss   = StainlessSteel304()
cube = BerylliumCopper()

# Parâmetros fixos da nanoestação
m_dac  = 0.0513 # Massa da Mini-DAC = 51.3 g
q_beam = 10e-3  # Aporte térmico do feixe síncrotron = 10 mW
L_sup  = 25e-3  # Comprimento dos suportes = 25 mm
E_ss   = youngs_modulus(ss, 20.0) # ~210 GPa a frio

# Faixa de frequências naturais desejadas (de 60 Hz a 350 Hz)
# Como fn = (1 / 2π) * √(k_axial / m_dac) e k_axial = E * A / L,
# podemos calcular a área transversal necessária para cada frequência:
fn_values = range(60.0, 350.0, length=40)
areas_sup = [(2π * fn)^2 * m_dac * L_sup / E_ss for fn in fn_values]

# Casos de interface térmica da cordoalha (R_contato com o dedo frio a 4.2 K)
# 1. R_contato = 0.1 K/W -> Excelente contato com folha de Índio bem prensada (>10 MPa)
# 2. R_contato = 0.5 K/W -> Contato típico com folha de Índio (~5 MPa)
# 3. R_contato = 2.0 K/W -> Contato imperfeito ou degradado
# 4. R_contato = 6.0 K/W -> Contato seco sem Índio
contact_cases = [
    (R = 0.1, label = "R_contato = 0.1 K/W (Índio Ótimo)", color = :blue),
    (R = 0.5, label = "R_contato = 0.5 K/W (Índio Típico)", color = :teal),
    (R = 2.0, label = "R_contato = 2.0 K/W (Índio Degradado)", color = :orange),
    (R = 6.0, label = "R_contato = 6.0 K/W (Contato Seco)", color = :red)
]

# Matriz para armazenar as curvas TDAC
T_dac_results = zeros(length(contact_cases), length(fn_values))

for (c_idx, case) in enumerate(contact_cases)
    println("Calculando caso: $(case.label)...")
    for (f_idx, A_sup) in enumerate(areas_sup)
        
        # Montagem da rede térmica com a geometria atual do suporte
        n1 = ThermalNode("Cryostat Cold Head", Inf, cu; is_fixed=true, fixed_temp=4.2)
        n2 = ThermalNode("Cu Braid", 0.134, cu; heat_load=0.0)
        n3 = ThermalNode("Mini-DAC (CuBe)", m_dac, cube; heat_load=q_beam)
        n4 = ThermalNode("SS304 Support", 0.332, ss; heat_load=0.0)
        n5 = ThermalNode("40K Radiation Shield", Inf, ss; is_fixed=true, fixed_temp=40.0)
        
        # Elos:
        l1 = ContactLink("ColdHead - Braid Contact", 1, 2, case.R)
        l2 = ConductionLink("Copper Braid Conduction", 2, 3, 12e-6, 60e-3, cu) # Braid de 12 mm^2
        l3 = ConductionLink("DAC Support Strut", 3, 4, A_sup, L_sup, ss)       # Geometria variável
        l4 = ConductionLink("Support to Shield", 4, 5, 25e-6, 20e-3, ss)
        l5 = RadiationLink("Shield Radiation to DAC", 5, 3, 10e-4, 0.05)
        
        sys = ThermalSystem([n1, n2, n3, n4, n5], [l1, l2, l3, l4, l5])
        
        res = solve_steady_state(sys; tol=1e-6)
        T_dac_results[c_idx, f_idx] = res.temperatures[3]
    end
end

println("Simulações concluídas com sucesso. Gerando figura com CairoMakie...")

# =============================================================================
# GERAÇÃO DA FIGURA COM CAIROMAKIE
# =============================================================================
fig = Figure(size = (1000, 650), fontsize = 14, fonts = (; regular = "sans-serif", bold = "sans-serif"))

# Painel Principal: Curvas de Trade-off
ax = Axis(fig[1, 1],
    title = "Compromisso Termo-Mecânico: Rigidez Estrutural vs Temperatura Criogênica (Linha EMA)",
    xlabel = "Primeira Frequência Natural da Montagem fn [Hz]  (Rigidez)",
    ylabel = "Temperatura de Equilíbrio da Mini-DAC T_DAC [K]",
    limits = ((50, 360), (4.0, 12.0))
)

# 1. Região de Requisitos (Envelope de Projeto da Linha EMA)
# Requisito mecânico: fn >= 150 Hz
vspan!(ax, 50.0, 150.0, color = (:red, 0.08), label = "Região Mecanicamente Inviável (fn < 150 Hz)")
vlines!(ax, [150.0], color = :red, linestyle = :dash, linewidth = 2)
text!(ax, 153.0, 11.5, text = "Corte Mecânico: fn ≥ 150 Hz", color = :red, fontsize = 12)

# Requisito térmico: TDAC <= 6.0 K
hlines!(ax, [6.0], color = :darkred, linestyle = :dot, linewidth = 2)
text!(ax, 55.0, 6.2, text = "Meta Criogênica: TDAC ≤ 6.0 K", color = :darkred, fontsize = 12)

# Destacar a Zona de Projeto Viável
poly!(ax, Point2f[(150, 4.0), (360, 4.0), (360, 6.0), (150, 6.0)],
      color = (:green, 0.12), label = "Envelope de Projeto Viável")
text!(ax, 200.0, 4.5, text = "ZONA ÓTIMA DE PROJETO\n(fn ≥ 150 Hz & TDAC ≤ 6 K)", color = :darkgreen, fontsize = 13, font = :bold)

# 2. Plotar as curvas para cada condição de contato
for (c_idx, case) in enumerate(contact_cases)
    lines!(ax, fn_values, T_dac_results[c_idx, :],
           label = case.label, color = case.color, linewidth = 2.5)
    scatter!(ax, fn_values[1:5:end], T_dac_results[c_idx, 1:5:end],
             color = case.color, markersize = 8)
end

# Legenda bem posicionada
axislegend(ax, position = :lt, framevisible = true, bgcolor = (:white, 0.85))

# Garante que a pasta 'figures' existe
figures_dir = joinpath(@__DIR__, "..", "figures")
mkpath(figures_dir)
output_path = joinpath(figures_dir, "tradeoff_ema_nanostation.png")

save(output_path, fig, px_per_unit = 2)
println("✓ Gráfico salvo com sucesso em: ", output_path)
