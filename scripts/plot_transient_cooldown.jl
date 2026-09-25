"""
Simulação Transiente de Resfriamento da Nanoestação Criogênica EMA (300 K -> 4.2 K)
==================================================================================

Este script executa a integração transiente não-linear da rede térmica de 5 nós da estação EMA
utilizando o solver L-estável de Rosenbrock com calor específico c_p(T) variável do NIST.

Painéis da Figura Gerada:
1. Painel A (Visão Global 0 a 3 Horas): Curvas de temperatura T(t) de todos os componentes
   conforme o criocooler resfria de 300 K até a base criogênica.
2. Painel B (Regime Criogênico Final e Gradientes): Zoom nos últimos 60 minutos (T < 50 K),
   destacando a temperatura de estabilização da Mini-DAC sob carga térmica de feixe síncrotron (10 mW).
3. Painel C (Física Criogênica: Por que LTI Falha?): Variação de mais de 3 ordens de grandeza
   no calor específico c_p(T) dos materiais (Cobre, CuBe, Inox 304) entre 300 K e 4.2 K.
"""

using CairoMakie
using LinearAlgebra

using CryoThermal

println("Configurando modelo térmico de 5 nós da estação EMA...")

# Materiais
cu   = CopperOFHC()
ss   = StainlessSteel304()
cube = BerylliumCopper()

# 1. Definição dos Nós Térmicos
node_cold   = ThermalNode("Cabeçote Frio Criocooler", Inf, cu; is_fixed=true, fixed_temp=4.2)
node_braid  = ThermalNode("Cordoalha de Cobre", 0.020, cu)
node_dac    = ThermalNode("Mini-DAC (CuBe)", 0.051, cube; heat_load=0.010) # 10 mW feixe
node_supp   = ThermalNode("Suporte Inox 304", 0.015, ss)
node_shield = ThermalNode("Escudo Térmico 40 K", Inf, cu; is_fixed=true, fixed_temp=40.0)

# 2. Definição dos Elos de Ligação
# Interface sapata/cabeçote com elo dinâmico de Índio (resistência de contato variável com a temperatura)
joint_in = IndiumBoltedJoint(10e-4; num_bolts=2, bolt_diameter=3e-3, torque=0.5, has_indium=true)
l_cont   = IndiumContactLink("Sapata Dedo Frio - Cordoalha (Índio)", 1, 2, joint_in)
l_braid   = ConductionLink("Cordoalha Cu (50 mm² x 50 mm)", 2, 3, 50e-6, 0.050, cu)
l_rad_dac = RadiationLink("Radiação Escudo -> DAC", 5, 3, 0.0020, 0.10)
l_supp    = ConductionLink("Hastes Suporte Inox (DAC -> Meio)", 3, 4, 3.5e-6, 0.030, ss)
l_shld    = ConductionLink("Hastes Suporte Inox (Meio -> Escudo)", 4, 5, 3.5e-6, 0.010, ss)

sys_ema = ThermalSystem([node_cold, node_braid, node_dac, node_supp, node_shield],
                        [l_cont, l_braid, l_rad_dac, l_supp, l_shld])

# 3. Dinâmica das Fronteiras (Resfriamento Exponencial do Criocooler)
tau_cooler = 1800.0 # constante de tempo do criocooler = 30 min (típico de GM / Pulse Tube)
fn_cold(t)   = exponential_cryocooler_cooldown(t; T_start=300.0, T_final=4.2, tau=tau_cooler)
fn_shield(t) = exponential_cryocooler_cooldown(t; T_start=300.0, T_final=40.0, tau=tau_cooler)

boundary_map = Dict{Int, Function}(1 => fn_cold, 5 => fn_shield)

t_final_sim = 3.0 * 3600.0 # 3 horas
println("Iniciando integração temporal transiente não-linear (0 a 3 horas)...")
@time res_trans = solve_transient(sys_ema, (0.0, t_final_sim);
                                  T_init=300.0,
                                  boundary_fns=boundary_map,
                                  dt_init=2.0,
                                  tol=1e-2)

println("Simulação concluída com sucesso! Passos adaptativos: ", length(res_trans.times))
t_hours = res_trans.times ./ 3600.0
T_matrix = res_trans.temperatures

# Temperaturas finais aos 3 horas
T_cold_end   = T_matrix[1, end]
T_braid_end  = T_matrix[2, end]
T_dac_end    = T_matrix[3, end]
T_supp_end   = T_matrix[4, end]
T_shield_end = T_matrix[5, end]

println("\n=== Temperaturas Finais ao Fim do Cooldown (t = 3.0 h) ===")
println("Cabeçote Frio: ", round(T_cold_end, digits=3), " K")
println("Cordoalha Cu:  ", round(T_braid_end, digits=3), " K")
println("Mini-DAC:      ", round(T_dac_end, digits=3), " K")
println("Suporte SS304: ", round(T_supp_end, digits=3), " K")
println("Escudo 40K:    ", round(T_shield_end, digits=3), " K")
println("Gradiente Térmico (DAC - ColdHead): ", round(T_dac_end - T_cold_end, digits=3), " K")

# =============================================================================
# GERAÇÃO DA FIGURA VETORIAL COM CAIROMAKIE
# =============================================================================
fig = Figure(size = (1500, 900), fontsize = 14)

# Paleta de cores consistente
c_cold   = RGBf(0.0, 0.45, 0.75)  # Azul criocooler
c_braid  = RGBf(0.85, 0.45, 0.1)  # Cobre metálico
c_dac    = RGBf(0.80, 0.0, 0.15)  # Vermelho DAC
c_supp   = RGBf(0.40, 0.40, 0.45) # Cinza aço inox
c_shield = RGBf(0.15, 0.65, 0.35) # Verde escudo

# -----------------------------------------------------------------------------
# PAINEL A: RESFRIAMENTO GLOBAL (0 a 3 Horas)
# -----------------------------------------------------------------------------
ax1 = Axis(fig[1, 1:2],
    title = "A) Resfriamento Transiente Global da Nanoestação EMA (300 K → 4.2 K)",
    xlabel = "Tempo de Resfriamento [horas]",
    ylabel = "Temperatura [K]",
    titlesize = 16,
    xgridvisible = true,
    ygridvisible = true
)

lines!(ax1, t_hours, T_matrix[1, :], color = c_cold,   linewidth = 3, linestyle = :dash, label = "Dedo Frio Criocooler (Fonte 4.2 K)")
lines!(ax1, t_hours, T_matrix[5, :], color = c_shield, linewidth = 2.5, linestyle = :dot,  label = "Escudo de Radiação (Fonte 40 K)")
lines!(ax1, t_hours, T_matrix[2, :], color = c_braid,  linewidth = 2.5, label = "Cordoalha Cu OFHC (Sapata c/ Índio)")
lines!(ax1, t_hours, T_matrix[3, :], color = c_dac,    linewidth = 3.5, label = "Mini-DAC CuBe (Amostra + Feixe 10 mW)")
lines!(ax1, t_hours, T_matrix[4, :], color = c_supp,   linewidth = 2,   label = "Suporte Estrutural Inox 304")

# Linhas de referência física
hlines!(ax1, [77.0], color = :gray60, linestyle = :dashdot, linewidth = 1)
text!(ax1, 0.1, 82.0, text = "Patamar N₂ Líquido (77 K)", color = :gray40, fontsize = 12)

hlines!(ax1, [4.2], color = :blue, linestyle = :dashdot, linewidth = 1)
text!(ax1, 0.1, 9.0, text = "Patamar He Líquido (4.2 K)", color = :blue, fontsize = 12)

axislegend(ax1, position = :rt, framevisible = true, backgroundcolor = (:white, 0.9))

# -----------------------------------------------------------------------------
# PAINEL B: POUSO CRIOGÊNICO E REGIME ESTACIONÁRIO (Últimos 90 minutos)
# -----------------------------------------------------------------------------
ax2 = Axis(fig[2, 1],
    title = "B) Pouso Criogênico e Gradiente Térmico Final (t ≥ 1.0 h)",
    xlabel = "Tempo [horas]",
    ylabel = "Temperatura [K]",
    titlesize = 15,
    xgridvisible = true,
    ygridvisible = true
)

xlims!(ax2, 1.0, 3.0)
ylims!(ax2, 3.0, 50.0)

lines!(ax2, t_hours, T_matrix[1, :], color = c_cold,   linewidth = 2.5, linestyle = :dash, label = "Dedo Frio (4.2 K)")
lines!(ax2, t_hours, T_matrix[5, :], color = c_shield, linewidth = 2.5, linestyle = :dot,  label = "Escudo Térmico (40 K)")
lines!(ax2, t_hours, T_matrix[2, :], color = c_braid,  linewidth = 2.5, label = "Cordoalha Cu")
lines!(ax2, t_hours, T_matrix[3, :], color = c_dac,    linewidth = 3.5, label = "Mini-DAC (CuBe)")
lines!(ax2, t_hours, T_matrix[4, :], color = c_supp,   linewidth = 2,   label = "Suporte Inox")

# Anotação do gradiente de temperatura
delta_T_final = T_dac_end - T_cold_end
text!(ax2, 1.8, 12.0,
      text = "Estabilização da Mini-DAC:\nT_DAC = $(round(T_dac_end, digits=2)) K\nΔT_cold = +$(round(delta_T_final, digits=2)) K\n(Sob Feixe 10 mW + Sapata c/ Índio)",
      color = c_dac, fontsize = 12)

axislegend(ax2, position = :rt, framevisible = true, backgroundcolor = (:white, 0.85))

# -----------------------------------------------------------------------------
# PAINEL C: POR QUE O MODELO LTI FALHA? VARIAÇÃO DE cp(T)
# -----------------------------------------------------------------------------
ax3 = Axis(fig[2, 2],
    title = "C) Por que Modelos Lineares (LTI) Falham: cp(T) do NIST",
    xlabel = "Temperatura [K]",
    ylabel = "Calor Específico cp [J/(kg·K)]",
    yscale = log10,
    xscale = log10,
    titlesize = 15,
    xgridvisible = true,
    ygridvisible = true
)

T_range = 10.0 .^ range(log10(4.0), log10(300.0), length=150)
cp_cu   = [CryoThermal.specific_heat(cu, t) for t in T_range]
cp_cube = [CryoThermal.specific_heat(cube, t) for t in T_range]
cp_ss   = [CryoThermal.specific_heat(ss, t) for t in T_range]

lines!(ax3, T_range, cp_cu,   color = c_braid, linewidth = 3, label = "Cobre OFHC (Cordoalha)")
lines!(ax3, T_range, cp_cube, color = c_dac,   linewidth = 3, label = "CuBe Alloy 25 (Mini-DAC)")
lines!(ax3, T_range, cp_ss,   color = c_supp,  linewidth = 3, label = "Aço Inox AISI 304 (Suporte)")

# Anotações destacando a queda de 3 a 4 ordens de grandeza
text!(ax3, 5.0, 0.05,
      text = "A 4.2 K:\ncp(Cu) ≈ 0.11 J/(kg·K)\n(Queda de >3500x vs 300 K!)",
      color = c_braid, fontsize = 12)

text!(ax3, 30.0, 180.0,
      text = "A 300 K:\ncp ≈ 400-500 J/(kg·K)\n(Inércia térmica clássica)",
      color = :gray20, fontsize = 12)

axislegend(ax3, position = :rb, framevisible = true, backgroundcolor = (:white, 0.85))

# Salva a figura
output_path = joinpath(@__DIR__, "..", "figures", "transient_cooldown_study.png")
save(output_path, fig, px_per_unit = 2)
println("Figura salva com sucesso em: ", output_path)
