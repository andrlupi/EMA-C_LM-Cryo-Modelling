"""
Suite de Testes e Validação Física do Modelo Criogênico da Linha EMA
====================================================================

OBJETIVOS DESTA SUÍTE:
1. Validar a termofísica do NIST: k, cp, E, ΔL/L para Cu, SS304, CuBe, Be e Ti-6Al-4V.
2. Validar limites físicos sub-Kelvin (Terceira Lei: k -> 0, cp -> 0) e guarda para T > 350 K.
3. Testar a ausência de alocações na heap (zero allocations) nas rotinas do inner-loop:
   thermal_conductivity, specific_heat e thermal_conductivity_integral.
4. Validar o cálculo da Transformada de Kirchhoff via quadratura de Gauss-Legendre de 7 pontos.
5. Testar elos dinâmicos IndiumContactLink e ExchangeGasLink acoplados a ThermalSystem.
6. Testar o solver de Newton-Raphson com Jacobiano exato ForwardDiff e busca linear retrógrada.
7. Simular o modelo real de 5 massas da nanoestação EMA em regime estacionário e transiente.
"""

using Test
using LinearAlgebra
using ForwardDiff
using CryoThermal

@testset "CryoThermal Framework Tests" begin
    
    # =========================================================================
    # TESTE 1: PROPRIEDADES DOS MATERIAIS CRIOGÊNICOS DO NIST
    # =========================================================================
    @testset "1. NIST Material Properties & New Materials" begin
        cu   = CopperOFHC()
        ss   = StainlessSteel304()
        cube = BerylliumCopper()
        be   = Beryllium()
        ti   = TitaniumTi6Al4V()
        
        # Teste de sanidade: as propriedades devem ser estritamente positivas em toda a faixa criogênica
        for T in [4.2, 20.0, 77.0, 300.0]
            k_cu   = thermal_conductivity(cu, T)
            k_ss   = thermal_conductivity(ss, T)
            k_cube = thermal_conductivity(cube, T)
            k_be   = thermal_conductivity(be, T)
            k_ti   = thermal_conductivity(ti, T)
            
            @test k_cu > 0
            @test k_ss > 0
            @test k_cube > 0
            @test k_be > 0
            @test k_ti > 0
            
            # FÍSICA: Cu >> SS304 e Ti-6Al-4V é ainda mais isolante que o Inox a 4.2 K
            @test k_cu > k_ss
            if T == 4.2
                @test k_ti < k_ss # Ti-6Al-4V tem condutividade térmica menor que SS304 a 4.2 K (~0.12 vs ~0.29 W/m·K)
            end
            
            cp_cu   = specific_heat(cu, T)
            cp_ss   = specific_heat(ss, T)
            cp_cube = specific_heat(cube, T)
            cp_be   = specific_heat(be, T)
            cp_ti   = specific_heat(ti, T)
            
            @test cp_cu > 0
            @test cp_ss > 0
            @test cp_cube > 0
            @test cp_be > 0
            @test cp_ti > 0
            
            # Densidades positivas
            @test density(cu) > 0
            @test density(ss) > 0
            @test density(cube) > 0
            @test density(be) > 0
            @test density(ti) > 0
            
            # Módulos de elasticidade
            @test youngs_modulus(cu, T) > 0
            @test youngs_modulus(ss, T) > 0
            @test youngs_modulus(cube, T) > 0
            @test youngs_modulus(be, T) > 0
            @test youngs_modulus(ti, T) > 0
        end
        
        # TESTE DA TRANSFORMADA DE KIRCHHOFF:
        theta_cu = thermal_conductivity_integral(cu, 4.2, 40.0)
        theta_ss = thermal_conductivity_integral(ss, 4.2, 40.0)
        theta_ti = thermal_conductivity_integral(ti, 4.2, 40.0)
        
        @test theta_cu > 0
        @test theta_ss > 0
        @test theta_ti > 0
        @test theta_cu > 100 * theta_ss
        @test theta_ti < theta_ss # Ti conduz menos calor integrado que SS304 entre 4.2 K e 40 K
        
        # Teste de caso trivial T1 == T2
        @test thermal_conductivity_integral(cu, 10.0, 10.0) == 0.0
    end
    
    # =========================================================================
    # TESTE 2: CONTRAÇÃO TÉRMICA INTEGRADA (THERMAL EXPANSION)
    # =========================================================================
    @testset "2. Thermal Expansion (ΔL/L_293)" begin
        cu   = CopperOFHC()
        ss   = StainlessSteel304()
        cube = BerylliumCopper()
        be   = Beryllium()
        ti   = TitaniumTi6Al4V()
        
        for mat in [cu, ss, cube, be, ti]
            # Em T = 293 K, ΔL/L_293 deve ser nulo (ou ordens de magnitude menor que 1e-4)
            le_293 = thermal_expansion(mat, 293.0)
            @test abs(le_293) < 5e-4
            
            # A 4.2 K, os materiais sofrem contração térmica líquida (ΔL/L < 0)
            le_4k = thermal_expansion(mat, 4.2)
            @test le_4k < 0.0
            
            # A 77 K, a contração permanece estritamente negativa em relação ao ambiente (293 K)
            @test thermal_expansion(mat, 77.0) < thermal_expansion(mat, 293.0)
        end
        
        # Monotonicidade estrita entre 4.2 K e 77 K para materiais estruturais/condutores (Cu, SS304, CuBe, Ti-6Al-4V)
        # (O polinômio NIST do Berílio policristalino possui um leve platô empírico sub-70 K de Δ(ΔL/L) ~ 0.01%)
        for mat in [cu, ss, cube, ti]
            @test thermal_expansion(mat, 4.2) < thermal_expansion(mat, 77.0)
        end
    end
    
    # =========================================================================
    # TESTE 3: LIMITES FÍSICOS SUB-KELVIN E ALTA TEMPERATURA (> 350 K)
    # =========================================================================
    @testset "3. Sub-Kelvin & High-Temperature Physical Safeguards" begin
        cu   = CopperOFHC()
        ss   = StainlessSteel304()
        cube = BerylliumCopper()
        be   = Beryllium()
        ti   = TitaniumTi6Al4V()
        
        # 3.1 Limite Sub-Kelvin e Terceira Lei da Termodinâmica para todos os 5 materiais
        for mat in [cu, ss, cube, be, ti]
            @test thermal_conductivity(mat, 0.0) == 0.0
            @test specific_heat(mat, 0.0) == 0.0
            
            # Em 0.5 K, k e cp devem ser estritamente positivos e menores que a 1.0 K
            @test 0.0 < thermal_conductivity(mat, 0.5) < thermal_conductivity(mat, 1.0)
            @test 0.0 < specific_heat(mat, 0.5) < specific_heat(mat, 1.0)
            
            # Derivada estritamente positiva em T = 0 K (sensibilidade garantida para Newton-Raphson/ForwardDiff)
            @test ForwardDiff.derivative(t -> thermal_conductivity(mat, t), 0.0) > 0.0
            @test ForwardDiff.derivative(t -> specific_heat(mat, t), 0.0) > 0.0
            @test ForwardDiff.derivative(t -> thermal_conductivity(mat, t), 0.5) > 0.0
        end
        
        # Derivada estritamente positiva em T = 0 K para interfaces frias
        joint_test = IndiumBoltedJoint(10e-4)
        gap_test   = HeliumExchangeGasGap(2e-3, 10e-4)
        @test ForwardDiff.derivative(t -> contact_conductance_indium(joint_test, t), 0.0) > 0.0
        @test ForwardDiff.derivative(t -> gas_gap_conductance(gap_test, t), 0.0) > 0.0
        
        # 3.2 Teto de alta temperatura (T > 350 K) não deve divergir para infinito/NaN
        for mat in [cu, ss, cube, be, ti]
            k_high = thermal_conductivity(mat, 400.0)
            @test isfinite(k_high)
            @test k_high > 0.0
            
            cp_high = specific_heat(mat, 400.0)
            @test isfinite(cp_high)
            @test cp_high > 0.0

            # Diferenciação automática de segunda ordem (AD aninhada / números Dual compostos)
            d2_k = ForwardDiff.derivative(t -> ForwardDiff.derivative(x -> thermal_conductivity(mat, x), t), 400.0)
            d2_cp = ForwardDiff.derivative(t -> ForwardDiff.derivative(x -> specific_heat(mat, x), t), 400.0)
            @test isfinite(d2_k)
            @test isfinite(d2_cp)
            @test d2_k == 0.0
            @test d2_cp == 0.0
        end
    end
    
    # =========================================================================
    # TESTE 4: VALIDAÇÃO DE ZERO ALOCAÇÕES NA HEAP (INNER LOOP)
    # =========================================================================
    @testset "4. Zero-Allocation Benchmarks (Heap Invariance)" begin
        cu   = CopperOFHC()
        ss   = StainlessSteel304()
        cube = BerylliumCopper()
        be   = Beryllium()
        ti   = TitaniumTi6Al4V()
        
        all_materials = [cu, ss, cube, be, ti]
        
        # Aquecimento do compilador JIT (warmup) para todos os materiais
        for mat in all_materials
            thermal_conductivity(mat, 10.0)
            specific_heat(mat, 10.0)
            youngs_modulus(mat, 20.0)
            thermal_expansion(mat, 20.0)
            thermal_conductivity_integral(mat, 4.2, 40.0)
        end
        
        # Asserções de Alocação Exata = 0 Bytes para todos os 5 materiais
        for mat in all_materials
            @test (@allocated thermal_conductivity(mat, 10.0)) == 0
            @test (@allocated specific_heat(mat, 10.0)) == 0
            @test (@allocated youngs_modulus(mat, 20.0)) == 0
            @test (@allocated thermal_expansion(mat, 20.0)) == 0
            @test (@allocated thermal_conductivity_integral(mat, 4.2, 40.0)) == 0
        end
    end
    
    # =========================================================================
    # TESTE 5: EQUILÍBRIO CONDUÇÃO PURA EM SISTEMA DE 2 ELOS / 3 NÓS
    # =========================================================================
    @testset "5. Simple 2-Node Conduction Steady State" begin
        cu = CopperOFHC()
        ss = StainlessSteel304()
        
        node_cold = ThermalNode("Cold Head", Inf, cu; is_fixed=true, fixed_temp=4.2)
        node_mid  = ThermalNode("Sample", 0.05, cu; heat_load=0.0)
        node_hot  = ThermalNode("Shield", Inf, ss; is_fixed=true, fixed_temp=40.0)
        
        link1 = ConductionLink("Cu Braid", 1, 2, 1e-4, 0.05, cu)
        link2 = ConductionLink("SS Support", 2, 3, 0.5e-4, 0.05, ss)
        
        sys = ThermalSystem([node_cold, node_mid, node_hot], [link1, link2])
        
        result = solve_steady_state(sys)
        @test result.converged
        @test result.iterations <= 10
        
        T_mid = result.temperatures[2]
        @test 4.2 < T_mid < 40.0
        @test T_mid < 6.0
    end
    
    # =========================================================================
    # TESTE 6: ELOS DINÂMICOS INDIUMCONTACTLINK E EXCHANGEGASLINK EM THERMALSYSTEM
    # =========================================================================
    @testset "6. IndiumContactLink & ExchangeGasLink System Integration" begin
        cu = CopperOFHC()
        
        # 6.1 Teste com IndiumContactLink dinâmico acoplado
        joint = IndiumBoltedJoint(10e-4; num_bolts=2, bolt_diameter=3e-3, torque=0.8, has_indium=true)
        node1 = ThermalNode("ColdHead", Inf, cu; is_fixed=true, fixed_temp=4.2)
        node2 = ThermalNode("Block", 0.05, cu; heat_load=0.050) # 50 mW
        
        link_indium = IndiumContactLink("IndiumJoint", 1, 2, joint)
        sys_in = ThermalSystem([node1, node2], [link_indium])
        
        res_in = solve_steady_state(sys_in; tol=1e-7)
        @test res_in.converged
        @test res_in.temperatures[2] > 4.2
        # Com 50 mW e junta de Índio (condutância total ~ 1.16 W/K), o salto é minúsculo (< 0.1 K)
        @test res_in.temperatures[2] - 4.2 < 0.1
        # Fluxo de calor balanceado
        @test isapprox(abs(res_in.heat_flows["IndiumJoint"]), 0.050, atol=1e-6)
        
        # 6.2 Teste com ExchangeGasLink dinâmico acoplado
        gap_high_p = HeliumExchangeGasGap(2e-3, 10e-4; pressure=100.0) # 1 mbar
        gap_low_p  = HeliumExchangeGasGap(2e-3, 10e-4; pressure=1.0)   # 0.01 mbar
        
        link_gas_high = ExchangeGasLink("HeGasHigh", 1, 2, gap_high_p)
        link_gas_low  = ExchangeGasLink("HeGasLow", 1, 2, gap_low_p)
        
        sys_gas_high = ThermalSystem([node1, node2], [link_gas_high])
        sys_gas_low  = ThermalSystem([node1, node2], [link_gas_low])
        
        res_gas_high = solve_steady_state(sys_gas_high; tol=1e-7)
        res_gas_low  = solve_steady_state(sys_gas_low; tol=1e-7)
        
        @test res_gas_high.converged
        @test res_gas_low.converged
        # Menor pressão = menor condutância = bloco aquece mais sob os mesmos 50 mW
        @test res_gas_low.temperatures[2] > res_gas_high.temperatures[2]
    end
    
    # =========================================================================
    # TESTE 7: REGULARIZAÇÃO DA CONDUTÂNCIA DE ÍNDIO EM TEMPERATURA AMBIENTE
    # =========================================================================
    @testset "7. Indium Contact Conductance Regularization & Dry Contact" begin
        joint_in  = IndiumBoltedJoint(10e-4; num_bolts=2, bolt_diameter=3e-3, torque=0.8, has_indium=true)
        joint_dry = IndiumBoltedJoint(10e-4; num_bolts=2, bolt_diameter=3e-3, torque=0.8, has_indium=false)
        
        hc_4k   = contact_conductance_indium(joint_in, 4.2)
        hc_20k  = contact_conductance_indium(joint_in, 20.0)
        hc_300k = contact_conductance_indium(joint_in, 300.0)
        
        # Monotônico crescente
        @test hc_4k < hc_20k < hc_300k
        
        # A 300 K a condutância NÃO deve explodir para 300.000 W/(m²·K), mas sim saturar em faixa física
        @test hc_300k < 35000.0
        @test hc_300k > 15000.0
        
        # Continuidade C1: derivada em torno de 25 K deve ser positiva e finita
        d_trans = ForwardDiff.derivative(t -> contact_conductance_indium(joint_in, t), 25.0)
        @test d_trans > 0.0
        @test isfinite(d_trans)
        
        # Contato seco a 4.2 K é ordens de grandeza inferior
        hc_dry_4k = contact_conductance_indium(joint_dry, 4.2)
        @test hc_4k > 15 * hc_dry_4k
    end
    
    # =========================================================================
    # TESTE 8: MODELO REAL DE 5 MASSAS DA NANOESTAÇÃO DA LINHA EMA
    # =========================================================================
    @testset "8. Full 5-Node EMA Nano-Station Model" begin
        cu   = CopperOFHC()
        ss   = StainlessSteel304()
        cube = BerylliumCopper()
        
        n1 = ThermalNode("Cryostat Cold Head", Inf, cu; is_fixed=true, fixed_temp=4.2)
        n2 = ThermalNode("Cu Braid", 0.134, cu; heat_load=0.0)
        n3 = ThermalNode("Mini-DAC (CuBe)", 0.051, cube; heat_load=10e-3)
        n4 = ThermalNode("SS304 Support", 0.332, ss; heat_load=0.0)
        n5 = ThermalNode("40K Radiation Shield", Inf, ss; is_fixed=true, fixed_temp=40.0)
        
        l1 = ContactLink("ColdHead - Braid Contact", 1, 2, 0.5)
        l2 = ConductionLink("Copper Braid Conduction", 2, 3, 10e-6, 60e-3, cu)
        l3 = ConductionLink("DAC Support Strut", 3, 4, 15e-6, 25e-3, ss)
        l4 = ConductionLink("Support to Shield", 4, 5, 20e-6, 20e-3, ss)
        l5 = RadiationLink("Shield Radiation to DAC", 5, 3, 10e-4, 0.05)
        
        sys = ThermalSystem([n1, n2, n3, n4, n5], [l1, l2, l3, l4, l5])
        
        # Teste de convergência a partir de chute inicial neutro
        result = solve_steady_state(sys)
        @test result.converged
        @test result.iterations <= 10
        @test result.temperatures[3] < 10.0
        
        # Teste de robustez a chute inicial adverso (T_guess quente = 200 K)
        result_hot_guess = solve_steady_state(sys; T_guess=[4.2, 200.0, 200.0, 200.0, 40.0])
        @test result_hot_guess.converged
        @test isapprox(result_hot_guess.temperatures[3], result.temperatures[3], atol=1e-5)
    end
    
    # =========================================================================
    # TESTE 9: SOLVER TRANSIENTE NÃO-LINEAR (RESFRIAMENTO COM FRONTEIRAS DINÂMICAS)
    # =========================================================================
    @testset "9. Non-Linear Transient Cooldown Solver" begin
        cu = CopperOFHC()
        
        node_cold = ThermalNode("Cold Reservoir", Inf, cu; is_fixed=true, fixed_temp=4.2)
        node_mass = ThermalNode("Copper Sample", 0.10, cu; heat_load=0.0)
        link = ConductionLink("Cu Rod", 1, 2, 1e-4, 0.05, cu)
        sys_trans = ThermalSystem([node_cold, node_mass], [link])
        
        # Simula resfriamento de 100 K até próximo de 4.2 K por 300 segundos
        res_trans = solve_transient(sys_trans, (0.0, 300.0); T_init=100.0, dt_init=0.5, tol=1e-2)
        
        @test length(res_trans.times) > 10
        T_hist = res_trans.temperatures[2, :]
        @test T_hist[end] < T_hist[1]
        @test all(diff(T_hist) .<= 1e-6)
        @test T_hist[end] < 5.0
        
        # Teste com condição de contorno temporal (resfriamento exponencial do cabeçote)
        fn_cool(t) = exponential_cryocooler_cooldown(t; T_start=300.0, T_final=4.2, tau=600.0)
        res_bmap = solve_transient(sys_trans, (0.0, 60.0); T_init=300.0, cold_head_fn=fn_cool, dt_init=1.0)
        @test res_bmap.temperatures[1, end] < res_bmap.temperatures[1, 1]
    end
    
    # =========================================================================
    # TESTE 10: TRATAMENTO DE SINGULARIDADES, CASOS DE BORDA E REDE ROBUSTA
    # =========================================================================
    @testset "10. Edge Cases, Singularity Handling & Network Robustness" begin
        cu = CopperOFHC()
        ss = StainlessSteel304()
        
        # 10.1 Detecção graciosa de Jacobiano singular (nó desconectado/isolado)
        n_fixed = ThermalNode("FixedCold", Inf, cu; is_fixed=true, fixed_temp=4.2)
        n_iso   = ThermalNode("IsolatedNode", 0.05, cu; heat_load=0.010)
        sys_singular = ThermalSystem([n_fixed, n_iso]) # Nenhum elo conectando o nó livre
        
        # Não deve lançar SingularException não-capturada; deve reportar converged = false
        res_sing = solve_steady_state(sys_singular)
        @test !res_sing.converged
        
        # 10.2 Transient solver com zero nós livres (n_free == 0)
        sys_only_fixed = ThermalSystem([n_fixed])
        res_nofree = solve_transient(sys_only_fixed, (0.0, 10.0))
        @test length(res_nofree.times) >= 2
        @test res_nofree.temperatures[1, 1] == 4.2
        @test res_nofree.temperatures[1, end] == 4.2
        
        # 10.3 Construtor flexível de ThermalSystem (links vazios e padrão)
        sys_default = ThermalSystem([n_fixed])
        sys_empty_any = ThermalSystem([n_fixed], Any[])
        @test length(sys_default.nodes) == 1
        @test length(sys_default.links) == 0
        @test length(sys_empty_any.links) == 0
        
        # 10.4 energy_balance_residuals! in-place exportada e consistente
        T_full = [4.2, 10.0]
        res_buf = zeros(Float64, 2)
        node_a = ThermalNode("A", Inf, cu; is_fixed=true, fixed_temp=4.2)
        node_b = ThermalNode("B", 0.05, cu; heat_load=0.02)
        link_ab = ConductionLink("LinkAB", 1, 2, 1e-4, 0.05, cu)
        sys_ab = ThermalSystem([node_a, node_b], [link_ab])
        
        energy_balance_residuals!(res_buf, sys_ab, T_full)
        res_alloc = energy_balance_residuals(sys_ab, T_full)
        @test res_buf == res_alloc
    end
end
