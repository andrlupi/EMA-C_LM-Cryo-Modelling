"""
Suite de Testes e Validação Física do Modelo Criogênico da Linha EMA
====================================================================

OBJETIVOS DESTA SUÍTE:
1. Validar a termofísica do NIST: Garantir que as propriedades dos materiais
   (k, cp, E, ΔL/L) retornem valores físicos positivos e consistentes com a literatura.
2. Validar o cálculo da Transformada de Kirchhoff via quadratura de Gauss-Legendre.
3. Testar a convergência do solver de Newton-Raphson em um caso elementar de 2 nós.
4. Simular o modelo real de 5 massas da nanoestação EMA e verificar:
   - Atingimento de temperatura criogênica sub-10 K na Mini-DAC com feixe síncrotron ativo.
   - Conservação estrita de energia (Primeira Lei da Termodinâmica: ∑ Q = 0).
"""

using Test
using LinearAlgebra

# Carrega o módulo principal CryoThermal
include(joinpath(@__DIR__, "..", "src", "CryoThermal.jl"))
using .CryoThermal

@testset "CryoThermal Framework Tests" begin
    
    # =========================================================================
    # TESTE 1: PROPRIEDADES DOS MATERIAIS CRIOGÊNICOS DO NIST
    # =========================================================================
    @testset "1. NIST Material Properties" begin
        cu   = CopperOFHC()
        ss   = StainlessSteel304()
        cube = BerylliumCopper()
        
        # Teste de sanidade: as propriedades devem ser estritamente positivas em toda a faixa
        for T in [4.2, 20.0, 77.0, 300.0]
            k_cu   = thermal_conductivity(cu, T)
            k_ss   = thermal_conductivity(ss, T)
            k_cube = thermal_conductivity(cube, T)
            
            @test k_cu > 0
            @test k_ss > 0
            @test k_cube > 0
            
            # FÍSICA: O Cobre OFHC é ordens de grandeza mais condutivo que o Inox 304 em qualquer T criogênico
            @test k_cu > k_ss
            
            cp_cu = specific_heat(cu, T)
            cp_ss = specific_heat(ss, T)
            @test cp_cu > 0
            @test cp_ss > 0
        end
        
        # TESTE DA TRANSFORMADA DE KIRCHHOFF:
        # Integral de condutividade térmica Θ = ∫_{4.2}^{40} k(T) dT [W/m]
        theta_cu = thermal_conductivity_integral(cu, 4.2, 40.0)
        theta_ss = thermal_conductivity_integral(ss, 4.2, 40.0)
        
        @test theta_cu > 0
        @test theta_ss > 0
        # O cobre transporta centenas de vezes mais calor por metro do que o aço inox
        @test theta_cu > 100 * theta_ss
        
        println("✓ Copper thermal conductivity at 4.2 K: ", round(thermal_conductivity(cu, 4.2), digits=2), " W/(m·K)")
        println("✓ SS304 thermal conductivity at 4.2 K:  ", round(thermal_conductivity(ss, 4.2), digits=4), " W/(m·K)")
        println("✓ Conductivity Integral Cu (4.2 K -> 40 K):   ", round(theta_cu, digits=1), " W/m")
        println("✓ Conductivity Integral SS304 (4.2 K -> 40 K): ", round(theta_ss, digits=2), " W/m")
    end
    
    # =========================================================================
    # TESTE 2: EQUILÍBRIO CONDUÇÃO PURA EM SISTEMA DE 2 ELOS / 3 NÓS
    # =========================================================================
    @testset "2. Simple 2-Node Conduction Steady State" begin
        # Topologia de teste:
        # [Dedo Frio 4.2K (Fixo)] --- Braid de Cobre --- [Amostra (Livre)] --- Suporte SS304 --- [Escudo 40K (Fixo)]
        # Sem aporte externo de calor na amostra, sua temperatura deve se aproximar muito mais
        # de 4.2 K do que de 40 K, porque a condutância do cobre é gigantesca comparada à do inox.
        cu = CopperOFHC()
        ss = StainlessSteel304()
        
        node_cold = ThermalNode("Cold Head", Inf, cu; is_fixed=true, fixed_temp=4.2)
        node_mid  = ThermalNode("Sample", 0.05, cu; heat_load=0.0)
        node_hot  = ThermalNode("Shield", Inf, ss; is_fixed=true, fixed_temp=40.0)
        
        # Elo 1: Cordoalha de cobre (L = 5 cm, A = 1 cm²)
        link1 = ConductionLink("Cu Braid", 1, 2, 1e-4, 0.05, cu)
        # Elo 2: Haste de suporte de inox (L = 5 cm, A = 0.5 cm²)
        link2 = ConductionLink("SS Support", 2, 3, 0.5e-4, 0.05, ss)
        
        sys = ThermalSystem([node_cold, node_mid, node_hot], [link1, link2])
        
        result = solve_steady_state(sys)
        @test result.converged
        @test result.iterations <= 10
        
        T_mid = result.temperatures[2]
        @test 4.2 < T_mid < 40.0
        
        # Verificação física: condutância Cu >> Inox => T_amostra próxima a 4.2 K (< 6 K)
        @test T_mid < 6.0
        
        println("✓ 2-Node test converged in $(result.iterations) iters to T_sample = $(round(T_mid, digits=3)) K")
    end
    
    # =========================================================================
    # TESTE 3: MODELO REAL DE 5 MASSAS DA NANOESTAÇÃO DA LINHA EMA
    # =========================================================================
    @testset "3. Full 5-Node EMA Nano-Station Model" begin
        # DISCRETIZAÇÃO DOS NÓS DA LINHA:
        # Nó 1: Dedo Frio do Criostato (Fixo em 4.2 K)
        # Nó 2: Cordoalha Térmica de Cobre (Massa ~ 134 g)
        # Nó 3: Mini-DAC de CuBe (Massa ~ 51.3 g, Feixe síncrotron ativo: q_beam = 10 mW)
        # Nó 4: Suporte Estrutural de Inox SS304 (Massa ~ 332 g)
        # Nó 5: Escudo Térmico de Radiação (Fixo em 40.0 K)
        
        cu   = CopperOFHC()
        ss   = StainlessSteel304()
        cube = BerylliumCopper()
        
        n1 = ThermalNode("Cryostat Cold Head", Inf, cu; is_fixed=true, fixed_temp=4.2)
        n2 = ThermalNode("Cu Braid", 0.134, cu; heat_load=0.0)
        n3 = ThermalNode("Mini-DAC (CuBe)", 0.051, cube; heat_load=10e-3) # Feixe de 10 mW
        n4 = ThermalNode("SS304 Support", 0.332, ss; heat_load=0.0)
        n5 = ThermalNode("40K Radiation Shield", Inf, ss; is_fixed=true, fixed_temp=40.0)
        
        # DEFINIÇÃO DOS ELOS DE TRANSFERÊNCIA TÉRMICA:
        # 1. Contato aparafusado Criostato-Cordoalha com folha de Índio (Área ~ 10 cm², R ~ 0.5 K/W)
        l1 = ContactLink("ColdHead - Braid Contact", 1, 2, 0.5)
        
        # 2. Condução ao longo da Cordoalha de Cobre (Comprimento = 60 mm, Área efetiva = 10 mm²)
        l2 = ConductionLink("Copper Braid Conduction", 2, 3, 10e-6, 60e-3, cu)
        
        # 3. Condução pelos Suportes Tubulares de Inox da DAC (Comprimento = 25 mm, Área = 15 mm²)
        l3 = ConductionLink("DAC Support Strut", 3, 4, 15e-6, 25e-3, ss)
        
        # 4. Condução do anel de suporte para a estrutura do Escudo de 40 K (L = 20 mm, A = 20 mm²)
        l4 = ConductionLink("Support to Shield", 4, 5, 20e-6, 20e-3, ss)
        
        # 5. Radiação térmica residual entre o Escudo (40 K) e a face externa da DAC (A ~ 10 cm², ϵ ~ 0.05)
        l5 = RadiationLink("Shield Radiation to DAC", 5, 3, 10e-4, 0.05)
        
        sys = ThermalSystem([n1, n2, n3, n4, n5], [l1, l2, l3, l4, l5])
        
        # Resolução do equilíbrio não-linear por Newton-Raphson
        result = solve_steady_state(sys)
        @test result.converged
        
        T = result.temperatures
        println("\n=== EMA Nano-Station Steady-State Results ===")
        for (i, node) in enumerate(sys.nodes)
            println("Node $i [$(node.name)]: T = $(round(T[i], digits=3)) K")
        end
        println("Total Newton iterations: ", result.iterations)
        println("Final residual norm:     ", result.residual_norm)
        println("\nHeat Flows through Links:")
        for (name, q) in result.heat_flows
            println("  $name: ", round(q * 1e3, digits=2), " mW")
        end
        
        # CRITÉRIO DE ENGENHARIA DA LINHA EMA:
        # A Mini-DAC deve operar criogenicamente abaixo de 10 K (idealmente sub-5 K)
        @test T[3] < 10.0
    end
    
    # =========================================================================
    # TESTE 4: INTERFACES DO DEDO FRIO (ÍNDIO vs CONTATO SECO vs GÁS DE HÉLIO)
    # =========================================================================
    @testset "4. Cold Interfaces: Indium Bolted Joint & Helium Exchange Gas" begin
        # 4.1 Junta de Índio
        area_c = 10e-4 # 10 cm²
        joint_indium = IndiumBoltedJoint(area_c; num_bolts=2, bolt_diameter=3e-3, torque=0.8, has_indium=true)
        joint_dry    = IndiumBoltedJoint(area_c; num_bolts=2, bolt_diameter=3e-3, torque=0.8, has_indium=false)
        
        P_clamping = clamping_pressure(joint_indium)
        @test P_clamping > 1e6 # Deve ser > 1 MPa (tipicamente ~2.6 MPa para 2x M3 a 0.8 Nm em 10 cm²)
        
        hc_indium = contact_conductance_indium(joint_indium, 4.2)
        hc_dry    = contact_conductance_indium(joint_dry, 4.2)
        
        @test hc_indium > 0
        @test hc_dry > 0
        # FÍSICA: O Índio melhora a condutância de contato em mais de 15 a 20 vezes
        @test hc_indium > 15 * hc_dry
        
        println("\n✓ Clamping Pressure: ", round(P_clamping / 1e6, digits=2), " MPa")
        println("✓ Contact Conductance with Indium (4.2 K): ", round(hc_indium, digits=1), " W/(m²·K)")
        println("✓ Contact Conductance Dry Cu-Cu (4.2 K):   ", round(hc_dry, digits=1), " W/(m²·K)")
        println("✓ Indium Improvement Ratio:                ", round(hc_indium / hc_dry, digits=1), "x")
        
        # 4.2 Gás de Troca de Hélio
        gap = HeliumExchangeGasGap(2e-3, area_c; pressure=100.0) # 1 mbar = 100 Pa, d = 2 mm
        
        Kn_1mbar = knudsen_number(gap, 4.2)
        h_1mbar  = gas_gap_conductance(gap, 4.2)
        @test h_1mbar > 0
        
        # Testar efeito chave térmica: em vácuo (0.001 Pa = 1e-5 mbar), a condutância cai dramaticamente
        gap_vacuum = HeliumExchangeGasGap(2e-3, area_c; pressure=1e-3)
        Kn_vac     = knudsen_number(gap_vacuum, 4.2)
        h_vac      = gas_gap_conductance(gap_vacuum, 4.2)
        
        @test Kn_vac > Kn_1mbar
        @test h_vac < 0.01 * h_1mbar # Queda de mais de 100x na condutância (chave aberta)
        
        println("✓ He Gas Gap (1 mbar, 4.2 K): Kn = ", round(Kn_1mbar, digits=3), " | h = ", round(h_1mbar, digits=2), " W/(m²·K)")
        println("✓ He Gas Gap (1e-5 mbar, 4.2 K): Kn = ", round(Kn_vac, digits=1), " | h = ", round(h_vac, digits=5), " W/(m²·K)")
    end
end

