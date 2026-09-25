"""
Módulo NIST_Data
================

Carregamento e tratamento dos dados termofísicos do NIST (National Institute of Standards
and Technology - Cryogenic Material Properties Database).

POR QUE ESTE MÓDULO EXISTE?
Em física de baixas temperaturas, as propriedades térmicas (condutividade k, calor específico cp)
e mecânicas (módulo de Young E, expansão linear ΔL/L) sofrem variações de ordens de grandeza
entre 4 K e 300 K. O NIST compilou ajustes polinomiais empíricos amplamente aceitos como padrão
áureo na criogenia internacional.

ESTRUTURA DAS TABELAS CSV:
- Linhas: Coeficientes polinomiais a, b, c, d, e, f, g, h, i.
- Colunas:
  * TC / TC1..TC5 : Thermal Conductivity (Condutividade Térmica)
  * SH            : Specific Heat (Calor Específico)
  * YM / YM1..YM2 : Young's Modulus (Módulo de Elasticidade / Rigidez)
  * LE            : Linear Expansion (Contração Térmica ΔL/L relativo a 293 K)
  * EC            : Electrical Conductivity (Condutividade Elétrica)

Fonte oficial: https://trc.nist.gov/cryogenics/materials/materialproperties.htm
"""
module NIST_Data

using DelimitedFiles

export NISTMaterialData, load_nist_data

"""
    NISTMaterialData

Estrutura imutável que encapsula as matrizes numéricas e tuplas estáticas de coeficientes NIST para:
- `cu_c`: Cobre eletrolítico OFHC
- `be_c`: Berílio puro
- `cube_c`: Cobre-Berílio (Alloy 25)
- `ss304_c`: Aço Inoxidável austenítico AISI 304
- Tuplas estáticas NTuple correspondentes para avaliação polinomial com zero alocações na heap.
"""
struct NISTMaterialData
    cu_c::Matrix{Float64}
    be_c::Matrix{Float64}
    cube_c::Matrix{Float64}
    ss304_c::Matrix{Float64}

    # Tuplas estáticas NTuple pré-computadas (zero alocação com evalpoly)
    cu_tc_num::NTuple{5, Float64}
    cu_tc_den::NTuple{5, Float64}
    cu_sh::NTuple{9, Float64}
    cu_le::NTuple{5, Float64}

    ss304_tc::NTuple{9, Float64}
    ss304_sh::NTuple{8, Float64}
    ss304_ym::NTuple{5, Float64}
    ss304_le::NTuple{5, Float64}

    cube_tc::NTuple{8, Float64}
    cube_le::NTuple{5, Float64}

    be_sh::NTuple{9, Float64}
    be_le::NTuple{5, Float64}
end

"""
    load_nist_data(base_path)

Lê os arquivos CSV da pasta 'Material Properties' e converte de forma segura todas
as entradas em matrizes e tuplas estáticas de Float64.
"""
function load_nist_data(base_path::String = dirname(dirname(@__DIR__)))
    mat_dir = joinpath(base_path, "Material Properties")
    
    cu_path    = joinpath(mat_dir, "CU-TC1_TC2+TC3_TC4_TC5_SH_EC.csv")
    be_path    = joinpath(mat_dir, "BE-SH_LE1_LE2_LE3.csv")
    cube_path  = joinpath(mat_dir, "CuBe-TC_LE.csv")
    ss304_path = joinpath(mat_dir, "SS304-TC_SH_YM_YM_LE.csv")
    
    # Leitura bruta com DelimitedFiles
    cu_data    = readdlm(cu_path, ',')
    be_data    = readdlm(be_path, ',')
    cube_data  = readdlm(cube_path, ',')
    ss304_data = readdlm(ss304_path, ',')
    
    # Parser robusto para elementos de tabela
    function to_float(x)
        if x isa Number
            return Float64(x)
        end
        s = strip(string(x))
        if isempty(s)
            return 0.0 # Células sem coeficiente equivalem a zero no polinômio
        end
        v = tryparse(Float64, s)
        return v === nothing ? 0.0 : v
    end
    
    # Recorte das linhas 2 a 10 (coeficientes a..i) e colunas numéricas de cada material
    cu_c    = to_float.(cu_data[2:10, 2:8])
    be_c    = to_float.(be_data[2:10, 2:5])
    cube_c  = to_float.(cube_data[2:10, 2:3])
    ss304_c = to_float.(ss304_data[2:10, 2:6])
    
    # Extração de tuplas estáticas NTuple para avaliação com zero alocações na heap
    cu_tc_num = (cu_c[1, 1], cu_c[3, 1], cu_c[5, 1], cu_c[7, 1], cu_c[9, 1])
    cu_tc_den = (1.0, cu_c[2, 1], cu_c[4, 1], cu_c[6, 1], cu_c[8, 1])
    cu_sh     = ntuple(i -> cu_c[i, 6], 9)
    # NIST SRM 736 Thermal Contraction for Copper OFHC: ΔL/L_293 * 1e5 = a + bT + cT^2 + dT^3 + eT^4
    cu_le     = (-323.38965975, -0.23311604, 0.00728486, 2.17126e-6, -3.91206e-8)

    ss304_tc  = ntuple(i -> ss304_c[i, 1], 9)
    ss304_sh  = ntuple(i -> ss304_c[i, 2], 8)
    ss304_ym  = ntuple(i -> ss304_c[i, 4], 5) # Coluna YM2: válida de 4 K até 300 K (sem divergência de alta T)
    ss304_le  = ntuple(i -> ss304_c[i, 5], 5)

    cube_tc   = ntuple(i -> cube_c[i, 1], 8)
    cube_le   = ntuple(i -> cube_c[i, 2], 5)

    be_sh     = ntuple(i -> be_c[i, 1], 9)
    be_le     = ntuple(i -> be_c[i, 4], 5) # Coluna LE Polycrystalline
    
    return NISTMaterialData(
        cu_c, be_c, cube_c, ss304_c,
        cu_tc_num, cu_tc_den, cu_sh, cu_le,
        ss304_tc, ss304_sh, ss304_ym, ss304_le,
        cube_tc, cube_le,
        be_sh, be_le
    )
end

end # module NIST_Data
