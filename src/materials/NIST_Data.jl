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

Estrutura imutável que encapsula as matrizes numéricas de coeficientes NIST para:
- `cu_c`: Cobre eletrolítico OFHC
- `be_c`: Berílio puro
- `cube_c`: Cobre-Berílio (Alloy 25)
- `ss304_c`: Aço Inoxidável austenítico AISI 304
"""
struct NISTMaterialData
    cu_c::Matrix{Float64}
    be_c::Matrix{Float64}
    cube_c::Matrix{Float64}
    ss304_c::Matrix{Float64}
end

"""
    load_nist_data(base_path)

Lê os arquivos CSV da pasta 'Material Properties' e converte de forma segura todas
as entradas em matrizes de Float64.

POR QUE A FUNÇÃO `to_float` É NECESSÁRIA?
Os arquivos CSV do NIST contêm células vazias ("") quando um polinômio tem menos termos
que o grau máximo (grau 8), além de cabeçalhos e aspas. A leitura padrão com `readdlm`
retorna `SubString{String}` para células vazias ou textos. A função `to_float` sanitiza
isso, convertendo strings numéricas em `Float64` e qualquer célula vazia/inválida em `0.0`,
evitando exceções do tipo `ArgumentError: cannot parse "" as Float64`.
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
    
    return NISTMaterialData(cu_c, be_c, cube_c, ss304_c)
end

end # module NIST_Data
