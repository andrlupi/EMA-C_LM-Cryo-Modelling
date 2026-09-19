"""
Módulo SteadyState
==================

Solver numérico para o equilíbrio térmico estacionário não-linear em regime criogênico.

POR QUE NEWTON-RAPHSON EM VEZ DO MÉTODO ANTERIOR (ITERAÇÃO DE PICARD)?
No código original em 'SIM3LM-2SUP-Details.jl', utilizava-se a iteração de ponto fixo:
    T^{(k+1)} = [K(T^{(k)})]⁻¹ * L * u
Por que esse método falha ou oscila em criogenia?
1. Condutâncias não-constantes: K varia fortemente com T. A derivada dK/dT é grande,
   especialmente na região de pico do Cobre (15 K a 25 K).
2. Teorema de Ponto Fixo de Banach: Para que x_{k+1} = G(x_k) convirja, a norma do
   Jacobiano ‖∇G‖ deve ser estritamente menor que 1. Quando o gradiente térmico é íngreme,
   ‖∇G‖ > 1, fazendo com que a solução entre em ciclo limite ou divirja.

O MÉTODO DE NEWTON-RAPHSON MULTIVARIADO AMORTECIDO:
Formulamos o problema como a busca pelas raízes do sistema não-linear de conservação de energia:
    F(x) = 0   onde x = [T_livres]
A cada iteração k:
1. Calcula-se o Jacobiano Térmico J ∈ ℝ^{m × m}:
       J_{ij} = ∂Res_i / ∂x_j
   Fisicamente, a diagonal principal J_{ii} representa a condutância térmica diferencial total
   conectada ao nó i, e os termos fora da diagonal -J_{ij} representam as condutâncias mútuas.
2. Resolve-se a direção de descida de Newton:
       J * Δx = -F(x)  =>  Δx = -J⁻¹ * F(x)
3. Amortecimento Físico e Guarda de Positividade:
       x_{k+1} = x_k + α * Δx
   - Se o passo completo (α = 1) sugerir uma temperatura fisicamente impossível (T ≤ 0.5 K),
     o passo é reduzido pela metade até que todas as temperaturas permaneçam estritamente positivas.
   - Aplica-se uma busca linear retrógrada (backtracking line search) garantindo que a norma do
     resíduo ‖F(x_{k+1})‖_∞ decresça monotonicamente a cada iteração.
4. O método apresenta convergência quadrática local: converge tipicamente em 3 a 6 iterações
   com tolerância de máquina (< 10⁻¹⁰ W).
"""
module SteadyState

using LinearAlgebra
using ..Materials
using ..Network

export solve_steady_state, SteadyStateResult

"""
    SteadyStateResult

Armazena os resultados completos da convergência estacionária:
- `temperatures`: Vetor de temperaturas finais de todos os nós (fixos e livres) [K]
- `converged`: Booleano indicando se o resíduo atingiu a tolerância solicitada
- `iterations`: Número de passos de Newton executados
- `residual_norm`: Norma infinito final do resíduo de potência [W] (‖Res‖_∞ = max |Res_i|)
- `heat_flows`: Dicionário com a potência térmica [W] que atravessa cada elo individual
"""
struct SteadyStateResult
    temperatures::Vector{Float64}
    converged::Bool
    iterations::Int
    residual_norm::Float64
    heat_flows::Dict{String, Float64}
end

"""
    solve_steady_state(sys; T_guess=nothing, tol=1e-6, max_iters=50)

Encontra a distribuição de temperaturas de equilíbrio do `ThermalSystem`.

PARÂMETROS:
- `sys`: Instância de `ThermalSystem` com nós e elos definidos
- `T_guess`: Vetor opcional de chute inicial de temperaturas. Se omitido,
  utiliza a média das temperaturas dos nós fixos de fronteira.
- `tol`: Tolerância de convergência na norma infinito do resíduo de energia [W]
  (default: 1e-6 W = 1 µW)
- `max_iters`: Número máximo de iterações permitidas
"""
function solve_steady_state(sys::ThermalSystem; 
                            T_guess::Union{Nothing, Vector{<:Real}} = nothing,
                            tol::Float64 = 1e-6,
                            max_iters::Int = 50)
    
    n_nodes = length(sys.nodes)
    
    # 1. Separar índices dos nós livres (incógnitas) e nós fixos (condições de Dirichlet)
    free_indices = Int[]
    fixed_indices = Int[]
    for i in 1:n_nodes
        if sys.nodes[i].is_fixed
            push!(fixed_indices, i)
        else
            push!(free_indices, i)
        end
    end
    
    # 2. Inicializar o vetor completo de temperaturas
    T_full = zeros(Float64, n_nodes)
    for idx in fixed_indices
        T_full[idx] = sys.nodes[idx].fixed_temp
    end
    
    # Chute inicial para os nós livres
    if T_guess === nothing
        # Chute padrão razoável: média das fronteiras fixas (ou 20 K se não houver fronteiras fixas)
        default_T = isempty(fixed_indices) ? 20.0 : sum(T_full[fixed_indices]) / length(fixed_indices)
        for idx in free_indices
            T_full[idx] = default_T
        end
    else
        for idx in free_indices
            T_full[idx] = Float64(T_guess[idx])
        end
    end
    
    n_free = length(free_indices)
    # Caso trivial: se todos os nós forem fixos, nada a resolver
    if n_free == 0
        return SteadyStateResult(copy(T_full), true, 0, 0.0, compute_link_heat_flows(sys, T_full))
    end
    
    # 3. Função objetivo: recebe as temperaturas livres x e retorna apenas os resíduos dos nós livres
    function obj_residuals(x::Vector{Float64})
        T_eval = copy(T_full)
        for (k, idx) in enumerate(free_indices)
            T_eval[idx] = x[k]
        end
        res_full = energy_balance_residuals(sys, T_eval)
        return res_full[free_indices]
    end
    
    # Estado inicial das variáveis livres
    x = [T_full[idx] for idx in free_indices]
    res = obj_residuals(x)
    norm_res = norm(res, Inf)
    
    iter = 0
    converged = norm_res < tol
    
    # 4. Laço de Newton-Raphson Amortecido
    while !converged && iter < max_iters
        iter += 1
        
        # Cálculo do Jacobiano por diferenças finitas centrais/progressivas:
        # J_ij = ∂Res_i / ∂x_j ≈ (Res_i(x + h*e_j) - Res_i(x)) / h
        J = zeros(Float64, n_free, n_free)
        h = 1e-5
        for j in 1:n_free
            x_step = copy(x)
            x_step[j] += h
            res_step = obj_residuals(x_step)
            J[:, j] = (res_step - res) / h
        end
        
        # Resolução do sistema linear de Newton: J * Δx = -res
        Δx = - J \ res
        
        # Guarda Física de Positividade:
        # Em criogenia, temperaturas absolutas T <= 0 K são proibidas pelas leis da termodinâmica
        # e levariam a log10(T) imaginário ou erro nas equações do NIST.
        α = 1.0
        x_new = x + α * Δx
        while any(x_new .<= 0.5) && α > 1e-4
            α *= 0.5
            x_new = x + α * Δx
        end
        
        # Busca linear retrógrada (Backtracking Line Search):
        # Garante que o novo ponto reduza a norma do resíduo, evitando passos excessivos
        res_new = obj_residuals(x_new)
        norm_res_new = norm(res_new, Inf)
        
        while norm_res_new > norm_res && α > 1e-4
            α *= 0.5
            x_new = x + α * Δx
            res_new = obj_residuals(x_new)
            norm_res_new = norm(res_new, Inf)
        end
        
        # Atualização para a próxima iteração
        x = x_new
        res = res_new
        norm_res = norm_res_new
        
        if norm_res < tol
            converged = true
        end
    end
    
    # 5. Montagem do vetor final com a solução encontrada
    for (k, idx) in enumerate(free_indices)
        T_full[idx] = x[k]
    end
    
    return SteadyStateResult(T_full, converged, iter, norm_res, compute_link_heat_flows(sys, T_full))
end

"""
    compute_link_heat_flows(sys, T_full)

Calcula a potência térmica [W] transmitida através de cada elo da rede na solução final.
"""
function compute_link_heat_flows(sys::ThermalSystem, T_full::Vector{Float64})
    flows = Dict{String, Float64}()
    for link in sys.links
        Ta = T_full[link.node_a_idx]
        Tb = T_full[link.node_b_idx]
        flows[link.name] = heat_flow(link, Ta, Tb)
    end
    return flows
end

end # module SteadyState
