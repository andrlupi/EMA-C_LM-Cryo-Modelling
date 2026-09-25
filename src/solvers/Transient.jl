"""
Módulo Transient
================

Solver numérico de alta estabilidade para a dinâmica térmica transiente não-linear
(curva de resfriamento / cooldown) de sistemas criogênicos.

POR QUE O MODELO LINEAR ANTERIOR (LTI / lsim) FALHAVA?
No código legado, utilizava-se `ControlSystems.lsim(ss(A, B, C, D))` com matrizes constantes.
Em física de baixas temperaturas:
1. O calor específico de rede cp(T) varia por mais de 3 ordens de grandeza entre 300 K e 4.2 K (Debye: cp ∝ T³).
2. A condutividade térmica k(T) varia de forma fortemente não-linear com picos em 20-30 K.
3. Um modelo LTI com cp constante superestima drasticamente o tempo de resfriamento a frio ou
   subestima a inércia térmica a quente, sendo fisicamente inválido para transições 300 K -> 4 K.

O DESAFIO DA RIGIDEZ NUMÉRICA (STIFFNESS) EM CRIOGENIA:
Conforme o sistema se aproxima de 4.2 K, a capacidade térmica C_i = M_i * cp_i(T_i) despenca para
valores ínfimos (~10⁻⁴ a 10⁻² J/K). A constante de tempo térmica local:
    τ_i = C_i / G_total_i
colapsa de dezenas de segundos para milissegundos. Integradores explícitos clássicos (como RK4 ou Heun)
exigem Δt < 2*τ_min para manter estabilidade; se Δt exceder essa fração de milissegundo, a solução
oscila artificialmente em torno do equilíbrio.

SOLUÇÃO ADOTADA: INTEGRADOR LINEARMENTE IMPLÍCITO DE ROSENBROCK (L-ESTÁVEL)
Resolvemos a equação diferencial ordinária não-linear:
    C_i(T_i) * (dT_i / dt) = Res_i(T)
Discretizando de forma linearmente implícita (método de Rosenbrock / W-method):
    (diag(C) - Δt * J) * ΔT = Δt * Res
Onde:
- C é o vetor de capacidades térmicas [M_i * cp_i(T_i)] dos nós livres.
- J = ∂Res/∂T é a matriz Jacobiana calculada exatamente via ForwardDiff.
- A matriz A = (diag(C) - Δt * J) é estritamente definida positiva e bem condicionada.
"""
module Transient

using LinearAlgebra
using ForwardDiff
using ..Materials
using ..Network

export solve_transient,
       TransientResult,
       exponential_cryocooler_cooldown

"""
    TransientResult

Estrutura imutável com a série temporal completa da simulação de resfriamento:
- `times`: Vetor de instantes de tempo [s]
- `temperatures`: Matriz [N_nós × N_passos] de temperaturas [K]
- `node_names`: Nomes descritivos de cada nó térmico
"""
struct TransientResult
    times::Vector{Float64}
    temperatures::Matrix{Float64}
    node_names::Vector{String}
end

"""
    exponential_cryocooler_cooldown(t; T_start=300.0, T_final=4.2, tau=1800.0)

Modelo empírico exponencial representativo da curva de potência e resfriamento
do 2º estágio de um criocooler de ciclo fechado (Gifford-McMahon ou Pulse Tube):
    T_cold(t) = T_final + (T_start - T_final) * exp(-t / tau)
"""
function exponential_cryocooler_cooldown(t::Real; T_start::Real = 300.0, T_final::Real = 4.2, tau::Real = 1800.0)
    return T_final + (T_start - T_final) * exp(-t / tau)
end

"""
    solve_transient(sys, t_span; T_init=300.0, cold_head_fn=nothing, boundary_fns=nothing, dt_init=1.0, tol=1e-3, max_steps=50000)

Executa a integração transiente do resfriamento criogênico via método linearmente implícito adaptativo.
"""
function solve_transient(sys::ThermalSystem, t_span::Tuple{<:Real, <:Real};
                         T_init::Real = 300.0,
                         cold_head_fn::Union{Nothing, Function} = nothing,
                         boundary_fns::Union{Nothing, Dict{Int, <:Function}} = nothing,
                         dt_init::Real = 1.0,
                         tol::Float64 = 1e-3,
                         max_steps::Int = 50000)
    
    t_start, t_end = Float64(t_span[1]), Float64(t_span[2])
    n_nodes = length(sys.nodes)
    
    # Mapeamento de nós livres e fixos
    free_indices = Int[]
    fixed_indices = Int[]
    for i in 1:n_nodes
        if sys.nodes[i].is_fixed
            push!(fixed_indices, i)
        else
            push!(free_indices, i)
        end
    end
    n_free = length(free_indices)
    
    # Construção do mapa de funções de contorno
    b_fns = Dict{Int, Function}()
    if boundary_fns !== nothing
        for (k, v) in boundary_fns
            b_fns[k] = v
        end
    elseif cold_head_fn !== nothing && !isempty(fixed_indices)
        # Aplica a função ao primeiro nó fixo por convenção
        b_fns[fixed_indices[1]] = cold_head_fn
    end
    
    # Função auxiliar para atualizar temperaturas dos nós fixos em um tempo t
    function update_fixed_nodes!(T_vec::Vector{Float64}, t_val::Float64)
        for idx in fixed_indices
            if haskey(b_fns, idx)
                T_vec[idx] = Float64(b_fns[idx](t_val))
            else
                T_vec[idx] = sys.nodes[idx].fixed_temp
            end
        end
    end
    
    # Caso especial: se não há nós livres, o sistema é puramente ditado pelas condições de contorno
    if n_free == 0
        times = [t_start, t_end]
        T_matrix = zeros(Float64, n_nodes, 2)
        T_start_vec = zeros(Float64, n_nodes)
        update_fixed_nodes!(T_start_vec, t_start)
        T_end_vec = zeros(Float64, n_nodes)
        update_fixed_nodes!(T_end_vec, t_end)
        T_matrix[:, 1] = T_start_vec
        T_matrix[:, 2] = T_end_vec
        return TransientResult(times, T_matrix, [node.name for node in sys.nodes])
    end
    
    # Estado inicial
    T_current = fill(Float64(T_init), n_nodes)
    update_fixed_nodes!(T_current, t_start)
    
    # Histórico de gravação
    time_history = Float64[t_start]
    temp_history = Vector{Float64}[copy(T_current)]
    
    # Função que resolve o passo de Rosenbrock com ForwardDiff
    function take_rosenbrock_step(T_state::Vector{Float64}, t_val::Float64, h::Float64)
        T_eval = copy(T_state)
        update_fixed_nodes!(T_eval, t_val)
        
        # 1. Capacidades térmicas nos nós livres
        C_free = Float64[
            max(sys.nodes[idx].mass * specific_heat(sys.nodes[idx].material, max(T_eval[idx], 0.5)), 1e-6)
            for idx in free_indices
        ]
        
        # 2. Resíduos e Jacobiano exato via ForwardDiff nos nós livres
        T_fixed_eval = copy(T_eval)
        function free_res(x_free::AbstractVector{T}) where {T<:Real}
            T_work = Vector{T}(undef, n_nodes)
            for idx in fixed_indices
                T_work[idx] = T_fixed_eval[idx]
            end
            for (k, idx) in enumerate(free_indices)
                T_work[idx] = x_free[k]
            end
            return energy_balance_residuals(sys, T_work)[free_indices]
        end
        
        x_curr = Float64[T_eval[idx] for idx in free_indices]
        res_free = free_res(x_curr)
        J = ForwardDiff.jacobian(free_res, x_curr)
        
        # 3. Solução do sistema linear (C - h * J) * ΔT = h * Res com proteção contra singularidade
        A = Diagonal(C_free) - h * J
        b = h * res_free
        delta_T_free = try
            A \ b
        catch e
            if e isa LinearAlgebra.SingularException || e isa LinearAlgebra.LAPACKException
                @warn "Transient solver matrix (C - h*J) is singular at t = $t_val s."
                zeros(Float64, n_free)
            else
                rethrow(e)
            end
        end
        
        return delta_T_free
    end
    
    t = t_start
    dt = Float64(dt_init)
    step_count = 0
    
    # -------------------------------------------------------------------------
    # LAÇO PRINCIPAL DE INTEGRAÇÃO ADAPTATIVA
    # -------------------------------------------------------------------------
    while t < t_end && step_count < max_steps
        step_count += 1
        if t + dt > t_end
            dt = t_end - t
        end
        
        # Passo 1: Passo completo de tamanho dt de t até t + dt
        delta_full = take_rosenbrock_step(T_current, t + dt, dt)
        
        # Passo 2: Dois meios-passos de tamanho dt/2 com sincronização rigorosa de contorno
        dt_half = 0.5 * dt
        delta_half1 = take_rosenbrock_step(T_current, t + dt_half, dt_half)
        T_mid = copy(T_current)
        T_mid[free_indices] += delta_half1
        # Sincronização explícita das fronteiras em t + dt_half
        update_fixed_nodes!(T_mid, t + dt_half)
        
        delta_half2 = take_rosenbrock_step(T_mid, t + dt, dt_half)
        T_half = copy(T_mid)
        T_half[free_indices] += delta_half2
        # Sincronização explícita das fronteiras em t + dt
        update_fixed_nodes!(T_half, t + dt)
        
        # Estimativa de erro local na norma infinito
        T_predicted_full = copy(T_current)
        T_predicted_full[free_indices] += delta_full
        err_local = maximum(abs.(T_half[free_indices] - T_predicted_full[free_indices]))
        
        # Controle adaptativo do passo
        if err_local <= tol || dt <= 1e-4
            # Passo aceito: adotamos o resultado dos dois meios-passos
            t += dt
            T_current = copy(T_half)
            for i in free_indices
                T_current[i] = max(T_current[i], 0.5)
            end
            
            push!(time_history, t)
            push!(temp_history, copy(T_current))
            
            # Ajuste de passo para a próxima iteração
            factor = err_local > 0 ? 0.9 * (tol / err_local)^0.33 : 2.0
            dt = clamp(dt * factor, 1e-4, 120.0)
        else
            # Passo rejeitado: reduz o passo e repete
            factor = 0.9 * (tol / err_local)^0.33
            dt = max(dt * factor, 1e-4)
        end
    end
    
    # Montagem da matriz final [N_nós × N_passos]
    n_records = length(time_history)
    T_matrix = zeros(Float64, n_nodes, n_records)
    for step in 1:n_records
        T_matrix[:, step] = temp_history[step]
    end
    
    node_names = [node.name for node in sys.nodes]
    return TransientResult(time_history, T_matrix, node_names)
end

end # module Transient
