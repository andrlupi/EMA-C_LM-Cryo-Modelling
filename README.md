# EMA-C_LM-Cryo-Modelling: Modelagem Criogênica e Co-Design Termo-Mecânico

Modelagem térmica por parâmetros concentrados (*Lumped-Parameter Thermal Network*) e co-design termo-mecânico da **nanoestação da linha de luz EMA do Sirius / LNLS (CNPEM)**.

---

## 1. Visão Geral do Projeto

Na nanoestação da linha EMA, amostras em condições extremas de pressão e temperatura são montadas dentro de uma célula de alta pressão de diamante (**Mini-DAC**, em Cobre-Berílio ou Inox). O desafio de engenharia consiste no **compromisso termo-mecânico**:

* **Requisito Mecânico:** Frequência natural de ressonância lateral $f_n \ge 150\text{ Hz}$ para rejeitar vibrações do criocooler e do piso do sincrotron, garantindo estabilidade submicrométrica de feixe.
* **Requisito Térmico:** Manter a Mini-DAC abaixo de $6.0\text{ K}$ (idealmente sub-5 K) sob aporte térmico do feixe síncrotron de raios X ($q_{\text{beam}} = 10\text{ mW}$), condução mecânica dos suportes de sustentação e radiação do escudo de $40\text{ K}$.

---

## 2. Arquitetura Modular do Código (`src/`)

O projeto foi totalmente reestruturado em Julia como uma biblioteca física modular:

```text
EMA-C_LM-Cryo-Modelling/
├── Project.toml & Manifest.toml  # Ambiente reproduzível do Julia
├── .gitignore                    # Filtro de artefatos de sistema e IDEs
│
├── src/
│   ├── CryoThermal.jl            # Ponto de entrada e exportação da API pública
│   ├── materials/
│   │   ├── NIST_Data.jl          # Carregamento e parser seguro das tabelas do NIST
│   │   └── Materials.jl          # Cu OFHC, Inox 304, CuBe, Be e Ti-6Al-4V + Quadratura de Gauss
│   ├── network/
│   │   └── Network.jl            # Grafos térmicos: nós (ThermalNode) e elos (ConductionLink, etc.)
│   ├── interfaces/
│   │   └── ColdInterfaces.jl     # Junta com folha de Índio e Gás de Troca de Hélio (Sherman-Lees)
│   └── solvers/
│       └── SteadyState.jl        # Solver Newton-Raphson com busca linear e guarda de positividade
│
├── test/
│   └── runtests.jl               # Suíte completa com 40 testes unitários de validação física
│
├── scripts/
│   ├── plot_cold_interfaces.jl   # Estudo comparativo das interfaces do dedo frio
│   └── plot_codesign_tradeoff.jl # Estudo da fronteira de Pareto termo-mecânica
│
├── figures/                      # Figuras vetoriais de alta resolução geradas com CairoMakie
│   ├── cold_interfaces_study.png
│   └── codesign_pareto_study.png
│
└── Material Properties/          # Dados brutos empíricos do NIST (CSVs)
```

---

## 3. Fundamentação Física e Numérica

### A. Transformada de Kirchhoff (Condução Exata)
Em temperaturas criogênicas, a condutividade térmica $k(T)$ varia bruscamente. O fluxo de calor condutivo unidimensional é resolvido sem aproximação linear pela **Integral de Condutividade Térmica**:
$$Q_{a \to b} = \frac{A}{L} \int_{T_b}^{T_a} k(T) \, dT = \frac{A}{L} \big[ \Theta(T_a) - \Theta(T_b) \big]$$
A integral é calculada numericamente por **Quadratura de Gauss-Legendre de 7 pontos**, garantindo exatidão analítica para polinômios de grau até 13 com apenas 7 chamadas de função.

### B. Solver Não-Linear por Newton-Raphson Amortecido
A conservação de energia em regime permanente (Primeira Lei da Termodinâmica) estabelece:
$$R_i(\mathbf{T}) = Q_{\text{ext}, i} + \sum_{j} Q_{j \to i}(\mathbf{T}) = 0 \quad \forall i \in \text{nós livres}$$
O solver calcula o Jacobiano térmico $J_{ij} = \frac{\partial R_i}{\partial T_j}$ e aplica passos de Newton com guarda física de positividade ($T > 0.5\text{ K}$) e busca linear retrógrada (*backtracking*), atingindo convergência quadrática com resíduos menores que $10^{-10}\text{ W}$ em 4 a 6 iterações.

---

## 4. Principais Resultados de Engenharia

### A. Interface Fria: Cordoalha com Índio vs Gás de Troca de Hélio
* **Folha de Índio:** Sob aperto de parafusos com pressão $P \ge 3\text{ MPa}$, a folha de índio escoa plasticamente e aumenta a condutância de contato em **$21\times$** em relação ao contato metal-metal seco, mantendo a DAC em **$4.9\text{ K}$**.
* **Gás de Troca ($^4\text{He}$):** Elimina 100% das vibrações mecânicas da bomba, mas sua condutividade térmica limita a temperatura de equilíbrio a **$\sim 11.6\text{ K}$** para uma folga de $2\text{ mm}$. A cordoalha com índio é essencial para operação sub-5 K.

### B. Co-Design Termo-Mecânico e Seleção de Materiais
* **Flexão Lateral:** Para um tripé isostático de 3 hastes, a rigidez lateral $k_{\text{flex}} \propto \frac{12 E I}{L^3}$ domina a estabilidade.
* **Titânio vs Inox:** O Titânio Ti-6Al-4V possui uma Figura de Mérito ($\text{FOM} = \frac{E}{\int k dT}$) **$52\%$ superior** à do Inox 304, mantendo a DAC em **$4.42\text{ K}$** mesmo a $400\text{ Hz}$.
* **Recomendação de Fabricação:** Hastes de Inox 304 com $d = 1.5\text{ mm}$ e $L = 30\text{ mm}$ garantem com folga $f_n \approx 185\text{ Hz}$ e $T_{\text{DAC}} \approx 4.55\text{ K}$.

---

## 5. Como Executar

### Pré-requisitos
* Julia $\ge$ 1.10 (testado em Julia 1.12.3).

### Executar os Testes Unitários
```bash
julia --project=. test/runtests.jl
```

### Gerar os Gráficos Científicos (CairoMakie)
```bash
julia --project=. scripts/plot_cold_interfaces.jl
julia --project=. scripts/plot_codesign_tradeoff.jl
```
As figuras serão salvas na pasta `figures/`.
