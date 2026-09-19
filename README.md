# EMA-C_LM-Cryo-Modelling: Modelagem Criogênica e Co-Design Termo-Mecânico

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)
[![Julia Version](https://img.shields.io/badge/Julia-1.10%2B-purple.svg)](https://julialang.org)
[![Tests](https://img.shields.io/badge/Tests-44%2F44%20Passing-brightgreen.svg)](test/runtests.jl)

Modelagem térmica por parâmetros concentrados (*Lumped-Parameter Thermal Network*) e co-design termo-mecânico da **nanoestação da linha de luz EMA do Sirius / LNLS (CNPEM)**.

> 📄 **Relatório Técnico Completo:** Para a fundamentação matemática, deduções e análise física aprofundada, consulte o documento [`docs/RELATORIO_TECNICO_REDESIGN.md`](docs/RELATORIO_TECNICO_REDESIGN.md).

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
├── CITATION.cff                      # Metadados de citação científica (BibTeX / APA)
├── LICENSE                           # Licença BSD 3-Clause
├── Project.toml & Manifest.toml      # Ambiente reproduzível do Julia
├── README.md                         # Visão geral do repositório
│
├── docs/
│   └── RELATORIO_TECNICO_REDESIGN.md # Relatório técnico completo de engenharia
│
├── src/
│   ├── CryoThermal.jl                # Ponto de entrada e exportação da API pública
│   ├── materials/
│   │   ├── NIST_Data.jl              # Carregamento e parser seguro das tabelas do NIST
│   │   └── Materials.jl              # Cu OFHC, Inox 304, CuBe, Be e Ti-6Al-4V + Gauss-Legendre
│   ├── network/
│   │   └── Network.jl                # Grafos térmicos: nós (ThermalNode) e elos (ConductionLink, etc.)
│   ├── interfaces/
│   │   └── ColdInterfaces.jl         # Junta com folha de Índio e Gás de Troca de Hélio (Sherman-Lees)
│   └── solvers/
│       ├── SteadyState.jl            # Solver Newton-Raphson com busca linear e guarda de positividade
│       └── Transient.jl              # Solver transiente de Rosenbrock L-estável com passo adaptativo
│
├── test/
│   └── runtests.jl                   # Suíte completa com 44 testes unitários de validação física
│
├── scripts/
│   ├── plot_cold_interfaces.jl       # Estudo comparativo das interfaces do dedo frio
│   ├── plot_codesign_tradeoff.jl     # Estudo da fronteira de Pareto termo-mecânica
│   └── plot_transient_cooldown.jl    # Simulação temporal transiente (300 K -> 4.2 K)
│
├── figures/                          # Figuras vetoriais de alta resolução geradas com CairoMakie
│   ├── cold_interfaces_study.png     # Estudo das interfaces frias (Índio vs Gás He)
│   ├── codesign_pareto_study.png     # Estudo da fronteira de Pareto e mapa de projeto (d vs L)
│   └── transient_cooldown_study.png  # Dinâmica transiente global, zoom criogênico e curva de cp(T)
│
└── Material Properties/              # Dados brutos empíricos do NIST (CSVs)
```

---

## 3. Fundamentação Física e Numérica

### A. Transformada de Kirchhoff (Condução Exata)
Em temperaturas criogênicas, a condutividade térmica $k(T)$ varia bruscamente. O fluxo de calor condutivo unidimensional é resolvido sem aproximação linear pela **Integral de Condutividade Térmica**:
$$Q_{a \to b} = \frac{A}{L} \int_{T_b}^{T_a} k(T) \, dT = \frac{A}{L} \big[ \Theta(T_a) - \Theta(T_b) \big]$$
A integral é calculada numericamente por **Quadratura de Gauss-Legendre de 7 pontos**, garantindo exatidão analítica para polinômios de grau até 13 com apenas 7 chamadas de função.

### B. Solver Estacionário por Newton-Raphson Amortecido
A conservação de energia em regime permanente (Primeira Lei da Termodinâmica) estabelece:
$$R_i(\mathbf{T}) = Q_{\text{ext}, i} + \sum_{j} Q_{j \to i}(\mathbf{T}) = 0 \quad \forall i \in \text{nós livres}$$
O solver calcula o Jacobiano térmico $J_{ij} = \frac{\partial R_i}{\partial T_j}$ e aplica passos de Newton com guarda física de positividade ($T > 0.5\text{ K}$) e busca linear retrógrada (*backtracking*), atingindo convergência quadrática com resíduos menores que $10^{-10}\text{ W}$ em 4 a 6 iterações.

### C. Solver Transiente Não-Linear de Resfriamento ($300\text{ K} \to 4.2\text{ K}$)
Para superar as limitações dos modelos LTI com capacidade constante, o solver integra:
$$\frac{dT_i}{dt} = \frac{\text{Res}_i(\mathbf{T})}{M_i \cdot c_{p,i}(T_i)}$$
utilizando um método **linearmente implícito de Rosenbrock ($L$-estável)** com controle adaptativo de passo:
$$(\text{diag}(\mathbf{C}) - \Delta t \cdot \mathbf{J}) \cdot \Delta \mathbf{T} = \Delta t \cdot \mathbf{Res}$$
Essa formulação é incondicionalmente estável contra a rigidez numérica gerada pela queda de mais de $3500\times$ no calor específico $c_p(T)$ Debye do cobre entre $300\text{ K}$ e $4.2\text{ K}$.

---

## 4. Principais Resultados de Engenharia

### A. Interface Fria: Cordoalha com Índio vs Gás de Troca de Hélio
* **Folha de Índio:** Sob aperto de parafusos com pressão $P \ge 2.5\text{ MPa}$, a folha de índio escoa plasticamente e aumenta a condutância de contato em **$21\times$** em relação ao contato metal-metal seco, mantendo a DAC em **$4.93\text{ K}$** sob $10\text{ mW}$ de feixe.
* **Gás de Troca ($^4\text{He}$):** Elimina 100% das vibrações mecânicas da bomba, mas sua condutividade térmica limita a temperatura de equilíbrio a **$\sim 11.6\text{ K}$** para uma folga de $2\text{ mm}$. A cordoalha com índio é indispensável para operação sub-5 K.

### B. Co-Design Termo-Mecânico e Seleção de Materiais
* **Flexão Lateral:** Para um tripé isostático de 3 hastes, a rigidez lateral $k_{\text{flex}} \propto \frac{12 E I}{L^3}$ domina a estabilidade vibracional.
* **Titânio vs Inox:** O Titânio Ti-6Al-4V possui uma Figura de Mérito ($\text{FOM} = \frac{E}{\int k dT}$) **$52\%$ superior** à do Inox 304, mantendo a DAC em **$4.42\text{ K}$** mesmo a $400\text{ Hz}$.
* **Recomendação de Fabricação:** Hastes de Inox 304 com $d = 1.5\text{ mm}$ e $L = 30\text{ mm}$ garantem com folga $f_n \approx 185\text{ Hz}$ e $T_{\text{DAC}} \approx 4.55\text{ K}$.

### C. Tempo de Resfriamento da Estação
* O tempo de resfriamento total é ditado pela taxa de extração do criocooler ($\tau \approx 30\text{ min}$).
* O sistema cruza o patamar de $77\text{ K}$ (LN2) em **$\approx 42\text{ min}$** e estabiliza em **$5.00\text{ K}$** em **$\approx 2.5\text{ horas}$**, com um gradiente residual na sapata de apenas $+0.07\text{ K}$.

---

## 5. Como Executar

### Pré-requisitos
* Julia $\ge$ 1.10 (testado e validado em Julia 1.12.3).

### Executar os Testes Unitários
```bash
julia --project=. test/runtests.jl
```

### Gerar os Gráficos Científicos (CairoMakie)
```bash
# 1. Interfaces Frias (Cordoalha c/ Índio vs Gás He):
julia --project=. scripts/plot_cold_interfaces.jl

# 2. Co-Design Termo-Mecânico e Fronteira de Pareto:
julia --project=. scripts/plot_codesign_tradeoff.jl

# 3. Dinâmica Transiente de Resfriamento (300 K -> 4.2 K):
julia --project=. scripts/plot_transient_cooldown.jl
```
As figuras de alta resolução vetorial são geradas na pasta `figures/`.

---

## 6. Citação Acadêmica

Se você utilizar este código, modelos ou metodologias em seu trabalho, consulte o arquivo [`CITATION.cff`](CITATION.cff) ou utilize o botão **"Cite this repository"** do GitHub:

```bibtex
@software{Lupi_EMA_Cryo_2026,
  author = {Lupi, André},
  title = {{EMA-C_LM-Cryo-Modelling: Cryogenic Lumped-Mass Modeling and Thermo-Mechanical Co-Design for the EMA Beamline Nano-Station at Sirius/LNLS}},
  year = {2026},
  publisher = {GitHub},
  version = {1.0.0},
  license = {BSD-3-Clause}
}
```

---

## 7. Declaração de Uso de Inteligência Artificial (Conformidade UNICAMP CONSU-A-005/2026)

Em conformidade com a **[Deliberação CONSU-A-005/2026](https://www.pg.unicamp.br/norma/32327/0)** da UNICAMP (publicada no D.O.E. em 14/04/2026), declara-se expressamente que:
* **Ferramenta Utilizada:** Google Antigravity (modelos Gemini / Google DeepMind) foi utilizado como ferramenta auxiliar de programação em par (*pair programming*) para refatoração de código Julia, construção da suíte de testes unitários e auxílio na formatação de documentação.
* **Autoria e Supervisão Humana:** A concepção intelectual, os requisitos de engenharia criogênica da linha EMA, as condições de contorno físicas, os dados termofísicos do NIST e o julgamento crítico de todos os resultados são de autoria e responsabilidade integral do pesquisador humano responsável (conforme Art. 3º, II e Art. 5º).
* Para a declaração formal detalhada, consulte a [Seção 9 do Relatório Técnico](docs/RELATORIO_TECNICO_REDESIGN.md#9-declara%C3%A7%C3%A3o-de-uso-de-intelig%C3%AAncia-artificial-generativa-conformidade-unicamp-consu-a-0052026).

