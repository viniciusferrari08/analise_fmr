# Software de Análise FMR para Filmes Finos

Software Python abrangente para análise de **Ressonância Ferromagnética (FMR)** em filmes magnéticos finos. Este pacote implementa a equação de Kittel para extrair parâmetros magnéticos (Ms, Ha, Hk) de medidas experimentais de frequência vs. campo de ressonância, incluindo análise automática de largura de linha e ajuste espectral completo.

## Visão Geral

Este software foi desenvolvido como parte do projeto de pesquisa **PIIC/UFES 2025** "Introdução à Ressonância Magnética Eletrônica" sob supervisão do Prof. Thiago Eduardo Pedreira Bueno. Fornece ferramentas de modelagem teórica e análise computacional para estudar a dinâmica de magnetização em filmes finos ferromagnéticos usando espectroscopia FMR.

## Características

### ✅ Funcionalidades Principais
- **Análise FMR Completa**: Pipeline de processamento completo para dados experimentais
- **Extração Automática de Campo de Ressonância**: Ajuste avançado de derivada Lorentziana
- **Ajuste da Equação de Kittel**: Geometrias no plano e perpendicular
- **Análise de Largura de Linha**: Extração automática de larguras de linha FMR (ΔH)
- **Interface Gráfica**: Interface gráfica amigável para processamento de dados
- **Análise Estatística**: Propagação de erros e cálculos de R²

### 🧬 Implementação Física
- **Equação de Kittel**: Implementação completa para geometrias de filmes finos
- **Parâmetros Magnéticos**:
  - **Ms**: Magnetização de saturação
  - **Ha**: Campo de anisotropia no plano
  - **Hk**: Campo de anisotropia perpendicular
- **Geometrias**: Configurações perpendicular (θ = 90°) e no plano (θ = 0°)
- **Razão Giromagnética**: γ = 2.8×10¹⁰ Hz/T para spins eletrônicos

## Início Rápido

### Instalação

```bash
# Clone o repositório
git clone https://github.com/username/fmr-thin-film-analysis.git
cd fmr-thin-film-analysis

# Instale as dependências
pip install -r requirements.txt
```

### Uso

#### Interface Gráfica (Recomendada)
```bash
python fmr_gui.py
```

#### Interface de Linha de Comando
```bash
# Análise completa
python analise_fmr_completa.py

# Análise de espectro individual
python extrair_campo_ressonancia.py
```

#### API Python
```python
from extrair_campo_ressonancia import FMRSpectrumAnalyzer
from fmr_fitting import FMRFitting

# Carregar e processar espectros
analyzer = FMRSpectrumAnalyzer()
results = analyzer.process_measurement_folder("pasta_dados")

# Ajustar parâmetros magnéticos
fitter = FMRFitting()
params, fitted_values = fitter.fit_inplane_data(frequencies, fields)
```

## Fundamentação Teórica

### Implementação da Equação de Kittel

O software implementa a equação de Kittel completa para ressonância ferromagnética em filmes finos, conforme descrito no framework teórico do relatório final PIIC/UFES 2025.

#### Geometria Perpendicular (θ = 90°)
```
f = (γ/2π) × √[(Hr + Ha + Hk)(Hr + Ha + Hk + Ms)]
```

#### Geometria no Plano (θ = 0°)
```
f = (γ/2π) × √[(Hr + Ha)(Hr + Ha + Ms - Hk)]
```

Onde:
- **γ = 2.8×10¹⁰ Hz/T**: Razão giromagnética para spins eletrônicos
- **Hr**: Campo de ressonância (extraído dos espectros experimentais)
- **Ms**: Magnetização de saturação (parâmetro ajustado)
- **Ha**: Campo de anisotropia no plano (parâmetro ajustado)
- **Hk**: Campo de anisotropia perpendicular (parâmetro ajustado)

### Ajuste de Derivada Lorentziana

Para espectros de absorção derivativos (formato dχ"/dH), o software ajusta a forma matematicamente correta:

```
dL/dH = -2A(H-Hr) / [ΔH² × (1 + ((H-Hr)/ΔH)²)²] + offset
```

Esta é a derivada exata da Lorentziana normalizada `L(H) = A / (1 + ((H-Hr)/ΔH)²)`, que permite extração precisa do campo de ressonância (Hr) e largura de linha (ΔH) com controle de qualidade estatístico (limite R² > 0.7).

## Estrutura do Projeto

```
📁 Módulos Principais
├── extrair_campo_ressonancia.py    # Classe FMRSpectrumAnalyzer - processamento de espectros
├── fmr_fitting.py                  # Classe FMRFitting - ajuste da equação de Kittel
├── fmr_gui.py                      # Interface gráfica do usuário
├── analise_fmr_completa.py         # Pipeline de análise completa
├── requirements.txt                # Dependências Python
├── README.md                       # Este arquivo
└── CLAUDE.md                       # Guia de desenvolvimento Claude Code

📁 Estrutura de Dados
├── EAFExp2025/                     # Raiz dos dados experimentais
│   ├── NiFe_Cu_6nm/medidas selec/  # Medidas selecionadas NiFe/Cu
│   └── [outras amostras]/          # Conjuntos de dados experimentais adicionais
└── NiFe/                           # Arquivos de dados da amostra NiFe

📁 Arquivos de Saída
├── campos_ressonancia_*.txt        # Campos de ressonância extraídos
├── resultado_analise_completa.txt  # Relatório de análise completa
└── ajuste_fmr_final.png           # Visualização do ajuste final
```

## Formatos de Dados

### Arquivos de Dados Experimentais (.dat)
**Formato de 4 colunas separadas por espaço:**
```
freq(GHz)  corrente(mA)  campo(Oe)  sinal_fmr(u.a.)
5.0        100           300        -0.0234
5.0        100           301        -0.0189
5.0        100           302        -0.0145
...
```

### Formato de Entrada Simples (.txt, .csv)
**2 colunas frequência vs. campo:**
```
# Frequencia(GHz)  Campo(mT)
8.0               150.2
8.5               175.8
9.0               201.5
9.5               227.1
10.0              252.8
```

**Unidades Suportadas:**
- **Frequência**: Hz, MHz, GHz
- **Campo Magnético**: T, mT, G, Oe (conversão automática)

## Dependências

- **NumPy** ≥ 1.20.0: Cálculos numéricos e operações com arrays
- **SciPy** ≥ 1.7.0: Algoritmos de ajuste de curvas e otimização
- **Matplotlib** ≥ 3.5.0: Visualização e plotagem de dados
- **Pandas** ≥ 1.3.0: Manipulação e tratamento de dados
- **Tkinter**: Framework de GUI (normalmente incluído com Python)

```bash
pip install numpy>=1.20.0 scipy>=1.7.0 matplotlib>=3.5.0 pandas>=1.3.0
```

## Exemplo de Resultados

### Análise de Campos de Ressonância
```
=== EXTRAÇÃO DE CAMPOS DE RESSONÂNCIA FMR ===
Freq(GHz)    Hr(Oe)    Hr(mT)    ΔH(Oe)    ΔH(mT)    R²
5.0          355.8     35.6      17.8      1.8       0.95
6.0          505.8     50.6      21.1      2.1       0.97
7.0          676.3     67.6      24.5      2.4       0.94
8.0          863.5     86.3      28.2      2.8       0.96
9.0          1064.8    106.5     31.9      3.2       0.93
10.0         1280.7    128.1     35.6      3.6       0.95
11.0         1507.6    150.8     39.3      3.9       0.97
12.0         1749.0    174.9     43.0      4.3       0.94
```

### Ajuste de Parâmetros Magnéticos
```
=== RESULTADOS DO AJUSTE DA EQUAÇÃO DE KITTEL ===
Geometria: No Plano (θ = 0°)

Parâmetros Ajustados:
• Ms = 0.845 ± 0.012 T
• Ha = 15.2 ± 2.1 mT
• Hk = 98.7 ± 5.4 mT

Qualidade Estatística:
• R² = 0.998734
• N° pontos de dados = 8
• Faixa de frequência = 5-12 GHz
```

## Aplicações

Este software é projetado para:

- **Análise quantitativa de espectroscopia FMR**
- **Determinação de parâmetros magnéticos** em filmes finos
- **Caracterização de campos de anisotropia** no plano
- **Estudos de dinâmica de magnetização** em materiais ferromagnéticos



## Referências

- **Kittel, C.** "Introduction to Solid State Physics" - Capítulo 15: Ressonância Ferromagnética
- **Farle, M.** "Ferromagnetic resonance of ultrathin metallic layers" Rep. Prog. Phys. **61** 755 (1998)
- **Heinrich, B. & Cochran, J.F.** "Ultrathin metallic magnetic films" Adv. Phys. **42** 523 (1993)

