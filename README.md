# Software de Análise FMR para Filmes Finos

Software Python abrangente para análise de **Ressonância Ferromagnética (FMR)** em filmes magnéticos finos. Este pacote implementa a equação de Kittel simplificada para extrair parâmetros magnéticos (Ms, H_eff) de medidas experimentais de frequência vs. campo de ressonância, incluindo análise automática de largura de linha e ajuste espectral completo com interface gráfica interativa.

## Visão Geral

Este software foi desenvolvido como parte do projeto de pesquisa **PIIC/UFES 2025** "Introdução à Ressonância Magnética Eletrônica" sob supervisão do Prof. Thiago Eduardo Pedreira Bueno. Fornece ferramentas de modelagem teórica e análise computacional para estudar a dinâmica de magnetização em filmes finos ferromagnéticos usando espectroscopia FMR.

## Características

### ✅ Funcionalidades Principais
- **Interface Gráfica Completa**: GUI interativa com processamento automático
- **Extração Automática de Campo de Ressonância**: Ajuste avançado de derivada Lorentziana
- **Ajuste da Equação de Kittel Simplificada**: Modelo de 2 parâmetros para geometria no plano
- **Análise de Largura de Linha**: Extração automática de larguras de linha FMR (ΔH)
- **Processamento Automático**: Upload e análise imediatos sem necessidade de comandos manuais
- **Seleção Múltipla**: Gerenciamento avançado de arquivos (upload/remoção múltiplos)
- **Exportação Profissional**: Múltiplos formatos (PNG, PDF, SVG) com alta resolução
- **Análise Estatística**: Propagação de erros e cálculos de R²

### 🧬 Implementação Física
- **Equação de Kittel Simplificada**: ω₀ = γ√((Hr - H_eff)(Hr - H_eff + Ms))
- **Parâmetros Magnéticos**:
  - **Ms**: Magnetização de saturação (3 casas decimais)
  - **H_eff**: Campo efetivo incluindo anisotropias (1 casa decimal em mT)
- **Geometria**: Configuração no plano otimizada para estabilidade do ajuste
- **Razão Giromagnética**: γ = 2.8×10¹⁰ Hz/T para spins eletrônicos
- **Precisão**: Formatação inteligente baseada em erros estatísticos

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

## Tutorial Visual da Interface Gráfica

A interface gráfica `fmr_gui.py` oferece uma experiência completa para análise de dados FMR através de uma interface intuitiva com múltiplas abas especializadas.

### 1. Interface Principal e Upload de Dados

![Interface Principal](imagens/01interface%20grafica.png)

**Funcionalidades Principais:**
- **Upload Automático**: Botão "Upload Arquivos .dat" carrega e processa múltiplos espectros automaticamente
- **Gerenciamento de Arquivos**: Visualização organizada com seleção/remoção múltipla (Ctrl/Shift + clique)
- **Processamento Automático**: Extração imediata de Hr e ΔH após upload
- **Ajuste Automático de Kittel**: Executado automaticamente quando ≥3 espectros carregados
- **Visualização**:
  - "Mostrar Todos Espectros" - Grade com curvas de ajuste lorentziano
  - Checkboxes para Hr e ajuste individual
- **Exportação Avançada**: Salvamento seletivo de múltiplos gráficos em PNG/PDF/SVG

### 2. Análise Individual de Espectros FMR

![Espectro Individual](imagens/02%20ajuste%20curva%20de%20ressonancia.png)

**Características da Análise Espectral:**
- **Espectro Experimental**: Curva azul mostra dados originais em formato dχ"/dH
- **Campo de Ressonância**: Linha vertical vermelha tracejada indica Hr extraído
- **Ajuste Lorentziano**: Curva vermelha sólida mostra derivada de Lorentziana ajustada
- **Informações Precisas**: Display automático de Hr e largura de linha ΔH
- **Qualidade do Ajuste**: Controle de qualidade com R² para validação estatística

### 3. Ajuste da Equação de Kittel

![Ajuste de Kittel](imagens/04%20ajuste%20kittel.png)

**Implementação da Teoria FMR Simplificada:**
- **Dados Experimentais**: Pontos vermelhos mostram frequência vs. campo de ressonância
- **Ajuste Teórico**: Curva azul implementa equação de Kittel simplificada (2 parâmetros)
- **Parâmetros Magnéticos Extraídos**:
  - **Ms**: Magnetização de saturação (3 casas decimais em T)
  - **H_eff**: Campo efetivo incluindo anisotropias (1 casa decimal em mT)
- **Validação Estatística**: R² > 0.99 típico, demonstrando excelente qualidade do ajuste
- **Faixa Ampla**: Cobertura de 5-12 GHz para caracterização robusta
- **Ajuste Automático**: Executado automaticamente quando ≥3 espectros carregados

### 4. Análise de Largura de Linha vs Frequência

![Largura de Linha](imagens/03%20ajuste%20linear.png)

**Análise de Amortecimento Magnético:**
- **Dados Experimentais**: Pontos vermelhos mostram ΔH extraído de cada espectro
- **Ajuste Linear**: Equação y = 0.227x + 0.353 com R² = 0.9987
- **Interpretação Física**:
  - **Slope (0.227 mT/GHz)**: Relacionado ao amortecimento magnético intrínseco
  - **Intercept (0.353 mT)**: Contribuição de amortecimento inomogêneo
- **Qualidade Excepcional**: R² > 0.99 indica dados de alta qualidade experimental

### 5. Funcionalidade Avançada de Exportação

A interface oferece sistema de exportação inteligente que permite:

**Seleção Múltipla de Gráficos:**
- ✅ Espectro Individual (aba ativa)
- ✅ Largura de Linha vs Frequência (se processado)
- ✅ Ajuste de Kittel (se executado)

**Formatos Profissionais:**
- **PNG**: Alta resolução (300 DPI) para apresentações
- **PDF**: Formato vetorial para publicações científicas
- **SVG**: Editável para customização adicional

**Nomenclatura Organizada:**
- Arquivos salvos como: `nome_base_espectro.png`, `nome_base_kittel.pdf`, etc.
- Validação automática de dados disponíveis
- Feedback detalhado de sucesso/erro

### 6. Fluxo de Trabalho Típico

**Passo a Passo para Análise Completa:**

1. **Preparação dos Dados**:
   ```bash
   # Organizar arquivos .dat no formato correto
   # freq(GHz) corrente(mA) campo(Oe) sinal_fmr(u.a.)
   ```

2. **Carregamento na Interface**:
   - Executar `python fmr_gui.py`
   - Clicar "Upload Arquivos .dat"
   - Selecionar múltiplos arquivos experimentais

3. **Processamento Automático**:
   - **Automático após upload**: Software extrai Hr e ΔH de cada espectro imediatamente
   - **Ajuste de Kittel automático**: Executado quando ≥3 espectros carregados
   - Validação com R² > 0.7 para garantir qualidade

4. **Análise de Parâmetros Magnéticos**:
   - Navegar para aba "Ajuste de Kittel"
   - Visualizar Ms e H_eff extraídos automaticamente
   - Verificar qualidade do ajuste (R² > 0.99 típico)

5. **Análise de Amortecimento**:
   - Aba "Largura de Linha" mostra dependência ΔH vs frequência
   - Ajuste linear automático revela contribuições intrínsecas/extrínseças

6. **Exportação Profissional**:
   - Clicar "Salvar Gráfico"
   - Selecionar plots desejados (Kittel + Largura de Linha)
   - Escolher formato (PDF para publicações)
   - Definir nome base dos arquivos

**Tempo Típico de Análise**: 5-10 minutos para conjunto completo de 8-12 espectros

## Fundamentação Teórica

### Implementação da Equação de Kittel

O software implementa a equação de Kittel **simplificada** para ressonância ferromagnética em filmes finos, otimizada para estabilidade numérica e precisão do ajuste.

#### Modelo Simplificado para Geometria no Plano
```
ω₀ = γ√((Hr - H_eff)(Hr - H_eff + Ms))
```

Resolvendo para o campo de ressonância Hr:
```
Hr = -b + √(b² + 4(ω/γ)²) / 2
onde: b = Ms - 2H_eff
```

**Parâmetros:**
- **γ = 2.8×10¹⁰ Hz/T**: Razão giromagnética para spins eletrônicos
- **Hr**: Campo de ressonância (extraído dos espectros experimentais)
- **Ms**: Magnetização de saturação (parâmetro ajustado, 3 decimais)
- **H_eff**: Campo efetivo incluindo anisotropias (parâmetro ajustado, 1 decimal em mT)

**Vantagens do Modelo Simplificado:**
- Apenas 2 parâmetros → **erros menores** e ajuste mais estável
- Campo efetivo H_eff engloba todas as contribuições de anisotropia
- Convergência robusta com chutes iniciais simples (Ms=1T, H_eff=0)

### Ajuste de Derivada Lorentziana

Para espectros de absorção derivativos (formato dχ"/dH), o software ajusta a forma matematicamente correta:

```
dL/dH = -2A(H-Hr) / [ΔH² × (1 + ((H-Hr)/ΔH)²)²] + offset
```

Esta é a derivada exata da Lorentziana normalizada `L(H) = A / (1 + ((H-Hr)/ΔH)²)`, que permite extração precisa do campo de ressonância (Hr) e largura de linha (ΔH) com controle de qualidade estatístico (limite R² > 0.7).


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


