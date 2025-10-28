#!/usr/bin/env python3
"""
Interface Gráfica para Análise FMR
Permite upload de arquivos .dat e visualização dos espectros
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import os
from typing import List
import threading

# Importações científicas necessárias
import pandas as pd
from lmfit import Model, Parameters
import warnings

class FMRSpectrumAnalyzer:
    """Classe para analisar espectros FMR e extrair campos de ressonância."""

    def __init__(self):
        self.spectra_data = {}
        self.resonance_fields = {}

    def load_spectrum(self, filepath: str):
        """Carrega um espectro FMR de arquivo .dat"""
        try:
            df = pd.read_csv(filepath, sep=r'\s+', engine='python')
            df.columns = ['frequencia', 'corrente', 'campo_magnetico', 'sinal_fmr']
            frequency = df['frequencia'].iloc[0]
            campo = df['campo_magnetico'].values
            sinal = df['sinal_fmr'].values
            return campo, sinal, frequency
        except Exception as e:
            raise RuntimeError(f"Erro ao carregar espectro {filepath}: {e}")

    def _fit_derivative_lorentzian(self, campo: np.ndarray, sinal: np.ndarray):
        """Ajusta derivada de Lorentziana para extrair Hr com precisão usando lmfit."""
        def derivative_lorentzian(H, Hr, A, Delta_H, offset):
            # Forma correta da derivada da Lorentziana: dL/dH = -2A(H-Hr) / (ΔH × (1 + ((H-Hr)/ΔH)²)²)
            normalized_diff = (H - Hr) / Delta_H
            denominator = Delta_H * (1 + normalized_diff**2)**2
            return -2 * A * (H - Hr) / denominator + offset

        min_idx = np.argmin(sinal)
        max_idx = np.argmax(sinal)
        Hr_guess = (campo[min_idx] + campo[max_idx]) / 2
        A_guess = np.max(np.abs(sinal)) * 2
        Delta_H_guess = abs(campo[max_idx] - campo[min_idx]) / 2
        if Delta_H_guess < (campo.max() - campo.min()) * 0.02:
            Delta_H_guess = (campo.max() - campo.min()) * 0.05
        Delta_H_guess = max(Delta_H_guess, (campo.max() - campo.min()) * 0.01)
        Delta_H_guess = min(Delta_H_guess, (campo.max() - campo.min()) * 0.3)
        offset_guess = np.mean([sinal[0], sinal[-1]])

        try:
            # Criar modelo lmfit
            model = Model(derivative_lorentzian, independent_vars=['H'])

            # Definir parâmetros com bounds
            params = Parameters()
            params.add('Hr', value=Hr_guess, min=campo.min(), max=campo.max())
            params.add('A', value=A_guess,
                      min=0.1 * np.max(np.abs(sinal)),
                      max=10 * np.max(np.abs(sinal)))
            params.add('Delta_H', value=Delta_H_guess,
                      min=(campo.max() - campo.min()) * 0.005,
                      max=(campo.max() - campo.min()) * 0.3)
            params.add('offset', value=offset_guess,
                      min=-np.max(np.abs(sinal)),
                      max=np.max(np.abs(sinal)))

            # Fazer ajuste com lmfit
            result = model.fit(sinal, params=params, H=campo, method='leastsq')

            fitted_curve = result.best_fit
            r_squared = 1 - result.residual.var() / np.var(sinal)
            Hr_fitted = result.params['Hr'].value
            Delta_H_fitted = abs(result.params['Delta_H'].value)

            Hr_range = (campo.min(), campo.max())
            valid_fit = (r_squared > 0.3 and Hr_range[0] <= Hr_fitted <= Hr_range[1] and
                        Delta_H_fitted > (campo.max() - campo.min()) * 0.01 and
                        Delta_H_fitted < (campo.max() - campo.min()) * 0.5)

            if valid_fit:
                return Hr_fitted, Delta_H_fitted
            else:
                Hr_backup = self.find_resonance_field_derivative(campo, sinal, "zero_crossing")
                Delta_H_backup = (campo.max() - campo.min()) * 0.05
                return Hr_backup, Delta_H_backup

        except Exception:
            Hr_backup = self.find_resonance_field_derivative(campo, sinal, "zero_crossing")
            Delta_H_backup = (campo.max() - campo.min()) * 0.05
            return Hr_backup, Delta_H_backup

    def find_resonance_field_derivative(self, campo: np.ndarray, sinal: np.ndarray, method: str = "zero_crossing"):
        """Encontra o campo de ressonância em espectro derivativo."""
        if method == "zero_crossing":
            zero_crossings = np.where(np.diff(np.signbit(sinal)))[0]
            if len(zero_crossings) == 0:
                min_idx = np.argmin(sinal)
                return campo[min_idx]
            center_idx = len(campo) // 2
            distances = np.abs(zero_crossings - center_idx)
            closest_crossing = zero_crossings[np.argmin(distances)]
            if closest_crossing < len(campo) - 1:
                x1, x2 = campo[closest_crossing], campo[closest_crossing + 1]
                y1, y2 = sinal[closest_crossing], sinal[closest_crossing + 1]
                if y2 != y1:
                    hr = x1 - y1 * (x2 - x1) / (y2 - y1)
                    return hr
            return campo[closest_crossing]
        elif method == "peak":
            min_idx = np.argmin(sinal)
            return campo[min_idx]
        else:
            raise ValueError(f"Método não reconhecido: {method}")

    def _fit_and_plot_lorentzian(self, campo: np.ndarray, sinal: np.ndarray):
        """Ajusta derivada Lorentziana e retorna Hr e curva ajustada para plotagem usando lmfit."""
        def derivative_lorentzian(H, Hr, A, Delta_H, offset):
            # Forma correta da derivada da Lorentziana: dL/dH = -2A(H-Hr) / (ΔH × (1 + ((H-Hr)/ΔH)²)²)
            normalized_diff = (H - Hr) / Delta_H
            denominator = Delta_H * (1 + normalized_diff**2)**2
            return -2 * A * (H - Hr) / denominator + offset

        min_idx = np.argmin(sinal)
        max_idx = np.argmax(sinal)
        Hr_guess = (campo[min_idx] + campo[max_idx]) / 2
        A_guess = np.max(np.abs(sinal)) * 2
        Delta_H_guess = abs(campo[max_idx] - campo[min_idx]) / 2
        if Delta_H_guess < (campo.max() - campo.min()) * 0.02:
            Delta_H_guess = (campo.max() - campo.min()) * 0.05
        Delta_H_guess = max(Delta_H_guess, (campo.max() - campo.min()) * 0.01)
        Delta_H_guess = min(Delta_H_guess, (campo.max() - campo.min()) * 0.3)
        offset_guess = np.mean([sinal[0], sinal[-1]])

        try:
            # Criar modelo lmfit
            model = Model(derivative_lorentzian, independent_vars=['H'])

            # Definir parâmetros com bounds
            params = Parameters()
            params.add('Hr', value=Hr_guess, min=campo.min(), max=campo.max())
            params.add('A', value=A_guess,
                      min=0.1 * np.max(np.abs(sinal)),
                      max=10 * np.max(np.abs(sinal)))
            params.add('Delta_H', value=Delta_H_guess,
                      min=(campo.max() - campo.min()) * 0.005,
                      max=(campo.max() - campo.min()) * 0.3)
            params.add('offset', value=offset_guess,
                      min=-np.max(np.abs(sinal)),
                      max=np.max(np.abs(sinal)))

            # Fazer ajuste com lmfit
            result = model.fit(sinal, params=params, H=campo, method='leastsq')

            fitted_curve = result.best_fit
            Hr_fitted = result.params['Hr'].value
            Delta_H_fitted = abs(result.params['Delta_H'].value)
            return Hr_fitted, Delta_H_fitted, fitted_curve

        except:
            # Fallback: tentar sem bounds
            try:
                model = Model(derivative_lorentzian, independent_vars=['H'])
                params = Parameters()
                params.add('Hr', value=Hr_guess)
                params.add('A', value=A_guess)
                params.add('Delta_H', value=Delta_H_guess)
                params.add('offset', value=offset_guess)

                result = model.fit(sinal, params=params, H=campo, method='leastsq')

                fitted_curve = result.best_fit
                Hr_fitted = result.params['Hr'].value
                Delta_H_fitted = abs(result.params['Delta_H'].value)
                return Hr_fitted, Delta_H_fitted, fitted_curve
            except:
                # Se falhar completamente, retornar estimativas
                return Hr_guess, Delta_H_guess, sinal


class FMRFitting:
    """Classe para ajuste de dados de ressonância ferromagnética em filmes finos com lmfit."""

    def __init__(self):
        self.gamma = 2.8e10  # razão giromagnética (Hz/T)
        self.mu0 = 4 * np.pi * 1e-7  # permeabilidade magnética do vácuo (H/m)
        # Parâmetros da amostra (podem ser configurados)
        self.espessura = 30e-9  # m (padrão 30 nm)
        self.area = 25e-12  # m² (padrão 5mm x 5mm)

    def kittel_inplane(self, Hr, Ms, H_eff):
        """Equação de Kittel simplificada para geometria no plano.
        f = (γ/2π)√((Hr - H_eff)(Hr - H_eff + Ms))

        onde γ/2π = 2.8×10¹⁰ Hz/T

        Retorna frequência em Hz. Ajustamos frequências experimentais (Hz)
        contra esta função.
        """
        term = (Hr - H_eff) * (Hr - H_eff + Ms)
        # Proteger contra raiz quadrada de número negativo
        term = np.maximum(term, 0)
        return self.gamma * np.sqrt(term)

    def fit_inplane_data(self, frequency_data: np.ndarray, field_data: np.ndarray):
        """Ajusta dados experimentais usando lmfit com PROPAGAÇÃO AUTOMÁTICA de incertezas.

        Parâmetros:
        - frequency_data: array de frequências em Hz
        - field_data: array de campos de ressonância em T

        Retorna:
        - params_dict: dicionário com parâmetros ajustados e propagados
        - fitted_fields: array de campos ajustados em T
        """

        # Criar modelo lmfit
        model = Model(self.kittel_inplane, independent_vars=['Hr'])

        # Definir parâmetros básicos (ajustáveis)
        params = Parameters()
        params.add('Ms', value=0.8, min=0.1, max=2.0)  # Magnetização de saturação (T)
        params.add('H_eff', value=0.01, min=-0.1, max=0.5)  # Campo efetivo (T)

        # Parâmetros fixos da amostra
        params.add('espessura', value=self.espessura, vary=False)
        params.add('area', value=self.area, vary=False)
        params.add('mu0', value=self.mu0, vary=False)

        # ===================================================================
        # PARÂMETROS DERIVADOS (Propagação AUTOMÁTICA de incertezas!)
        # ===================================================================

        # 1. Magnetização total do filme (A·m)
        params.add('Ms_total', expr='Ms * espessura')

        # 2. Campo de desmagnetização (T)
        params.add('H_demag', expr='Ms')  # Para filme fino no plano

        # 3. Campo total efetivo (T)
        params.add('H_total', expr='H_eff + H_demag')

        # 4. Razão Ms/H_eff (adimensional)
        params.add('razao_Ms_Heff', expr='Ms / H_eff')

        # 5. Constante de anisotropia efetiva (J/m³)
        params.add('K_eff', expr='H_eff * Ms / 2')

        # 6. Energia de anisotropia por área (J/m²)
        params.add('K_superficie', expr='K_eff * espessura')

        # 7. Volume do filme (m³)
        params.add('volume', expr='area * espessura')

        # 8. Momento magnético total (A·m²)
        params.add('momento_total', expr='Ms * volume / mu0')

        try:
            # Fazer ajuste com lmfit
            # Igual ao teste_propagacao: ajustamos frequências (Hz) contra função que retorna ω (rad/s)
            result = model.fit(frequency_data, params=params, Hr=field_data,
                             method='leastsq')

            # Extrair parâmetros ajustados com incertezas
            Ms_val = result.params['Ms'].value
            Ms_err = result.params['Ms'].stderr if result.params['Ms'].stderr else 0
            H_eff_val = result.params['H_eff'].value
            H_eff_err = result.params['H_eff'].stderr if result.params['H_eff'].stderr else 0

            # Extrair parâmetros derivados (com incertezas propagadas!)
            H_total_val = result.params['H_total'].value
            H_total_err = result.params['H_total'].stderr if result.params['H_total'].stderr else 0

            K_eff_val = result.params['K_eff'].value
            K_eff_err = result.params['K_eff'].stderr if result.params['K_eff'].stderr else 0

            razao_val = result.params['razao_Ms_Heff'].value
            razao_err = result.params['razao_Ms_Heff'].stderr if result.params['razao_Ms_Heff'].stderr else 0

            # Valores ajustados (são frequências em Hz do modelo)
            fitted_freqs = result.best_fit

            # A GUI espera CAMPOS ajustados
            # Converter de volta: para cada frequência input, calcular o campo de ressonância
            # Equação inversa de Kittel: f = γ_Hz × √((Hr - H_eff)(Hr - H_eff + Ms))
            # Resolver: (Hr - H_eff)(Hr - H_eff + Ms) = (f/γ_Hz)²
            # Expandindo: Hr² + Hr(Ms - 2H_eff) + (H_eff² - H_eff*Ms - (f/γ_Hz)²) = 0
            fitted_fields = np.zeros_like(field_data)
            for i, freq in enumerate(frequency_data):
                # Coeficientes da equação quadrática ax² + bx + c = 0
                a = 1
                b = Ms_val - 2*H_eff_val
                c = H_eff_val**2 - H_eff_val*Ms_val - (freq / self.gamma)**2
                discriminant = b**2 - 4*a*c
                if discriminant >= 0:
                    fitted_fields[i] = (-b + np.sqrt(discriminant)) / (2*a)
                else:
                    fitted_fields[i] = field_data[i]  # Fallback

            # Retornar dicionário com TODOS os parâmetros
            params_dict = {
                # Parâmetros básicos ajustados
                'Ms': Ms_val,
                'Ms_err': Ms_err,
                'H_eff': H_eff_val,
                'H_eff_err': H_eff_err,

                # Parâmetros derivados (propagação automática!)
                'H_total': H_total_val,
                'H_total_err': H_total_err,
                'K_eff': K_eff_val,
                'K_eff_err': K_eff_err,
                'razao_Ms_Heff': razao_val,
                'razao_Ms_Heff_err': razao_err,

                # Objeto result completo para análise avançada
                'lmfit_result': result
            }

            return params_dict, fitted_fields

        except Exception as e:
            raise RuntimeError(f"Erro no ajuste de Kittel: {e}")

    def calculate_r_squared(self, y_data: np.ndarray, y_fit: np.ndarray):
        """Calcula o coeficiente de determinação R²."""
        ss_res = np.sum((y_data - y_fit) ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        return r_squared


class FMRGuiApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FMR Spectrum Analyzer - Interface Gráfica")
        self.root.geometry("1400x800")
        
        # Inicializar analisadores
        self.analyzer = FMRSpectrumAnalyzer()
        self.fmr_fitting = FMRFitting()
        self.loaded_files = []
        self.current_spectrum_data = {}
        self.kittel_params = None
        
        # Criar interface
        self.create_widgets()
        
    def create_widgets(self):
        """Cria todos os widgets da interface"""
        
        # Frame principal
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Frame para controles (lado esquerdo)
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        # Frame para gráficos (lado direito)
        plot_frame = ttk.Frame(main_frame)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # === CONTROLES ===
        
        # Título
        title_label = ttk.Label(control_frame, text="FMR Spectrum Analyzer", 
                               font=('Arial', 14, 'bold'))
        title_label.pack(pady=(0, 20))
        
        # Botão para upload de arquivos
        upload_btn = ttk.Button(control_frame, text="Upload Arquivos .dat",
                               command=self.upload_files)
        upload_btn.pack(fill=tk.X, pady=5)

        # Botão para remover arquivos selecionados
        remove_btn = ttk.Button(control_frame, text="Remover Arquivo(s) Selecionado(s)",
                               command=self.remove_selected_file)
        remove_btn.pack(fill=tk.X, pady=5)

        # Lista de arquivos carregados
        files_label = ttk.Label(control_frame, text="Arquivos Carregados:")
        files_label.pack(anchor=tk.W, pady=(20, 5))
        
        # Frame para listbox e scrollbar
        listbox_frame = ttk.Frame(control_frame)
        listbox_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Scrollbar
        scrollbar = ttk.Scrollbar(listbox_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Listbox de arquivos (com seleção múltipla)
        self.files_listbox = tk.Listbox(listbox_frame, yscrollcommand=scrollbar.set,
                                       height=8, width=30, selectmode=tk.EXTENDED)
        self.files_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.files_listbox.yview)
        
        # Bind para seleção
        self.files_listbox.bind('<<ListboxSelect>>', self.on_file_select)

        # Botão para mostrar todos os espectros
        show_all_btn = ttk.Button(control_frame, text="Mostrar Todos Espectros",
                                 command=self.show_all_spectra)
        show_all_btn.pack(fill=tk.X, pady=(20, 10))

        # Opções de visualização
        viz_frame = ttk.LabelFrame(control_frame, text="Visualização", padding=10)
        viz_frame.pack(fill=tk.X, pady=10)
        
        self.show_resonance_var = tk.BooleanVar(value=True)
        resonance_check = ttk.Checkbutton(viz_frame, text="Mostrar Hr", 
                                         variable=self.show_resonance_var,
                                         command=self.update_plot)
        resonance_check.pack(anchor=tk.W)
        
        self.show_fit_var = tk.BooleanVar(value=False)
        fit_check = ttk.Checkbutton(viz_frame, text="Mostrar Ajuste", 
                                   variable=self.show_fit_var,
                                   command=self.update_plot)
        fit_check.pack(anchor=tk.W)
        
        # Botões de exportação
        export_frame = ttk.LabelFrame(control_frame, text="Exportar", padding=10)
        export_frame.pack(fill=tk.X, pady=10)

        save_plot_btn = ttk.Button(export_frame, text="Salvar Gráfico",
                                  command=self.save_plot)
        save_plot_btn.pack(fill=tk.X, pady=2)
        
        # === ÁREA DE GRÁFICOS ===
        
        # Notebook para abas
        self.notebook = ttk.Notebook(plot_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Aba para espectro individual
        self.spectrum_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.spectrum_frame, text="Espectro Individual")

        # Aba para largura de linha
        self.linewidth_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.linewidth_frame, text="Largura de Linha")

        # Aba para ajuste de Kittel
        self.kittel_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.kittel_frame, text="Ajuste de Kittel")

        # Criar gráficos
        self.create_spectrum_plot()
        self.create_linewidth_plot()
        self.create_kittel_plot()
        
        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Pronto para carregar arquivos")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, 
                              relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
    
    def create_spectrum_plot(self):
        """Cria o gráfico para espectro individual"""
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.fig, self.spectrum_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Toolbar para navegação
        toolbar = tk.Frame(self.spectrum_frame)
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)

    def create_linewidth_plot(self):
        """Cria o gráfico para largura de linha vs frequência"""
        self.linewidth_fig = Figure(figsize=(8, 6), dpi=100)
        self.linewidth_ax = self.linewidth_fig.add_subplot(111)

        self.linewidth_canvas = FigureCanvasTkAgg(self.linewidth_fig, self.linewidth_frame)
        self.linewidth_canvas.draw()
        self.linewidth_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def create_kittel_plot(self):
        """Cria o gráfico para ajuste de Kittel"""
        self.kittel_fig = Figure(figsize=(8, 6), dpi=100)
        self.kittel_ax = self.kittel_fig.add_subplot(111)

        self.kittel_canvas = FigureCanvasTkAgg(self.kittel_fig, self.kittel_frame)
        self.kittel_canvas.draw()
        self.kittel_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def upload_files(self):
        """Upload de arquivos individuais"""
        filetypes = [('Arquivos DAT', '*.dat'), ('Todos os arquivos', '*.*')]
        files = filedialog.askopenfilenames(title="Selecionar arquivos .dat", 
                                           filetypes=filetypes)
        
        if files:
            self.load_files(files)
    
    def load_files(self, filepaths: List[str]):
        """Carrega arquivos na interface"""
        loaded_count = 0
        errors = []
        total_files = len(filepaths)

        for i, filepath in enumerate(filepaths, 1):
            # Atualizar status de progresso
            filename = os.path.basename(filepath)
            self.status_var.set(f"Carregando {i}/{total_files}: {filename}")
            self.root.update()  # Forçar atualização da interface

            if filepath not in [f['path'] for f in self.loaded_files]:
                try:
                    # Verificar se arquivo existe
                    if not os.path.exists(filepath):
                        raise FileNotFoundError(f"Arquivo não encontrado")

                    # Carregar espectro básico para obter informações
                    campo, sinal, freq = self.analyzer.load_spectrum(filepath)

                    file_info = {
                        'path': filepath,
                        'filename': os.path.basename(filepath),
                        'frequency': freq,
                        'campo': campo,
                        'sinal': sinal,
                        'processed': False
                    }

                    self.loaded_files.append(file_info)

                    # Adicionar à listbox
                    display_name = f"{freq:.1f} GHz - {os.path.basename(filepath)}"
                    self.files_listbox.insert(tk.END, display_name)
                    loaded_count += 1

                except Exception as e:
                    error_msg = f"{os.path.basename(filepath)}: {str(e)}"
                    errors.append(error_msg)

        # Processar automaticamente os novos arquivos
        if loaded_count > 0:
            self.status_var.set(f"Processando {loaded_count} arquivo(s)...")
            self.root.update()
            self.process_new_files()

        # Atualizar status final
        if loaded_count > 0:
            self.status_var.set(f"✓ {len(self.loaded_files)} arquivo(s) carregado(s) e processado(s)")
        else:
            self.status_var.set(f"{len(self.loaded_files)} arquivo(s) carregado(s)")

        # Mostrar erros se houver
        if errors:
            error_summary = f"Problemas ao carregar {len(errors)} arquivo(s):\n\n" + "\n".join(errors[:5])
            if len(errors) > 5:
                error_summary += f"\n... e mais {len(errors) - 5} erro(s)"
            messagebox.showwarning("Avisos de Carregamento", error_summary)

    
    def on_file_select(self, event=None):
        """Callback para seleção de arquivo na lista"""
        selection = self.files_listbox.curselection()
        if selection:
            index = selection[0]
            file_info = self.loaded_files[index]
            self.plot_spectrum(file_info)

    def remove_selected_file(self):
        """Remove arquivos selecionados da lista"""
        selection = self.files_listbox.curselection()
        if not selection:
            messagebox.showwarning("Aviso", "Selecione pelo menos um arquivo para remover")
            return

        # Preparar mensagem de confirmação
        if len(selection) == 1:
            file_info = self.loaded_files[selection[0]]
            confirm_msg = f"Deseja remover o arquivo:\n{file_info['filename']}?"
        else:
            confirm_msg = f"Deseja remover {len(selection)} arquivos selecionados?"

        # Confirmar remoção
        confirm = messagebox.askyesno("Confirmar Remoção", confirm_msg)
        if not confirm:
            return

        # Coletar frequências removidas e índices (ordenar em ordem reversa para remoção)
        removed_freqs = []
        removed_names = []
        indices_to_remove = sorted(selection, reverse=True)

        # Remover arquivos (de trás para frente para não afetar índices)
        for index in indices_to_remove:
            file_info = self.loaded_files[index]
            removed_freq = file_info['frequency']
            removed_freqs.append(removed_freq)
            removed_names.append(file_info['filename'])

            # Remover da lista de arquivos carregados
            self.loaded_files.pop(index)

            # Remover da listbox
            self.files_listbox.delete(index)

            # Remover dos dados do analyzer (se processado)
            if removed_freq in self.analyzer.spectra_data:
                del self.analyzer.spectra_data[removed_freq]
            if removed_freq in self.analyzer.resonance_fields:
                del self.analyzer.resonance_fields[removed_freq]

        # Limpar gráfico
        self.ax.clear()
        self.ax.text(0.5, 0.5, 'Selecione um arquivo para visualizar',
                    ha='center', va='center', transform=self.ax.transAxes)
        self.canvas.draw()

        # Atualizar plots de largura de linha e Kittel
        self.update_linewidth_plot()
        if self.kittel_params is not None:
            # Limpar ajuste de Kittel (necessário refazer)
            self.kittel_params = None
            self.kittel_ax.clear()
            self.kittel_ax.text(0.5, 0.5, 'Execute o ajuste de Kittel novamente',
                               ha='center', va='center', transform=self.kittel_ax.transAxes)
            self.kittel_canvas.draw()

        # Atualizar status
        if len(removed_names) == 1:
            self.status_var.set(f"Arquivo removido: {removed_names[0]}")
        else:
            self.status_var.set(f"{len(removed_names)} arquivos removidos")

    
    def plot_spectrum(self, file_info: dict):
        """Plota espectro individual"""
        self.ax.clear()
        
        # Plot básico
        self.ax.plot(file_info['campo'], file_info['sinal'], 'b-', 
                    linewidth=1.5, label='Espectro FMR experimental')
        
        # Mostrar campo de ressonância se processado
        if file_info['processed'] and self.show_resonance_var.get():
            if 'Hr' in file_info:
                hr_label = f'Hr = {file_info["Hr"]:.1f} Oe'
                if 'Delta_H' in file_info:
                    hr_label += f', ΔH = {file_info["Delta_H"]:.1f} Oe'
                self.ax.axvline(file_info['Hr'], color='green', linestyle='--', 
                               linewidth=2, label=hr_label)
        
        # Mostrar ajuste se solicitado
        if file_info['processed'] and self.show_fit_var.get():
            if 'fitted_curve' in file_info:
                fit_label = f'Ajuste Lorentziano'
                self.ax.plot(file_info['campo'], file_info['fitted_curve'], 'r-', 
                            linewidth=2, alpha=0.7, label=fit_label)
        
        self.ax.set_xlabel('Campo Magnético (Oe)')
        self.ax.set_ylabel('dP"/dH (u.a.)')
        self.ax.set_title(f'Espectro FMR - {file_info["frequency"]:.1f} GHz')
        self.ax.grid(True, alpha=0.3)
        self.ax.legend()
        
        self.fig.tight_layout()
        self.canvas.draw()
    
    def update_plot(self):
        """Atualiza o plot atual"""
        selection = self.files_listbox.curselection()
        if selection:
            index = selection[0]
            file_info = self.loaded_files[index]
            self.plot_spectrum(file_info)
    
    def process_new_files(self):
        """Processa apenas os arquivos não processados"""
        unprocessed_files = [f for f in self.loaded_files if not f['processed']]

        if not unprocessed_files:
            return

        for file_info in unprocessed_files:
            try:
                # Extrair Hr e Delta_H
                Hr, Delta_H = self.analyzer._fit_derivative_lorentzian(
                    file_info['campo'], file_info['sinal'])

                # Gerar curva ajustada para visualização
                try:
                    _, _, fitted_curve = self.analyzer._fit_and_plot_lorentzian(
                        file_info['campo'], file_info['sinal'])
                    file_info['fitted_curve'] = fitted_curve
                except:
                    # Se falhar, não salvar curva ajustada
                    pass

                # Atualizar informações do arquivo
                file_info['Hr'] = Hr
                file_info['Delta_H'] = Delta_H
                file_info['processed'] = True

                # Armazenar no analyzer para compatibilidade
                freq = file_info['frequency']
                self.analyzer.spectra_data[freq] = {
                    'campo': file_info['campo'],
                    'sinal': file_info['sinal'],
                    'Hr': Hr,
                    'Delta_H': Delta_H
                }
                self.analyzer.resonance_fields[freq] = Hr

            except Exception as e:
                print(f"Erro ao processar {file_info['filename']}: {e}")

        # Atualizar plots
        self.update_linewidth_plot()

        # Tentar ajuste de Kittel automaticamente se houver dados suficientes
        self.auto_fit_kittel()

    def process_all_files(self):
        """Processa todos os arquivos carregados (mantido para compatibilidade)"""
        if not self.loaded_files:
            messagebox.showwarning("Aviso", "Nenhum arquivo carregado")
            return

        self.status_var.set("Processando arquivos...")
        self.process_new_files()
        self.status_var.set(f"Processamento concluído: {len(self.loaded_files)} arquivos")
    
    def show_all_spectra(self):
        """Mostra todos os espectros processados em janela separada"""
        if not any(f['processed'] for f in self.loaded_files):
            messagebox.showwarning("Aviso", "Nenhum arquivo foi processado ainda")
            return
        
        # Criar nova janela
        all_window = tk.Toplevel(self.root)
        all_window.title("Todos os Espectros")
        all_window.geometry("1200x800")
        
        # Criar figura com subplots
        fig = Figure(figsize=(12, 8), dpi=100)
        
        processed_files = [f for f in self.loaded_files if f['processed']]
        n_plots = len(processed_files)
        
        if n_plots > 8:
            n_plots = 8
            processed_files = processed_files[:8]
        
        # Calcular layout de subplots
        cols = 4
        rows = (n_plots + cols - 1) // cols
        
        for i, file_info in enumerate(processed_files):
            ax = fig.add_subplot(rows, cols, i + 1)

            # Plot do espectro experimental
            ax.plot(file_info['campo'], file_info['sinal'], 'b-', linewidth=1, label='Dados')

            # Plot da curva ajustada (se disponível)
            if 'fitted_curve' in file_info:
                ax.plot(file_info['campo'], file_info['fitted_curve'], 'r-',
                       linewidth=1.5, alpha=0.7, label='Ajuste')

            # Linha vertical no Hr
            if 'Hr' in file_info:
                ax.axvline(file_info['Hr'], color='green', linestyle='--', linewidth=1)

            title = f'{file_info["frequency"]:.1f} GHz'
            if 'Delta_H' in file_info:
                title += f'\nΔH = {file_info["Delta_H"]:.1f} Oe'

            ax.set_title(title, fontsize=10)
            ax.set_xlabel('Campo (Oe)', fontsize=8)
            ax.set_ylabel('dχ"/dH', fontsize=8)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=8)
            ax.legend(fontsize=6, loc='best')
        
        fig.tight_layout()
        
        # Canvas
        canvas = FigureCanvasTkAgg(fig, all_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def save_plot(self):
        """Salva gráfico(s) com opções de seleção"""
        # Criar janela de seleção
        save_window = tk.Toplevel(self.root)
        save_window.title("Salvar Gráficos")
        save_window.geometry("400x300")
        save_window.transient(self.root)
        save_window.grab_set()

        # Centralizar janela
        save_window.update_idletasks()
        x = (save_window.winfo_screenwidth() // 2) - (400 // 2)
        y = (save_window.winfo_screenheight() // 2) - (300 // 2)
        save_window.geometry(f"400x300+{x}+{y}")

        # Frame principal
        main_frame = ttk.Frame(save_window, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Título
        title_label = ttk.Label(main_frame, text="Selecione os gráficos para salvar:",
                               font=('Arial', 12, 'bold'))
        title_label.pack(pady=(0, 20))

        # Variáveis de controle para checkboxes
        save_spectrum_var = tk.BooleanVar(value=True)
        save_linewidth_var = tk.BooleanVar(value=False)
        save_kittel_var = tk.BooleanVar(value=False)

        # Checkboxes
        spectrum_check = ttk.Checkbutton(main_frame, text="Espectro Individual",
                                        variable=save_spectrum_var)
        spectrum_check.pack(anchor=tk.W, pady=5)

        linewidth_check = ttk.Checkbutton(main_frame, text="Largura de Linha vs Frequência",
                                         variable=save_linewidth_var)
        linewidth_check.pack(anchor=tk.W, pady=5)

        kittel_check = ttk.Checkbutton(main_frame, text="Ajuste de Kittel",
                                      variable=save_kittel_var)
        kittel_check.pack(anchor=tk.W, pady=5)

        # Verificar se os dados necessários estão disponíveis
        if not any(f['processed'] for f in self.loaded_files):
            linewidth_check.config(state='disabled')

        if self.kittel_params is None:
            kittel_check.config(state='disabled')
        else:
            save_kittel_var.set(True)  # Marcar por padrão se disponível

        # Se há dados de largura de linha, marcar por padrão
        if any(f.get('processed') and 'Delta_H' in f for f in self.loaded_files):
            save_linewidth_var.set(True)

        # Frame para formato de arquivo
        format_frame = ttk.LabelFrame(main_frame, text="Formato de arquivo", padding="10")
        format_frame.pack(fill=tk.X, pady=(20, 10))

        format_var = tk.StringVar(value="png")
        ttk.Radiobutton(format_frame, text="PNG (recomendado)", variable=format_var, value="png").pack(anchor=tk.W)
        ttk.Radiobutton(format_frame, text="PDF (vetorial)", variable=format_var, value="pdf").pack(anchor=tk.W)
        ttk.Radiobutton(format_frame, text="SVG (vetorial)", variable=format_var, value="svg").pack(anchor=tk.W)

        # Frame para botões
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=(20, 0))

        def save_selected():
            """Salva os gráficos selecionados"""
            plots_to_save = []

            if save_spectrum_var.get():
                plots_to_save.append(("espectro", "Espectro Individual", self.fig))
            if save_linewidth_var.get():
                plots_to_save.append(("largura_linha", "Largura de Linha", self.linewidth_fig))
            if save_kittel_var.get():
                plots_to_save.append(("kittel", "Ajuste de Kittel", self.kittel_fig))

            if not plots_to_save:
                messagebox.showwarning("Aviso", "Selecione pelo menos um gráfico para salvar")
                return

            # Escolher diretório base
            directory = filedialog.askdirectory(title="Escolher pasta para salvar os gráficos")
            if not directory:
                return

            # Escolher nome base para os arquivos
            base_name = simpledialog.askstring("Nome dos arquivos",
                                              "Nome base para os arquivos:",
                                              initialvalue="fmr_analise")
            if not base_name:
                base_name = "fmr_analise"

            # Salvar cada gráfico
            saved_files = []
            extension = f".{format_var.get()}"

            for plot_id, plot_name, figure in plots_to_save:
                try:
                    filename = os.path.join(directory, f"{base_name}_{plot_id}{extension}")
                    figure.savefig(filename, dpi=300, bbox_inches='tight')
                    saved_files.append(f"• {plot_name}: {os.path.basename(filename)}")
                except Exception as e:
                    messagebox.showerror("Erro", f"Erro ao salvar {plot_name}:\n{str(e)}")
                    return

            # Mostrar sucesso
            success_msg = f"Gráficos salvos com sucesso em:\n{directory}\n\n" + "\n".join(saved_files)
            messagebox.showinfo("Sucesso", success_msg)
            save_window.destroy()

        def cancel_save():
            save_window.destroy()

        # Botões
        ttk.Button(button_frame, text="Salvar", command=save_selected).pack(side=tk.RIGHT, padx=(5, 0))
        ttk.Button(button_frame, text="Cancelar", command=cancel_save).pack(side=tk.RIGHT)

    def update_linewidth_plot(self):
        """Atualiza o plot de largura de linha vs frequência"""
        processed_files = [f for f in self.loaded_files if f['processed'] and 'Delta_H' in f]

        if len(processed_files) < 2:
            self.linewidth_ax.clear()
            self.linewidth_ax.text(0.5, 0.5, 'Processamento insuficiente\n(mín. 2 espectros)',
                                  ha='center', va='center', transform=self.linewidth_ax.transAxes)
            self.linewidth_canvas.draw()
            return

        # Extrair dados
        frequencies = []
        linewidths = []

        for file_info in processed_files:
            frequencies.append(file_info['frequency'])
            linewidths.append(file_info['Delta_H'] * 0.1)  # Oe -> mT

        frequencies = np.array(frequencies)
        linewidths = np.array(linewidths)

        # Ordenar por frequência
        sort_idx = np.argsort(frequencies)
        frequencies = frequencies[sort_idx]
        linewidths = linewidths[sort_idx]

        self.linewidth_ax.clear()

        # Plot dos dados
        self.linewidth_ax.scatter(frequencies, linewidths, color='red', s=80, alpha=0.8,
                                 label='Largura de linha experimental', zorder=3)

        # Ajuste linear
        try:
            coeffs = np.polyfit(frequencies, linewidths, 1)
            freq_fit = np.linspace(frequencies.min(), frequencies.max(), 100)
            linewidth_fit = np.polyval(coeffs, freq_fit)

            self.linewidth_ax.plot(freq_fit, linewidth_fit, 'b-', linewidth=3,
                                  label='Ajuste linear', zorder=2)

            # Calcular R²
            linewidth_pred = np.polyval(coeffs, frequencies)
            ss_res = np.sum((linewidths - linewidth_pred) ** 2)
            ss_tot = np.sum((linewidths - np.mean(linewidths)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)

            # Adicionar texto com informações do ajuste
            # Usar mais casas decimais se R² muito próximo de 1
            r2_format = f'{r_squared:.4f}' if r_squared > 0.999 else f'{r_squared:.3f}'
            textstr = f'Ajuste Linear:\ny = {coeffs[0]:.2f}x + {coeffs[1]:.2f}\nR² = {r2_format}'
            props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
            self.linewidth_ax.text(0.05, 0.95, textstr, transform=self.linewidth_ax.transAxes,
                                  fontsize=10, verticalalignment='top', bbox=props)

        except Exception as e:
            print(f"Erro no ajuste linear: {e}")

        self.linewidth_ax.set_xlabel('Frequência (GHz)', fontsize=12)
        self.linewidth_ax.set_ylabel('Largura de Linha ΔH (mT)', fontsize=12)
        self.linewidth_ax.set_title('Largura de Linha vs Frequência', fontsize=14)
        self.linewidth_ax.legend()
        self.linewidth_ax.grid(True, alpha=0.3)

        # Formatar eixo x com uma casa decimal
        from matplotlib.ticker import FormatStrFormatter
        self.linewidth_ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        self.linewidth_fig.tight_layout()
        self.linewidth_canvas.draw()

    def auto_fit_kittel(self):
        """Tenta fazer ajuste de Kittel automaticamente (sem mensagens)"""
        processed_files = [f for f in self.loaded_files if f['processed'] and 'Hr' in f]

        if len(processed_files) < 3:
            return  # Silenciosamente retorna se não há dados suficientes

        try:
            # Preparar dados para ajuste
            frequencies = []
            fields_oe = []

            for file_info in processed_files:
                frequencies.append(file_info['frequency'] * 1e9)  # GHz -> Hz
                fields_oe.append(file_info['Hr'])  # Oe

            frequencies = np.array(frequencies)
            fields_tesla = np.array(fields_oe) * 1e-4  # Oe -> T

            # Ordenar dados por frequência (independente da ordem de upload)
            sort_idx = np.argsort(frequencies)
            frequencies = frequencies[sort_idx]
            fields_tesla = fields_tesla[sort_idx]

            # Fazer ajuste no plano
            self.kittel_params, fitted_fields = self.fmr_fitting.fit_inplane_data(frequencies, fields_tesla)

            # Atualizar plot de Kittel
            self.update_kittel_plot(frequencies, fields_tesla, fitted_fields)

        except Exception as e:
            print(f"Erro no ajuste automático de Kittel: {e}")

    def fit_kittel_data(self):
        """Executa ajuste de Kittel e atualiza plot (mantido para compatibilidade)"""
        processed_files = [f for f in self.loaded_files if f['processed'] and 'Hr' in f]

        if len(processed_files) < 3:
            messagebox.showwarning("Aviso", "Necessário pelo menos 3 espectros processados para ajuste de Kittel")
            return

        try:
            # Preparar dados para ajuste
            frequencies = []
            fields_oe = []

            for file_info in processed_files:
                frequencies.append(file_info['frequency'] * 1e9)  # GHz -> Hz
                fields_oe.append(file_info['Hr'])  # Oe

            frequencies = np.array(frequencies)
            fields_tesla = np.array(fields_oe) * 1e-4  # Oe -> T

            # Ordenar dados por frequência (independente da ordem de upload)
            sort_idx = np.argsort(frequencies)
            frequencies = frequencies[sort_idx]
            fields_tesla = fields_tesla[sort_idx]

            # Fazer ajuste no plano com lmfit
            self.kittel_params, fitted_fields = self.fmr_fitting.fit_inplane_data(frequencies, fields_tesla)

            # Atualizar plot de Kittel
            self.update_kittel_plot(frequencies, fields_tesla, fitted_fields)

            # Mostrar parâmetros com propagação de incertezas
            ms_str = f"{self.kittel_params['Ms']:.3f} ± {self.kittel_params['Ms_err']:.3f} T"
            heff_str = f"{self.kittel_params['H_eff']*1000:.1f} ± {self.kittel_params['H_eff_err']*1000:.1f} mT"

            # Parâmetros derivados (PROPAGADOS automaticamente!)
            htotal_str = f"{self.kittel_params['H_total']*1000:.1f} ± {self.kittel_params['H_total_err']*1000:.1f} mT"
            keff_str = f"{self.kittel_params['K_eff']:.2e} ± {self.kittel_params['K_eff_err']:.2e} J/m³"
            razao_str = f"{self.kittel_params['razao_Ms_Heff']:.1f} ± {self.kittel_params['razao_Ms_Heff_err']:.1f}"

            r2 = self.fmr_fitting.calculate_r_squared(fields_tesla, fitted_fields)

            params_text = f"""Ajuste de Kittel com lmfit - Propagação Automática de Incertezas

PARÂMETROS AJUSTADOS:
  Ms     = {ms_str}
  H_eff  = {heff_str}
  R²     = {r2:.4f}

PARÂMETROS DERIVADOS (propagação automática!):
  H_total      = {htotal_str}
  K_eff        = {keff_str}
  Ms/H_eff     = {razao_str}

O lmfit propagou automaticamente as incertezas usando
a matriz de covariância completa!"""

            messagebox.showinfo("Ajuste de Kittel - lmfit", params_text)

        except Exception as e:
            messagebox.showerror("Erro", f"Erro no ajuste de Kittel:\n{str(e)}")

    def _format_with_error(self, value, error, unit=''):
        """Formata valor ± erro com casas decimais apropriadas (máximo 2 casas)"""
        if error == 0:
            return f"{value:.2f} ± {error:.2f} {unit}"

        # Determinar casas decimais baseado na magnitude do erro
        # Queremos mostrar 1 dígito significativo no erro
        log_error = np.log10(abs(error))

        # Número de casas decimais para ter 1 dígito significativo no erro
        if error >= 1:
            error_decimals = 0  # Inteiro se erro >= 1
        else:
            # Casas decimais = posição do primeiro dígito significativo
            error_decimals = int(-np.floor(log_error))

        # Limitar a no máximo 2 casas decimais
        decimals = min(error_decimals, 2)

        format_str = f"{{:.{decimals}f}}"
        return f"{format_str.format(value)} ± {format_str.format(error)} {unit}"

    def update_kittel_plot(self, frequencies, fields_experimental, fields_fitted):
        """Atualiza o plot do ajuste de Kittel com parâmetros propagados"""
        self.kittel_ax.clear()

        # Converter unidades para plot
        freq_ghz = frequencies / 1e9
        field_mt_exp = fields_experimental * 1000
        field_mt_fit = fields_fitted * 1000

        # Plot dados experimentais e ajuste
        self.kittel_ax.scatter(freq_ghz, field_mt_exp, color='red', s=80, alpha=0.8,
                              label='Dados experimentais', zorder=3)
        self.kittel_ax.plot(freq_ghz, field_mt_fit, 'b-', linewidth=3,
                           label='Ajuste teórico (Kittel - lmfit)', zorder=2)

        self.kittel_ax.set_xlabel('Frequência (GHz)', fontsize=12)
        self.kittel_ax.set_ylabel('Campo de Ressonância (mT)', fontsize=12)
        self.kittel_ax.set_title('Ajuste de Kittel com Propagação de Incertezas (lmfit)', fontsize=14)
        self.kittel_ax.legend()
        self.kittel_ax.grid(True, alpha=0.3)

        # Formatar eixo x com uma casa decimal
        from matplotlib.ticker import FormatStrFormatter
        self.kittel_ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        # Adicionar texto com parâmetros se disponível
        if self.kittel_params:
            r2 = self.fmr_fitting.calculate_r_squared(fields_experimental, fields_fitted)

            # Parâmetros básicos ajustados
            ms_str = f"{self.kittel_params['Ms']:.3f} ± {self.kittel_params['Ms_err']:.3f} T"
            heff_str = f"{self.kittel_params['H_eff']*1000:.1f} ± {self.kittel_params['H_eff_err']*1000:.1f} mT"

            # Parâmetros derivados (PROPAGADOS!)
            htotal_str = f"{self.kittel_params['H_total']*1000:.1f} ± {self.kittel_params['H_total_err']*1000:.1f} mT"
            keff_str = f"{self.kittel_params['K_eff']:.2e} ± {self.kittel_params['K_eff_err']:.2e} J/m³"

            textstr = f"""Parâmetros Ajustados:
Ms = {ms_str}
H_eff = {heff_str}

Propagados (lmfit):
H_total = {htotal_str}
K_eff = {keff_str}

R² = {r2:.4f}"""

            props = dict(boxstyle='round', facecolor='lightblue', alpha=0.9)
            self.kittel_ax.text(0.02, 0.98, textstr, transform=self.kittel_ax.transAxes,
                               fontsize=9, verticalalignment='top', bbox=props)

        self.kittel_fig.tight_layout()
        self.kittel_canvas.draw()

def main():
    """Função principal para executar a GUI"""
    root = tk.Tk()
    FMRGuiApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()