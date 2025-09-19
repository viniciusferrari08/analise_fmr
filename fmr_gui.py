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
from scipy.optimize import curve_fit
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
        """Ajusta derivada de Lorentziana para extrair Hr com precisão."""
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
            Hr_range = (campo.min(), campo.max())
            A_range = (0.1 * np.max(np.abs(sinal)), 10 * np.max(np.abs(sinal)))
            Delta_H_range = ((campo.max() - campo.min()) * 0.005, (campo.max() - campo.min()) * 0.3)
            offset_range = (-np.max(np.abs(sinal)), np.max(np.abs(sinal)))
            bounds = ([Hr_range[0], A_range[0], Delta_H_range[0], offset_range[0]],
                     [Hr_range[1], A_range[1], Delta_H_range[1], offset_range[1]])

            popt, _ = curve_fit(derivative_lorentzian, campo, sinal,
                               p0=[Hr_guess, A_guess, Delta_H_guess, offset_guess],
                               bounds=bounds, maxfev=20000)

            fitted_curve = derivative_lorentzian(campo, *popt)
            r_squared = 1 - np.sum((sinal - fitted_curve)**2) / np.sum((sinal - np.mean(sinal))**2)
            Hr_fitted = popt[0]
            Delta_H_fitted = abs(popt[2])

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
        """Ajusta derivada Lorentziana e retorna Hr e curva ajustada para plotagem."""
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
            Hr_range = (campo.min(), campo.max())
            A_range = (0.1 * np.max(np.abs(sinal)), 10 * np.max(np.abs(sinal)))
            Delta_H_range = ((campo.max() - campo.min()) * 0.005, (campo.max() - campo.min()) * 0.3)
            offset_range = (-np.max(np.abs(sinal)), np.max(np.abs(sinal)))
            bounds = ([Hr_range[0], A_range[0], Delta_H_range[0], offset_range[0]],
                     [Hr_range[1], A_range[1], Delta_H_range[1], offset_range[1]])

            popt, _ = curve_fit(derivative_lorentzian, campo, sinal,
                               p0=[Hr_guess, A_guess, Delta_H_guess, offset_guess],
                               bounds=bounds, maxfev=20000)
        except:
            popt, _ = curve_fit(derivative_lorentzian, campo, sinal,
                               p0=[Hr_guess, A_guess, Delta_H_guess, offset_guess],
                               maxfev=20000)

        fitted_curve = derivative_lorentzian(campo, *popt)
        Hr_fitted = popt[0]
        Delta_H_fitted = abs(popt[2])
        return Hr_fitted, Delta_H_fitted, fitted_curve


class FMRFitting:
    """Classe para ajuste de dados de ressonância ferromagnética em filmes finos."""

    def __init__(self):
        self.gamma = 2.8e10  # razão giromagnética (Hz/T)
        self.mu0 = 4 * np.pi * 1e-7  # permeabilidade magnética do vácuo (H/m)

    def resonance_field_inplane(self, frequency: np.ndarray, Ms: float, Ha: float = 0, Hk: float = 0):
        """Campo de ressonância para configuração no plano (theta = 0°)."""
        omega = 2 * np.pi * frequency
        a = 1
        b = 2 * Ha + Ms - Hk
        c = Ha**2 + Ha * (Ms - Hk) - (omega / self.gamma)**2
        discriminant = b**2 - 4*a*c
        if np.any(discriminant < 0):
            warnings.warn("Discriminante negativo encontrado. Verifique os parâmetros.")
        Hr = (-b + np.sqrt(np.maximum(discriminant, 0))) / (2*a)
        return Hr

    def fit_inplane_data(self, frequency_data: np.ndarray, field_data: np.ndarray):
        """Ajusta dados experimentais para configuração no plano."""
        def model_func(freq, Ms, Ha, Hk):
            return self.resonance_field_inplane(freq, Ms, Ha, Hk)

        p0 = [1.0, 0.01, 0.1]  # Ms, Ha, Hk

        try:
            popt, pcov = curve_fit(model_func, frequency_data, field_data, p0=p0)
            perr = np.sqrt(np.diag(pcov))
            fitted_values = model_func(frequency_data, *popt)

            params = {
                'Ms': popt[0], 'Ha': popt[1], 'Hk': popt[2],
                'Ms_err': perr[0], 'Ha_err': perr[1], 'Hk_err': perr[2]
            }

            return params, fitted_values

        except Exception as e:
            raise RuntimeError(f"Erro no ajuste: {e}")

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
        
        
        # Lista de arquivos carregados
        files_label = ttk.Label(control_frame, text="Arquivos Carregados:")
        files_label.pack(anchor=tk.W, pady=(20, 5))
        
        # Frame para listbox e scrollbar
        listbox_frame = ttk.Frame(control_frame)
        listbox_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Scrollbar
        scrollbar = ttk.Scrollbar(listbox_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Listbox de arquivos
        self.files_listbox = tk.Listbox(listbox_frame, yscrollcommand=scrollbar.set,
                                       height=8, width=30)
        self.files_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.files_listbox.yview)
        
        # Bind para seleção
        self.files_listbox.bind('<<ListboxSelect>>', self.on_file_select)
        
        # Botões de análise
        analysis_frame = ttk.LabelFrame(control_frame, text="Análise", padding=10)
        analysis_frame.pack(fill=tk.X, pady=(20, 10))
        
        process_btn = ttk.Button(analysis_frame, text="Processar Todos", 
                                command=self.process_all_files)
        process_btn.pack(fill=tk.X, pady=2)
        
        show_all_btn = ttk.Button(analysis_frame, text="Mostrar Todos Espectros",
                                 command=self.show_all_spectra)
        show_all_btn.pack(fill=tk.X, pady=2)

        fit_kittel_btn = ttk.Button(analysis_frame, text="Ajustar Kittel",
                                   command=self.fit_kittel_data)
        fit_kittel_btn.pack(fill=tk.X, pady=2)
        
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

        # Aba para resultados
        self.results_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.results_frame, text="Resultados")
        
        # Criar gráficos
        self.create_spectrum_plot()
        self.create_linewidth_plot()
        self.create_kittel_plot()

        # Criar tabela de resultados
        self.create_results_table()
        
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
        
    def create_results_table(self):
        """Cria a tabela de resultados"""
        # Frame para tabela
        table_frame = ttk.Frame(self.results_frame)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Treeview para tabela
        columns = ('Frequência (GHz)', 'Hr (Oe)', 'Hr (mT)', 'ΔH (Oe)', 'ΔH (mT)', 'Arquivo')
        self.results_tree = ttk.Treeview(table_frame, columns=columns, show='headings')
        
        # Configurar colunas
        for col in columns:
            self.results_tree.heading(col, text=col)
            self.results_tree.column(col, width=100)
        
        # Scrollbars
        v_scrollbar = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.results_tree.yview)
        h_scrollbar = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        # Pack tudo
        self.results_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
    
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

        # Atualizar status final
        if loaded_count > 0:
            self.status_var.set(f"✓ {len(self.loaded_files)} arquivo(s) carregado(s) ({loaded_count} novos)")
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
                self.ax.axvline(file_info['Hr'], color='red', linestyle='--', 
                               linewidth=2, label=hr_label)
        
        # Mostrar ajuste se solicitado
        if file_info['processed'] and self.show_fit_var.get():
            if 'fitted_curve' in file_info:
                fit_label = f'Ajuste Lorentziano'
                self.ax.plot(file_info['campo'], file_info['fitted_curve'], 'r-', 
                            linewidth=2, alpha=0.7, label=fit_label)
        
        self.ax.set_xlabel('Campo Magnético (Oe)')
        self.ax.set_ylabel('dχ"/dH (u.a.)')
        self.ax.set_title(f'Espectro FMR - {file_info["frequency"]} GHz')
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
    
    def process_all_files(self):
        """Processa todos os arquivos carregados"""
        if not self.loaded_files:
            messagebox.showwarning("Aviso", "Nenhum arquivo carregado")
            return
        
        def process_thread():
            self.status_var.set("Processando arquivos...")
            
            # Limpar tabela de resultados
            for item in self.results_tree.get_children():
                self.results_tree.delete(item)
            
            for i, file_info in enumerate(self.loaded_files):
                try:
                    self.status_var.set(f"Processando {i+1}/{len(self.loaded_files)}: {file_info['filename']}")
                    
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
                    
                    # Adicionar à tabela de resultados
                    freq = file_info['frequency']
                    Hr_mt = Hr * 0.1  # Oe -> mT
                    Delta_H_mt = Delta_H * 0.1  # Oe -> mT
                    
                    self.results_tree.insert('', 'end', values=(
                        f"{freq:.1f}",
                        f"{Hr:.1f}",
                        f"{Hr_mt:.1f}",
                        f"{Delta_H:.1f}",
                        f"{Delta_H_mt:.1f}",
                        file_info['filename']
                    ))
                    
                    # Armazenar no analyzer para compatibilidade
                    self.analyzer.spectra_data[freq] = {
                        'campo': file_info['campo'],
                        'sinal': file_info['sinal'],
                        'Hr': Hr,
                        'Delta_H': Delta_H
                    }
                    self.analyzer.resonance_fields[freq] = Hr
                    
                except Exception as e:
                    print(f"Erro ao processar {file_info['filename']}: {e}")
            
            self.status_var.set(f"Processamento concluído: {len(self.loaded_files)} arquivos")

            # Atualizar plots
            self.update_linewidth_plot()

        # Executar em thread separada para não travar a interface
        thread = threading.Thread(target=process_thread)
        thread.daemon = True
        thread.start()
    
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
            
            ax.plot(file_info['campo'], file_info['sinal'], 'b-', linewidth=1)
            
            if 'Hr' in file_info:
                ax.axvline(file_info['Hr'], color='red', linestyle='--', linewidth=1.5)
            
            title = f'{file_info["frequency"]} GHz'
            if 'Delta_H' in file_info:
                title += f'\nΔH = {file_info["Delta_H"]:.1f} Oe'
            
            ax.set_title(title, fontsize=10)
            ax.set_xlabel('Campo (Oe)', fontsize=8)
            ax.set_ylabel('dχ"/dH', fontsize=8)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=8)
        
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
                                  label=f'Ajuste linear (slope = {coeffs[0]:.3f} mT/GHz)', zorder=2)

            # Calcular R²
            linewidth_pred = np.polyval(coeffs, frequencies)
            ss_res = np.sum((linewidths - linewidth_pred) ** 2)
            ss_tot = np.sum((linewidths - np.mean(linewidths)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)

            # Adicionar texto com informações do ajuste
            textstr = f'Ajuste Linear:\ny = {coeffs[0]:.3f}x + {coeffs[1]:.3f}\nR² = {r_squared:.4f}'
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

        self.linewidth_fig.tight_layout()
        self.linewidth_canvas.draw()

    def fit_kittel_data(self):
        """Executa ajuste de Kittel e atualiza plot"""
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

            # Fazer ajuste no plano
            self.kittel_params, fitted_fields = self.fmr_fitting.fit_inplane_data(frequencies, fields_tesla)

            # Atualizar plot de Kittel
            self.update_kittel_plot(frequencies, fields_tesla, fitted_fields)

            # Mostrar parâmetros
            params_text = f"""Parâmetros Ajustados (Geometria No Plano):
Ms = {self.kittel_params['Ms']:.3f} ± {self.kittel_params['Ms_err']:.3f} T
Ha = {self.kittel_params['Ha']*1000:.1f} ± {self.kittel_params['Ha_err']*1000:.1f} mT
Hk = {self.kittel_params['Hk']*1000:.1f} ± {self.kittel_params['Hk_err']*1000:.1f} mT"""

            r2 = self.fmr_fitting.calculate_r_squared(fields_tesla, fitted_fields)
            params_text += f"\nR² = {r2:.6f}"

            messagebox.showinfo("Ajuste de Kittel", params_text)

        except Exception as e:
            messagebox.showerror("Erro", f"Erro no ajuste de Kittel:\n{str(e)}")

    def update_kittel_plot(self, frequencies, fields_experimental, fields_fitted):
        """Atualiza o plot do ajuste de Kittel"""
        self.kittel_ax.clear()

        # Converter unidades para plot
        freq_ghz = frequencies / 1e9
        field_mt_exp = fields_experimental * 1000
        field_mt_fit = fields_fitted * 1000

        # Plot dados experimentais e ajuste
        self.kittel_ax.scatter(freq_ghz, field_mt_exp, color='red', s=80, alpha=0.8,
                              label='Dados experimentais', zorder=3)
        self.kittel_ax.plot(freq_ghz, field_mt_fit, 'b-', linewidth=3,
                           label='Ajuste teórico (Kittel)', zorder=2)

        self.kittel_ax.set_xlabel('Frequência (GHz)', fontsize=12)
        self.kittel_ax.set_ylabel('Campo de Ressonância (mT)', fontsize=12)
        self.kittel_ax.set_title('Ajuste de Kittel - Geometria No Plano', fontsize=14)
        self.kittel_ax.legend()
        self.kittel_ax.grid(True, alpha=0.3)

        # Adicionar texto com parâmetros se disponível
        if self.kittel_params:
            r2 = self.fmr_fitting.calculate_r_squared(fields_experimental, fields_fitted)
            textstr = f"""Parâmetros ajustados:
Ms = {self.kittel_params['Ms']:.3f} ± {self.kittel_params['Ms_err']:.3f} T
Ha = {self.kittel_params['Ha']*1000:.1f} ± {self.kittel_params['Ha_err']*1000:.1f} mT
Hk = {self.kittel_params['Hk']*1000:.1f} ± {self.kittel_params['Hk_err']*1000:.1f} mT

R² = {r2:.6f}"""

            props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
            self.kittel_ax.text(0.02, 0.98, textstr, transform=self.kittel_ax.transAxes,
                               fontsize=10, verticalalignment='top', bbox=props)

        self.kittel_fig.tight_layout()
        self.kittel_canvas.draw()

def main():
    """Função principal para executar a GUI"""
    root = tk.Tk()
    FMRGuiApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()