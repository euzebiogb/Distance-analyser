# -*- coding: utf-8 -*-
"""
Este script gera um vídeo de heatmap temporal a partir de um arquivo PDB de dinâmica molecular.

O heatmap mostra a distância entre cada átomo de um ligante especificado e o 
ponto médio da cadeia lateral de cada resíduo da proteína, frame a frame.

O mapa de cores padrão é azul (perto) para vermelho (longe).

Pré-requisitos:
- Python 3
- Bibliotecas: MDAnalysis, numpy, matplotlib, seaborn, Pillow
- FFmpeg: necessário para salvar o vídeo em formato .mp4.

Como executar no terminal:
# Análise padrão
python temporal_heatmap.py -i seu_arquivo.pdb -l LIG -o video_normal.mp4

# Análise de DIFERENÇA (mostra a mudança em relação ao primeiro frame)
python temporal_heatmap.py -i seu_arquivo.pdb -l LIG -o video_diff.mp4 --diff

# Análise de ALTA VARIABILIDADE (destaca interações e rótulos mais dinâmicos)
python temporal_heatmap.py -i seu_arquivo.pdb -l LIG -o video_highvar.mp4 --highvar

# Análise com EXAGERO VISUAL (destaca mais as mudanças de distância)
python temporal_heatmap.py -i seu_arquivo.pdb -l LIG -o video_exagerado.mp4 --exaggerate

# Combinando DIFERENÇA com ALTA VARIABILIDADE
python temporal_heatmap.py -i arq.pdb -l LIG -o video_combo.mp4 --diff --highvar
"""

import argparse
import warnings
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns

# Suprime avisos comuns do MDAnalysis para uma saída mais limpa
warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')

def parse_arguments():
    """Analisa os argumentos da linha de comando."""
    parser = argparse.ArgumentParser(
        description="Cria um vídeo de heatmap de distância entre um ligante e resíduos de proteína a partir de um arquivo PDB.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        dest='pdb_file',
        help="Caminho para o arquivo PDB de entrada contendo a trajetória."
    )
    parser.add_argument(
        '-l', '--ligand',
        required=True,
        dest='ligand_name',
        help="Nome do resíduo do ligante (ex: LIG)."
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        dest='output_file',
        help="Nome do arquivo de saída (ex: heatmap_movie.mp4 ou animacao.gif)."
    )
    parser.add_argument(
        '--format',
        default='mp4',
        choices=['mp4', 'gif'],
        dest='output_format',
        help="Formato do arquivo de saída. 'mp4' requer ffmpeg. 'gif' requer a biblioteca Pillow."
    )
    parser.add_argument(
        '-s', '--sample',
        type=int,
        default=100,
        dest='sample_percentage',
        help="Porcentagem de frames para amostrar (1-100). O último frame é sempre incluído. (padrão: 100)"
    )
    parser.add_argument(
        '--highvar',
        action='store_true',
        help="Ativa o modo de alta variabilidade, usando um método estatístico (IQR) para mostrar apenas as interações que mais variam com o tempo."
    )
    parser.add_argument(
        '--exaggerate',
        action='store_true',
        help="Ativa o modo de exagero visual. Usa uma faixa de cores mais estreita (percentil 5-95)\npara destacar mais as variações de distância."
    )
    parser.add_argument(
        '--diff',
        action='store_true',
        help="Ativa o modo de diferença. O heatmap mostra a mudança na distância em relação ao primeiro frame."
    )
    return parser.parse_args()

def create_heatmap_video(pdb_file, ligand_name, output_file, output_format, sample_percentage, high_variability_mode, exaggerate_mode, diff_mode):
    """
    Função principal para carregar dados, calcular distâncias e gerar o vídeo.
    """
    try:
        print(f"Carregando o universo do arquivo: {pdb_file}...")
        universe = mda.Universe(pdb_file)
    except Exception as e:
        print(f"Erro ao carregar o arquivo PDB: {e}")
        return

    # --- 1. Seleção dos Átomos e Resíduos ---
    print("Realizando seleções de átomos e resíduos...")
    try:
        ligand_atoms = universe.select_atoms(f"resname {ligand_name}")
        if not ligand_atoms.n_atoms:
            print(f"Erro: Nenhum átomo encontrado para o nome do ligante '{ligand_name}'. Verifique o nome do resíduo.")
            return
            
        protein_residues = universe.select_atoms("protein").residues
        protein_residues = [res for res in protein_residues if res.resname != ligand_name]
        
        if not protein_residues:
            print("Erro: Nenhum resíduo de proteína encontrado no arquivo.")
            return
    except Exception as e:
        print(f"Erro durante a seleção de átomos: {e}")
        return

    # --- 2. Amostragem de Frames ---
    n_frames = len(universe.trajectory)
    if not (1 <= sample_percentage <= 100):
        print("Aviso: A porcentagem de amostragem deve estar entre 1 e 100. Usando 100%.")
        sample_percentage = 100
    
    num_samples = max(1, int(n_frames * sample_percentage / 100.0))
    frame_indices = np.unique(np.linspace(0, n_frames - 1, num=num_samples, dtype=int))

    # --- 3. Pré-cálculo de todas as distâncias ---
    print(f"Pré-calculando distâncias para {len(frame_indices)} de {n_frames} frames ({sample_percentage}%)... Isso pode levar um tempo.")
    
    all_distances_list = []
    for i, frame_idx in enumerate(frame_indices):
        print(f"  - Processando frame original {frame_idx + 1}/{n_frames}", end='\r')
        ts = universe.trajectory[frame_idx]
        distance_matrix = np.zeros((len(ligand_atoms), len(protein_residues)))
        lig_coords = ligand_atoms.positions

        for j, residue in enumerate(protein_residues):
            sidechain_atoms = residue.atoms.select_atoms("not backbone and not name H*")
            if len(sidechain_atoms) == 0:
                ref_point = residue.atoms.select_atoms("name CA").center_of_geometry()
            else:
                ref_point = sidechain_atoms.center_of_geometry()
            
            distances = np.linalg.norm(lig_coords - ref_point, axis=1)
            distance_matrix[:, j] = distances
        
        all_distances_list.append(distance_matrix)
    
    print("\nCálculo de distâncias concluído.")

    all_distances = np.array(all_distances_list)
    
    # --- 4. Processamento dos dados e Escala de Cores ---
    processed_data = all_distances
    cbar_label = 'Distância (Å)'

    if diff_mode:
        print("\n--- MODO DE DIFERENÇA ATIVO ---")
        if exaggerate_mode:
            print("Aviso: O modo --exaggerate é ignorado quando --diff está ativo.")
        
        reference_distances = all_distances[0]
        processed_data = all_distances - reference_distances
        cbar_label = 'Diferença de Distância (Å)'
        
        # Escala de cores simétrica em torno de zero
        limit = np.max(np.abs(processed_data))
        vmin, vmax = -limit, limit
        print(f"Intervalo da diferença (min-max): {-limit:.2f} Å - {limit:.2f} Å")
        print("--------------------------------")

    else: # Modo normal (distância absoluta)
        vmin, vmax = np.min(all_distances), np.max(all_distances)
        print(f"Intervalo de distância global (min-max): {vmin:.2f} Å - {vmax:.2f} Å")
        
        if exaggerate_mode:
            print("\n--- MODO DE EXAGERO VISUAL ATIVO ---")
            vmin = np.percentile(all_distances, 5)
            vmax = np.percentile(all_distances, 95)
            print(f"A escala de cores foi ajustada para o intervalo de percentil 5-95: {vmin:.2f} Å - {vmax:.2f} Å.")
            print("--------------------------------------")

    # --- 5. Geração da Máscara de Variabilidade (se ativado) ---
    variability_mask = None
    if high_variability_mode:
        if len(all_distances) < 2:
            print("Aviso: O modo de alta variabilidade requer pelo menos 2 frames. Desativando o modo.")
        else:
            print("\n--- MODO DE ALTA VARIABILIDADE ATIVO ---")
            # A variabilidade é sempre calculada com base nas distâncias originais
            variability = np.std(all_distances, axis=0)
            q1 = np.percentile(variability, 25)
            q3 = np.percentile(variability, 75)
            iqr = q3 - q1
            cutoff_value = q3 + 1.5 * iqr
            variability_mask = variability < cutoff_value
            num_shown = np.sum(variability >= cutoff_value)
            percent_shown = (num_shown / variability.size) * 100
            print(f"Análise de variabilidade (Desvio Padrão) via IQR:")
            print(f"  - Limiar de alta variabilidade (Q3 + 1.5*IQR): {cutoff_value:.2f} Å")
            print(f"Destacando {num_shown} de {variability.size} interações ({percent_shown:.1f}% do total).")
            print("-----------------------------------------")

    # --- 6. Geração da Animação ---
    ligand_atom_names = [atom.name for atom in ligand_atoms]
    protein_residue_labels = [f"{res.resname}{res.resid}" for res in protein_residues]

    if high_variability_mode and variability_mask is not None:
        print("Limpando os rótulos do eixo X para mostrar apenas os resíduos de alta variabilidade...")
        is_low_variability_residue = np.all(variability_mask, axis=0)
        protein_residue_labels = [
            label if not is_low_var else "" 
            for label, is_low_var in zip(protein_residue_labels, is_low_variability_residue)
        ]

    fig = plt.figure(figsize=(18, 10))
    ax = fig.add_axes([0.1, 0.25, 0.8, 0.65])
    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.65])

    def update(frame):
        ax.clear()
        cbar_ax.clear()
        
        original_frame_number = frame_indices[frame]
        display_matrix = processed_data[frame]
        
        sns.heatmap(display_matrix, ax=ax, cbar_ax=cbar_ax, cmap="coolwarm", 
                    vmin=vmin, vmax=vmax, mask=variability_mask,
                    cbar_kws={'label': cbar_label})
        
        title = f"Heatmap de Distância Temporal | Ligante ({ligand_name}) vs Proteína | Frame da Trajetória: {original_frame_number + 1}"
        if diff_mode: title += " (Modo de Diferença)"
        if high_variability_mode: title += " (Modo de Alta Variabilidade)"
        if exaggerate_mode and not diff_mode: title += " (Modo de Exagero Visual)"
        ax.set_title(title)
        
        ax.set_xlabel("Resíduo da Proteína (Ponto Médio da Cadeia Lateral)")
        ax.set_ylabel(f"Átomos do Ligante ({ligand_name})")
        
        ax.set_xticks(np.arange(len(protein_residue_labels)) + 0.5)
        ax.set_xticklabels(protein_residue_labels, rotation=90, fontsize=8)
        ax.set_yticks(np.arange(len(ligand_atom_names)) + 0.5)
        ax.set_yticklabels(ligand_atom_names, rotation=0, fontsize=8)
        
        print(f"  - Renderizando frame de saída {frame + 1}/{len(all_distances)}", end='\r')

    print(f"\nGerando animação e salvando em '{output_file}'...")
    ani = animation.FuncAnimation(fig, update, frames=len(all_distances), interval=150)
    
    try:
        if output_format == 'mp4':
            writer = animation.FFMpegWriter(fps=10, metadata=dict(artist='Python Script'), bitrate=1800)
            ani.save(output_file, writer=writer)
        elif output_format == 'gif':
            writer = animation.PillowWriter(fps=10)
            ani.save(output_file, writer=writer)
        print(f"\n\nArquivo salvo com sucesso em: {output_file}")
    except FileNotFoundError:
        if output_format == 'mp4':
            print("\n\nErro: 'ffmpeg' não encontrado. Por favor, instale o FFmpeg e certifique-se de que ele está no PATH do seu sistema para salvar o vídeo.")
        else:
            print("\n\nErro: Ocorreu um erro ao salvar o arquivo. Verifique se as dependências estão instaladas.")
    except Exception as e:
        print(f"\n\nOcorreu um erro ao salvar o arquivo: {e}")

if __name__ == "__main__":
    args = parse_arguments()
    create_heatmap_video(
        args.pdb_file, 
        args.ligand_name, 
        args.output_file, 
        args.output_format, 
        args.sample_percentage,
        args.highvar,
        args.exaggerate,
        args.diff
    )

