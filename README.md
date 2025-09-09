# Heatmap de Distância Temporal para Dinâmica Molecular
Este repositório contém um script Python (temporal_heatmap.py) para gerar vídeos de heatmap a partir de trajetórias de dinâmica molecular. A ferramenta visualiza a evolução da distância entre os átomos de um ligante e os resíduos de uma proteína ao longo do tempo, oferecendo múltiplos modos de análise para destacar interações dinâmicas.

## Funcionalidades
Heatmap de Distância: Gera um vídeo mostrando as distâncias absolutas (Å) entre o ligante e a proteína.

Modo de Diferença (--diff): Mostra a variação da distância em relação ao primeiro frame, destacando aproximações (azul) e afastamentos (vermelho).

Modo de Alta Variabilidade (--highvar): Utiliza um filtro estatístico (IQR) para exibir apenas as interações que mais variam ao longo do tempo, limpando o gráfico e as legendas para focar nos eventos dinâmicos mais relevantes.

Modo de Exagero Visual (--exaggerate): Ajusta a escala de cores para "exagerar" visualmente as pequenas variações de distância, tornando-as mais perceptíveis.

Amostragem de Frames (-s): Permite analisar uma porcentagem da trajetória para gerar resultados mais rapidamente.

Múltiplos Formatos (--format): Suporte para saída em .mp4 (requer FFmpeg) ou .gif.

## 1. Pré-requisitos
Software
Python 3.6+

FFmpeg (necessário para salvar vídeos em formato .mp4). Para instalar em sistemas baseados em Debian/Ubuntu:

`sudo apt-get update && sudo apt-get install ffmpeg`

Bibliotecas Python
Você pode instalar todas as bibliotecas necessárias com o pip:

pip install MDAnalysis numpy matplotlib seaborn Pillow

## 2. Preparação do Arquivo PDB
Antes de executar o script, a trajetória (formato .xtc) precisa ser convertida para um arquivo PDB multi-frame. É crucial corrigir problemas de "salto" de caixa de simulação (periodic boundary conditions) para que as distâncias sejam calculadas corretamente.

Utilize o gmx trjconv do GROMACS com a opção -pbc nojump.

Comando:

`gmx trjconv -s PROD.tpr -f output_cpu.xtc -pbc nojump -o fitted.pdb -n index.ndx`

-s PROD.tpr: Seu arquivo de run input.

-f output_cpu.xtc: Sua trajetória.

-o fitted.pdb: O arquivo PDB de saída que será usado pelo script.

-n index.ndx: Seu arquivo de índice. Ao ser solicitado, selecione um grupo que contenha o ligante e a proteína para que ambos sejam incluídos no arquivo de saída.

## 3. Como Usar o Script
Execute o script temporal_heatmap.py a partir do seu terminal.

## Exemplos de Uso
Análise Padrão:
Gera um heatmap de distância com todos os frames.

`python temporal_heatmap.py -i fitted.pdb -l LIG -o movie_normal.mp4`

Análise de Diferença:
Mostra como as distâncias mudam em relação ao primeiro frame.

python temporal_heatmap.py -i fitted.pdb -l LIG -o movie_diff.mp4 --diff

Análise de Alta Variabilidade:
Foca apenas nas interações mais dinâmicas, limpando o gráfico.

`python temporal_heatmap.py -i fitted.pdb -l LIG -o movie_highvar.mp4 --highvar`

Combinando Modos:
Você pode combinar as opções para análises mais complexas. Por exemplo, analisar a diferença apenas das interações de alta variabilidade:

`python temporal_heatmap.py -i fitted.pdb -l LIG -o movie_diff_highvar.mp4 --diff --highvar`

Amostragem + Exagero Visual:
Analisa 20% dos frames com ênfase visual nas variações e salva como GIF:

`python temporal_heatmap.py -i fitted.pdb -l LIG -o movie_exaggerate.gif -s 20 --exaggerate --format gif`
