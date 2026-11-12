"""Report the authorized longer-transient repetitions, preserving previous results."""
from pathlib import Path
import csv,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def stats(path):
    with (ROOT/path).open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
    x=np.array([float(r['f_C']) for r in rows])
    return x.mean(),x.std(ddof=1)/np.sqrt(len(x)),np.mean([float(r['all_D']) for r in rows])
def main():
    data=json.loads((ROOT/'results/transient_refinement/completed.json').read_text())
    lines=['# Repetições direcionadas do transiente', '',
        'Foram selecionados somente os pontos ainda sinalizados que participavam dos cortes em b_A ou b_B, ou das grades E3–E6. Pontos reutilizados foram contados uma vez.', '',
        'Os seis pontos distintos foram repetidos com T0 = 20.000, mantendo os demais parâmetros, 30 realizações e a janela de medida de 2.000 rodadas. As saídas T0 = 5.000 e 10.000 permanecem intactas. Nenhum parâmetro novo ou extensão de ruído foi simulado.', '',
        '| Uso | Núcleo | b | ε | f_C em T0=10.000 | f_C em T0=20.000 | Todos-D final | Fração com diferença >0,02 |',
        '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in data:
        task=r['task'];p=task['parameters'];old=stats(task['old_source']);new=stats(r['new_source'])
        k=p['kernel']+'('+ (f"{p['a']:g}; {p['m']}" if p['kernel']=='PP' else str(p['M']))+')'
        lines.append(f"| {', '.join(task['used_by'])} | {k} | {p['b']:g} | {p['epsilon']:g} | {old[0]:.5f} ± {old[1]:.5f} | {new[0]:.5f} ± {new[1]:.5f} | {100*new[2]:.1f}% | {100*r['meta']['nonstationary_fraction']:.1f}% |")
    lines+=['','As incertezas são erros-padrão entre realizações. Cinco pontos permanecem sinalizados; E(17), b=1,30, ε=0 passou na checagem. Os novos arquivos foram selecionados nas tabelas, figuras e números do manuscrito, inclusive no reaproveitamento entre E5 e E6.', '',
        'As médias continuam a descrever o horizonte simulado. Não se considera comprovada a estacionariedade dos pontos ainda sinalizados, nem se interpreta ausência de absorção em uma amostra como ausência de risco.', '',
        'O enquadramento mantém os efeitos separados de a e m, a dispersão entre núcleos com G próximo e a ausência de ótimo intermediário apenas para ε ≤ 0,1. Nenhuma nova rodada foi lançada após estas repetições.', '',
        'Proveniência detalhada: `design.json`, `completed.json` e os arquivos `T0_20000.csv` indicados neles.']
    (ROOT/'results/transient_refinement/report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
if __name__=='__main__':main()
