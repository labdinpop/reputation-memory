# S3 — comparação isolada da seleção

Autorização explícita dos autores: mensagem de 5 de outubro de 2026, item 1.b.

Seleção: fator r_j / soma r_k da Eq. (5) de Li et al. (2025), PDF p. 3. Na soma nula, seleção uniforme por convenção explicitamente registrada. Foram mantidos K=0.1, payoff normalizado, inicialização verdadeira, observação antes da imitação, sincronismo e ausência de erros. O modelo-base foi preservado intacto; a cópia isolada está em `src/model_s3.py`.

Grade completa: 7 valores de M × 30 realizações = 210; L=100, b_A=1.15, epsilon=0; T0=5000, T1=2000, amostragem a cada 10 rodadas. Nenhum ponto exigiu repetição pelo critério de deriva.

| M | Média f_C | Erro-padrão | Todos C | Todos D |
|---|---|---|---|---|
| 1 | 1.000000 | 0.000000 | 30/30 | 0/30 |
| 2 | 1.000000 | 0.000000 | 30/30 | 0/30 |
| 3 | 1.000000 | 0.000000 | 30/30 | 0/30 |
| 5 | 1.000000 | 0.000000 | 30/30 | 0/30 |
| 9 | 1.000000 | 0.000000 | 30/30 | 0/30 |
| 17 | 0.619040 | 0.000263 | 0/30 | 0/30 |
| 33 | 0.441221 | 0.000306 | 0/30 | 0/30 |

M=1,2,3,5,9 empatam em cooperação completa. M=17 e M=33 apresentam cooperação menor. A seleção proporcional aumenta a cooperação relativa à mistura original, mas não produz um máximo intermediário superior ao endpoint M=1 nesta grade. As diferenças de sensibilidade de aprendizagem, payoff e demais convenções dos artigos de Li continuam fora desta intervenção isolada.

A rotina reutiliza a saída antecipada do modelo-base para estados absorventes sem erros antes da janela de medida; f_C permanece exatamente constante nesses casos. Não são salvas trajetórias completas. Os arquivos por realização contêm sementes, parâmetros e medidas do protocolo. A verificação de sementes e identidade está em `qa/S3/verification.json`; custo medido em `qa/S3/runtime.json`.
