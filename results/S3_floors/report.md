# S3: piso na seleção proporcional

420 novas realizações concluídas; 210 realizações anteriores de piso zero preservadas. Todas usam L=100, b_A=1.15, epsilon=0, R=30, T0=5000, T1=2000, amostragem a cada 10 rodadas. Nenhum ponto exigiu extensão de transiente.

A seleção é (r_j+eta)/sum(r_k+eta). Para eta=0 e soma nula, a escolha é uniforme. Só a função de seleção difere da cópia anterior; o modelo principal permanece intacto.

Li (2026), PDF p. 4 Eq. (6), p. 7 seção 5.2 item 3: piso 10^-6. Portanto, 1e-6 e o piso de [17] são um único tratamento. Li (2025), p. 3 Eq. (5), não especifica piso nem convenção para divisão por zero.

| eta | M | f_C | Erro-padrão | Fração todos-C |
|---|---|---|---|---|
| 0 | 1 | 1.000 | 0.000000 | 1 |
| 0 | 2 | 1.000 | 0.000000 | 1 |
| 0 | 3 | 1.000 | 0.000000 | 1 |
| 0 | 5 | 1.000 | 0.000000 | 1 |
| 0 | 9 | 1.000 | 0.000000 | 1 |
| 0 | 17 | 0.619 | 0.000263 | 0 |
| 0 | 33 | 0.441 | 0.000306 | 0 |
| 1e-06 | 1 | 1.000 | 0.000000 | 1 |
| 1e-06 | 2 | 1.000 | 0.000000 | 1 |
| 1e-06 | 3 | 1.000 | 0.000000 | 1 |
| 1e-06 | 5 | 1.000 | 0.000000 | 1 |
| 1e-06 | 9 | 1.000 | 0.000000 | 1 |
| 1e-06 | 17 | 0.619 | 0.000288 | 0 |
| 1e-06 | 33 | 0.441 | 0.000293 | 0 |
| 0.01 | 1 | 1.000 | 0.000000 | 1 |
| 0.01 | 2 | 1.000 | 0.000000 | 1 |
| 0.01 | 3 | 1.000 | 0.000000 | 1 |
| 0.01 | 5 | 1.000 | 0.000000 | 1 |
| 0.01 | 9 | 0.861 | 0.000615 | 0 |
| 0.01 | 17 | 0.597 | 0.000259 | 0 |
| 0.01 | 33 | 0.418 | 0.000304 | 0 |

O platô de cooperação total vai até M=9 para eta=0 e 10^-6, e até M=5 para eta=0.01. Em M=9, eta=0.01 produz f_C=0.861. Logo, a extensão do platô depende do piso. Nenhum caso produz um máximo interior superior ao valor em M=1.

Quando r_j=0 e a soma das reputações vizinhas é positiva, a probabilidade de selecionar j é eta/(sum r_k+4 eta). Com piso desprezível ela é quase nula; se todas as reputações são zero, a escolha é uniforme. A interpretação do platô no artigo explicita essa condição.

Scripts: src/run_s3_floors.py (profile/run/status), src/analyze_s3_floors.py; configuração config/S3_floors.json; sementes e medidas por realização em results/S3_floors/points/. O custo medido está em qa/S3_floors/runtime.json.
