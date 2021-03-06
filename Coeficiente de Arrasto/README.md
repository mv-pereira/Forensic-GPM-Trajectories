# Estimativas de Coeficiente de Arrasto Cd baseado em perda energética.

O programa *Cd_projeteis.c* recebe dados fornecidos e retorna o Coeficiente de Arrasto adimensional Cd para cálculo de trajetórias GPMT.

O programa requisita as seguintes informações para funcionamento:
- **Velocidade inicial** (em m/s) do projétil;
  - Estimada com base nas características do projétil ou outras considerações;
- **Distancia** (em m) o Cd será calculado;
  - Valores fornecidos pelo fabricante. Pela CBC, temos a tabela: *Projeteis_CBC.md*;
- **Velocidade final** (em m/s) do projétil;
- **Massa** do projétil (em g);
  - Se o projétil perdeu massa na impactação é necessário conhecer a massa original;
- **Diâmetro** do projétil (em mm);
  - Este diâmetro se refere ao diâmetro do projétil após passagem pelo cano da arma, pode ser medido diretamente do projétil ou de um outro de referência;

O programa inicia o cálculo com um Cd estimado de **0.2** e velocidade inicial fornecida.

Após **x** m, se o projétil terminar com mais velocidade que o fornecido pelo fabricante, o programa incrementa Cd (aumenta o arrasto) e recalcula.

Se após **x** m o projétil o terminar com menor velocidade, o programa diminui Cd (diminui o arrasto) e recalcula.
O cálculo é repetido até que o **Δv < ε** (velocidade final estimada pelo programa e a velocidade final fornecida pelo fabricante).

Abaixo, um fluxograma (https://www.zenflowchart.com) do procedimento adotado pelo programa.
<img src="https://github.com/mv-pereira/Forensic-GPM-Trajectories/blob/main/Coeficiente%20de%20Arrasto/flowchart_Cd.png" alt="1" width="822" height="439">


## Exemplo

Para um projétil encamisado .38 foi calculado:

[mario@g7 2.1]$ ./Cd

>Programa para cálculo do Coeficiente de Arrasto Cd (adimensional) para projéteis de dimensões conhecidas utilizando parâmetros fornecidos pelos fabricantes de projétis.
>Digite a velocidade inicial (em m/s) do projétil:
>230
>
>Digite para qual distancia (em m) o Cd será calculado:
>100
>
>Digite a velocidade final (em m/s) do projétil:
>211
>
>Digite a massa (em g) do projétil:
>10.24
>
>Digite o diâmetro (em mm) do projétil:
>8.82
>
>Considerando uma perda de velocidade de 19.00 m/s em 100.00 m,
>O coeficiente de Arrasto vale: **Cd=0.235800**.
