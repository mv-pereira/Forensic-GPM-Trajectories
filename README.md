# Forensic-GPM-Trajectories
Estimativa de origem de disparos subsônicos baseado em Trajetórias de Massas Pontuais Generalizadas a partir de Projéteis Dinamicamente estáveis em toda a trajetória utilizando a linguagem de programação C.

Para uma revisão breve do conteúdo proposto, ver: **GPMT.md**

O programa requisita as seguintes informações para funcionamento:
- Ângulo **ϕ** (em °) com a horizontal no que o projétil impactou;
- Ângulo **γ** (em °) com o Norte (azimute) no que o projétil impactou;
- **Velocidade inicial** (em m/s) do projétil;
  - Estimada com base nas características do projétil ou outras considerações;
- **Massa** do projétil (em g);
  - Se o projétil perdeu massa na impactação é necessário conhecer a massa original;
- **Diâmetro** do projétil (em mm);
  - Este diâmetro se refere ao diâmetro do projétil após passagem pelo cano da arma, pode ser medido diretamente do projétil ou de um outro de referência;
- Coeficiente de Arrasto (**Cd**);
  - Aproximadamente 0.2 para projéteis subsônicos ou pode-se utilizar o programa **Forensic-GPM-Trajectories/Coeficiente de Arrasto/Cd_projeteis.c**;
- Se o projetil é **dextrogiro** ou **levogiro**;
  - Essa informação permitirá estimar o *wind drift* do projétil com base na fórmula de Bryan Litz e The Overwatch (https://theoverwatch.wixsite.com/theoverwatch/post/spin-drift);
- Módulo da **velocidade do vento** (em km/h) e sua **direção** com relação ao Norte (em °);
  - Essas informações podem ser obtidas de diversas fontes, com atenção para estimativa da época do disparo;
    - https://earth.nullschool.net/
    - https://www.windy.com/
    - http://pt.windfinder.com
    - http://www.climatempo.com.br
- **Latitude** e **Longitude** decimal do disparo (em °);
- **Azimute** estimado do disparo;
  - Geralmente próximo ao Azimute da impactação;

O programa utilizarmos as soluções encontradas para o caso sem arrasto como ponto de partida, pois podemos estimar a condição de ângulo inicial do disparo θ0.

Caso o projétil com arrasto chegue no ponto final com uma angulação menor que ϕ, ou seja, ϕf < ϕ, incrementamos o ângulo incial de disparo θ0 com um fator e recalculamos toda a trajetória. Por outro lado, se ao fim ϕf > ϕ, reduzimos θ0 com um fato e recalculamos toda a trajetória.

O critério de parada do cálculo é que ![image](https://user-images.githubusercontent.com/86118560/122674610-4e1c7500-d1ac-11eb-90d4-e1afafdb3a1f.png), sendo um δϕ um valor arbitrário aceitável. 

O mesmo cálculo é realizado para o Azimute.

## Exemplo:
Vamos supor que um disparo atingiu o 11° andar de um edifício no final da Avenida Caxangá, bairro da Madalena, Recife/PE.

Ao medir as informações, a equipe obteve os seguintes valores:

y0 = 36 m; ϕ = 8°; γ = 20°;

Velocidade do Projétil .38 = 305 m/s;

Cd (.38) = 0.2;

Massa (.38) = 8 g;

Diâmetro (.38) = 8.82 mm;

Velocidade do vento, direção e Latitude são dadas pela carta dos ventos.


<img src="https://user-images.githubusercontent.com/86118560/122675558-5676af00-d1b0-11eb-9eeb-a7482cd5cb26.png" alt="1" width="711" height="400">

Valores da velocidade do vento em Recife em 31/03/2020 as 15:00.
Fonte: https://pt.windfinder.com/#16/-8.0570/-34.9108/report


[mario@g7 2.2]$ gcc Trajetoria_Balistica.c -o teste -lm

[mario@g7 2.2]$ ./teste  

Digite a altura (em metros) da impactação em relação ao disparo: 36 
A impactação ocorreu em trajetória descendente? 
1 - SIM 2 - NÃO. 
1 

Digite o ângulo ϕ com a horizontal (em °): 8 

Digite o ângulo γ com o Norte (em °) (azimute): 20 

Digite a velocidade inicial (em m/s) do projétil: 305 

Digite a massa do projétil (em g): 8 

Digite o diâmetro do projétil já disparado (em mm): 8.82 

Digite o Coeficiente de Arrasto Cd (aproximadamente 0.2 em casos subsônicos) do projétil: 0.2 

O projetil é dextrogiro ou levogiro? 1 - Dextrogiro.    2 - Levogiro. 
1 

Digite o módulo da velocidade do vento (em km/h): 17 

Digite a direção do vento em relação ao Norte (em °): 131 

Digite a latitude decimal do disparo (em °), Valor precisa ser negativo, se o disparo ocorrer He
misfério Sul: -8.0570 

Digite a longitude decimal do disparo (em °): -34.9106 

Digite o azimute estimado do disparo (em °), medido em sentido horário em relação ao Norte: 20 

O ângulo θ do início do disparo com a horizontal considerado a partir do solo e em um sistema se
m arrasto ou outras correcoes vale: 9.42° 

Considerando um sistema sem arrasto: x- = 234.937 m e x+ = 2836.528 m 

Considerando um sistema com arrasto, os cálculos terminaram com os seguintes valores: 
Foram efetuados 261 cálculos de trajetória. 
Downrange Total = 957.355 m. 
Altura de impactação = 36.000 m. 
Desvios para esquerda = 8.644 m. 
Ângulo θ (inicial) do disparo = 7.93º. 
Azimute inicial do disparo = 21.11º. 

A trajetória teve outro desvio devido ao spindrift de, aproximadamente, 226 cm para direta, não 
incluidos nos cálculos. 

Considerando um disparo a partir do solo, o projétil partiu, aproximadamente, das coordenadas: -
8.065051 N/S, -34.913595 L/O. 

O projétil levou cerca de 5.1 segundos do momento do disparo à impactação. 
Possuía velocidade final de 125.22 m/s e energia cinética de 62.72 J.

Energia cinética de alguns projéteis para comparação: 
![image](https://user-images.githubusercontent.com/86118560/122675812-7f4b7400-d1b1-11eb-9e46-b584f982daf7.png)

Argumentos no **Gnuplot** para visualização da trajetória:

set xrange [0:900]  

set yrange [900:0]

set zrange [0:100] 

set parametric 

set urange [0:1000]

set vrange [0:50]

set style line 1 lw 3

splot "data" u 2:4:3 with lines ls 1, u,0,v with lines

<img src="https://user-images.githubusercontent.com/86118560/122675566-6098ad80-d1b0-11eb-9cd1-c935540e481e.png" alt="1" width="711" height="400">


<img src="https://user-images.githubusercontent.com/86118560/122675568-62fb0780-d1b0-11eb-942b-951683c8bb39.png" alt="1" width="711" height="400">
