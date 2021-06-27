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
    - https://google-earth-pro.gosur.com
    - https://www.windy.com/
    - http://pt.windfinder.com
    - http://www.climatempo.com.br
- **Latitude** e **Longitude** decimal da impactação (em °);
- **Azimute** estimado do disparo;
  - Geralmente próximo ao Azimute da impactação;

O programa utiliza as soluções encontradas para o caso sem arrasto como ponto de partida, pois podemos estimar a condição de ângulo inicial do disparo **θ0**.

Caso o projétil com arrasto chegue no ponto final com uma angulação menor que **ϕ**, ou seja, **ϕf < ϕ**, incrementamos o ângulo incial de disparo **θ0** com um fator e recalculamos toda a trajetória. Por outro lado, se ao fim **ϕf > ϕ**, reduzimos **θ0** com um fator e recalculamos toda a trajetória.

O critério de parada do cálculo é que ![image](https://user-images.githubusercontent.com/86118560/122674610-4e1c7500-d1ac-11eb-90d4-e1afafdb3a1f.png), sendo um **δϕ** um valor arbitrário aceitável. 

O mesmo cálculo é realizado para o Azimute.

## Exemplo:
Vamos supor que um projétil de arma de fogo atingiu o Edifício Maria Juliana, localizado na Av. Boa Viagem, 4398, a 89 m do solo, no último pavimento.

Para Localização e imagens vamos utilizar o https://earth.google.com/web/.

![01](https://user-images.githubusercontent.com/86118560/123481188-c1255180-d5d9-11eb-9a16-ed9b1921209c.jpg)

Ao fazer a linha de tiro da janela ao interior da residência, foram medidos os seguintes parâmetros:
| Parâmetros da Impactação        | Medidas      |
|---------------------------------|--------------|
| Altura da Impactação            | 89 m         |
| Ângulo ϕ com a Horizontal       | 4 °          |
| Trajetória Era Descendente?     | Não          |
| Ângulo γ com o Norte (azimute)  | 183 °        |
| Latitude Decimal da Impactação  | -8.127727 °  |
| Longitude Decimal da Impactação | -34.898383 ° |

Após recuperação do projétil, localizado no interior da residência, ficou constatado tratar-se de um projétil CBC, NTA .38 SPL EOPP 158gr.
|  Parametros do Projétil           | Medidas    |
|-----------------------------------|------------|
| Massa do projétil                 | 10.240 g   |
| Velocidade inicial do projétil    | 230 m/s    |
| Diâmetro do projétil já disparado | 8.82 mm    |
| Coeficiente de Arrasto Cd         | 0.235800   |
| Dextrogiro ou Levogiro            | Dextrogiro |

> O Coeficiente de Arrasto do projétil pode ser calculado utilizando o programa Forensic-GPM-Trajectories/Coeficiente de Arrasto/Cd_projeteis.c

> É possível rodar online, copiando o código para: https://www.onlinegdb.com/online_c_compiler

Após informações sobre a época do fato (o disparo hipotético ocorreu no dia 01/01/21 as 03h da madrugada), foram registradas as condições do vento no local. Nesse exemplo foi utilizado o https://earth.nullschool.net/pt/

![02](https://user-images.githubusercontent.com/86118560/123483006-7c4eea00-d5dc-11eb-8635-68d998fc2b81.jpg)

![03](https://user-images.githubusercontent.com/86118560/123483012-7eb14400-d5dc-11eb-8f80-4463bfac3b76.jpg)

![04](https://user-images.githubusercontent.com/86118560/123483016-807b0780-d5dc-11eb-863f-d4a5ec3bcf2f.jpg)

| Parametros do Vento | Medidas   |
|---------------------|-----------|
| Velocidade do Vento | 10 km/h   |
| Direção do Vento    | 100 °     |

O programa foi aberto, alimentado com os parâmetros fornecendo o seguinte resultado parcial:
>O ângulo θ do início do disparo com a horizontal considerado a partir do solo e em um sistema sem arrasto ou outras correcoes vale: 12.04°
>
>Considerando um sistema sem arrasto: x- = 559.208 m e x+ = 1646.637 m
>Existe alguma edificacao entre o impacto e a possível origem do disparo no solo de onde possa ter partido o tiro?
>|              | Latitude  | Longitude  |
>|--------------|-----------|------------|
>| Impacto      | -8.127727 | -34.898383 |
>| Origem (NMM) | -8.122955 | -34.898146 |

Ao examinar as coordenadas de possível origem observou-se que havia uma edificação próxima.

![05](https://user-images.githubusercontent.com/86118560/123485984-9939ec00-d5e1-11eb-9b85-d6e80a7a91b7.jpg)

![06](https://user-images.githubusercontent.com/86118560/123486256-0d748f80-d5e2-11eb-88e7-65ce055537c6.jpg)

Após responder Sim (1) ao programa, foi inserido as coordenadas e altura dessa edificação:
| Edificacao | Características |
|------------|-----------------|
| Latitude   | -8.123596 °     |
| Longitude  | -34.898176 °    |
| Altura     | 60 m            |

O programa conclui que, ao invés da origem estimada, partiu, de fato, da edificação indicada, com o segunite resultado:

| Resultado                     | Valores        |
|-------------------------------|----------------|
| Downrange Total               | 458.996 m      |
| Altura de impactação          | 88.974 m       |
| Desvios para direita          | 1.292 m        |
| Ângulo θ (inicial) do disparo | 11.71º         |
| Azimute inicial do disparo    | 182.75º        |
| Latitude do disparo           | -8.123609 N/S  |
| Longitude do disparo          | -34.898176 L/O |
| Altura do Disparo             | 20.87 m        |

Assim, procurou-se as coordenadas e altura fornecida pelo programa:

![07](https://user-images.githubusercontent.com/86118560/123486950-63960280-d5e3-11eb-8d8f-ac8c8d8aece2.jpg)

A região que ocorreu o disparo é dada na posição do mouse.

O programa fornece, por fim, a estimativa de duração do tempo de deslocamento do projétil, Velocidade Final e Energia Cinética para comparação com valores nominais de outros calibres.

| Resultado                     | Valores        |
|-------------------------------|----------------|
| Tempo de Trajetória           | 2.5 segundos   |
| Velocidade Final              | 151.45 m/s     |
| Energia Cinética              | 117.44 J       |

> É importante destacar que, caso haja ainda outras edificações ainda mais próximas da impactação é necessário rodar novamente o programa com essas novas coordenadas e altura até não sobrar dúvidas.

Vamos supor que havia uma suspeita que o projétil pudesse ser disparado a partir do Edf. D Pedro I, na Rua Dos Navegantes, 768, dado encontrar-se no caminho entre a origem e a impactação.

| Edificacao | Características |
|------------|-----------------|
| Latitude   | -8.125346 °     |
| Longitude  | -34.898260 °    |
| Altura     | 55 m            |

![01_1](https://user-images.githubusercontent.com/86118560/123550644-e5ea0800-d744-11eb-8602-18e65dd0fa21.jpg)

Com esses parâmetros o programa retorna o seguinte resultado:
> O projétil passou por cima da edificação fornecida.

Considerando que partiu, de fato, da primeira edificação fornecida, ao analisar o Output `data` do programa, próximo às coordenadas fornecidas para a segunda edificação, é possivel observar que:

| Tempo  | Latitude        | Longitude        | Altura      |
|--------|-----------------|------------------|-------------|
| …      | …               | …                | …           |
| 0.937 s| -8.125346086077 | -34.898260159499 | 57.005053 m |
| **0.937 s**| **-8.125346256627** | **-34.898260167993** | **57.008145 m** |
| 0.937 s| -8.125346427174 | -34.898260176488 | 57.011236 m |
| …      | …               | …                | …           |

Ou seja, na região de interesse, o projétil passou 2 metros acima da edificação.
