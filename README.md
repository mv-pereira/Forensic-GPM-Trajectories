# Forensic-GPM-Trajectories
Estimativa de origem de disparos subsônicos baseado em Trajetórias de Massas Pontuais Generalizadas a partir de Projéteis Dinamicamente estáveis em toda a trajetória utilizando a linguagem de programação C (Compilado Ubuntu 22.04.3 LTS, Kernel: 5.15.0-35-generic.).

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
- **Twist rate** do Projétil;
  - Baseado no cano do armamento;
  - http://www.taurususa.com
- **Comprimento** do projétil;
- Opção de exibir a contribuição do efeito Coriolis em separado;

O programa utiliza as soluções encontradas para o caso sem arrasto como ponto de partida, pois podemos estimar a condição de ângulo inicial do disparo **θ0**.

Caso o projétil com arrasto chegue no ponto final com uma angulação menor que **ϕ**, ou seja, **ϕf < ϕ**, incrementamos o ângulo incial de disparo **θ0** com um fator e recalculamos toda a trajetória. Por outro lado, se ao fim **ϕf > ϕ**, reduzimos **θ0** com um fator e recalculamos toda a trajetória.

O critério de parada do cálculo é que ![image](https://user-images.githubusercontent.com/86118560/122674610-4e1c7500-d1ac-11eb-90d4-e1afafdb3a1f.png), sendo um **δϕ** um valor arbitrário aceitável. 

O mesmo cálculo é realizado para o Azimute.

## Exemplo:
Vamos supor que um projétil de arma de fogo atingiu o Edifício Maria Juliana, localizado na Av. Boa Viagem, 4398, a 89 m do solo, no último pavimento.

Para Localização e imagens vamos utilizar o https://earth.google.com/web/.

![01](https://user-images.githubusercontent.com/86118560/123481188-c1255180-d5d9-11eb-9a16-ed9b1921209c.jpg)

![00](https://user-images.githubusercontent.com/86118560/123556326-f7d9a400-d760-11eb-84ce-71f57044ad95.jpg)

Após recuperação do projétil, localizado no interior da residência, ficou constatado tratar-se de um projétil CBC, NTA .38 SPL EOPP 158gr.
Ao fazer a linha de tiro da janela ao interior da residência, foram medidos os seguintes parâmetros:
Arquivo **input**:
```
# Este arquivo é usado para inserir dados para o cálculo de trajetórias de projéteis.
# Substitua os valores a seguir pelos de interesse para a sua análise.
# Cada linha contém um rótulo descritivo e um valor. Altere apenas o valor.
# Por exemplo, se a 'Altura de Impacto' é 100 metros, a linha deve ser 'Altura de Impacto (m): 100'.
# O programa calcula a trajetória de projéteis levando em conta fatores como gravidade, vento, etc.

Altura de Impacto (m): 89
Ângulo Phi (graus): 4
Ângulo Azimute (graus): 183
Latitude (graus): -8.127727
Longitude (graus): -34.898383
Velocidade Inicial (m/s): 230
Massa do Projétil (g): 10.24
Diâmetro do Projétil (mm): 8.82
Coeficiente de Arrasto: 0.235800
Rotação do Projétil (1-Dextrogiro, 2-Levogiro): 1
Twist rate do Projétil: 1:16.5
Comprimento do Projétil (mm): 16.6
Velocidade do Vento (km/h): 10
Direção do Vento (graus): 100
Calcular Coriolis Separadamente: 1
```


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

O programa foi executado, considerando o arqvuio **edificacao** com número de edificações igual a **0**.
```
# Dados da Edificação
# Substitua os valores pelos dados da edificação de interesse.
# Cada linha contém um rótulo descritivo e um valor. Altere apenas o valor.

Numero de edificacoes: 0
Latitude (graus): 0
Longitude (graus): 0
Altura (metros): 0
```

### Resultado do programa

>------------------------------------------------------------------------------------------------------------
>
>Leitura e ordenamento das edificações a partir da mais próxima da impactação.
>
>------------------------------------------------------------------------------------------------------------
>
>Sem edificações e considerando o arrasto, os cálculos terminaram com os seguintes valores:
>
>Tempo de deslocamento total do projétil: 3.0 segundos.
>
>Downrange Total = 526.3 m.
>
>Ângulo θ (inicial) do disparo = 13.5º.
>
>Ângulo Az (inicial) = 182.6º.
>
>Coordenadas Georgráficas do disparo ao NMM (N/S, L/O): -8.122989, -34.898140.

Ao examinar as coordenadas de possível origem observou-se que havia diversas edificações próximas.


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

> Sugestão de Compilação

> gcc Trajetoria_Balistica.c projetil.c -o Trajetoria -lm
