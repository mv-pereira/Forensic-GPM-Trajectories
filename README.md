# Forensic-GPM-Trajectories
Estimativa de origem de disparos subsônicos baseado em Trajetórias de Massas Pontuais Generalizadas a partir de Projéteis Dinamicamente estáveis em toda a trajetória utilizando a linguagem de programação C (Sistema usado: Ubuntu 22.04.3 LTS, Kernel: 5.15.0-35-generic.).

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
- **Espessura** do vitral impactado (em mm), **0** se não houver;
  - Perda de energia no vitral realizado a partir de ajuste quadrático. Conferir [Discussão](https://github.com/mv-pereira/Forensic-GPM-Trajectories/blob/186e238272a6d0f465e6bf1fb8d2b4c36012155d/9mm%20Balistics/Redu%C3%A7%C3%A3o%20de%20Velocidade.ipynb).
- **Azimute** estimado do disparo;
  - Geralmente próximo ao Azimute da impactação;
- **Twist rate** do Projétil;
  - Baseado no cano do armamento;
  - http://www.taurususa.com
- **Comprimento** do projétil;
- Opção de exibir a contribuição do efeito Coriolis em separado;
- **Temperatura** da Região (ºC);
  - https://earth.nullschool.net/
- **Altura média** da Cidade ou Região;
  - https://pt-br.topographic-map.com/map-6r951/Recife/?center=-8.07657%2C-34.96124&zoom=12

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
>Downrange Total = 524.5 m.
>
>Ângulo θ (inicial) do disparo = 13.6º.
>
>Ângulo Az (inicial) = 183º.
>
>Coordenadas Georgráficas do disparo ao NMM (N/S, L/O): -8.123007, -34.898122.

Ao examinar as coordenadas de possível origem observou-se que havia diversas edificações próximas.

![disparo_01](https://github.com/mv-pereira/Forensic-GPM-Trajectories/assets/86118560/7098c0e0-0aff-46e2-9502-1ae93acceb33)

Ao editar o arquivo **edificacao** com número de edificações igual a **5**:
```
# Dados da Edificação
# Substitua os valores pelos dados da edificação de interesse.
# Cada linha contém um rótulo descritivo e um valor. Altere apenas o valor.

Numero de edificacoes: 5
Latitude (graus): -8.123596,-8.124629,-8.125496,-8.126696,-8.127302
Longitude (graus): -34.898176,-34.898046,-34.898242,-34.898255,-34.898351
Altura (metros): 60,50,55,70,51
```

O programa conclui que, ao invés da origem estimada, partiu, de fato, de uma das edificações indicadas, com o segunite aviso:
>------------------------------------------------------------------------------------------------------------
>
>Testando a edificação #2.
>
>Coordenadas: (N/S, L/O): -8.124629, -34.898046, Altura 50.00 m, Distância até Impacto: 346.84 m.
>
>A distância lateral aproximada que o projétil partindo do solo (ou nível do mar) passa deste ponto vale: 19.62 m.
>
>------------------------------------------------------------------------------------------------------------
>
>O projétil provavelmente partiu desta edificação, no entanto, o projétil passa a mais de cinco metros (5 m) de distância do ponto fornecido para a edificação.
>
>Considere escolher um ponto para esta edificação mais próximo das coordenadas (N/S) -8.124620, -34.898224.
>
>Caso inexista, considere a remoção desta edificação e reinicie os cálculos.

Como de fato a edificação número 2 inserida distava muito da linha de tiro foi removida e o cálculo foi reiniciado (**lembrando de atualizar para 4 o número de edificações e remover a latitude, longitude e altura correspondente**).

```
# Dados da Edificação
# Substitua os valores pelos dados da edificação de interesse.
# Cada linha contém um rótulo descritivo e um valor. Altere apenas o valor.

Numero de edificacoes: 4
Latitude (graus): -8.123596,-8.125496,-8.126696,-8.127302
Longitude (graus): -34.898176,-34.898242,-34.898255,-34.898351
Altura (metros): 60,55,70,51
```

Após remoção e recálculo, houve o seguinte resultado:
```
------------------------------------------------------------------------------------------------------------
Testando a edificação #2.
Coordenadas: (N/S, L/O): -8.125496, -34.898242, Altura 55.00 m, Distância até Impacto: 248.82 m.
A distância lateral aproximada que o projétil partindo do solo (ou nível do mar) passa deste ponto vale: 1.94 m.
------------------------------------------------------------------------------------------------------------
O projétil passou por cima desta edificação.



------------------------------------------------------------------------------------------------------------
Testando a edificação #1.
Coordenadas: (N/S, L/O): -8.123596, -34.898176, Altura 60.00 m, Distância até Impacto: 460.40 m.
A distância lateral aproximada que o projétil partindo do solo (ou nível do mar) passa deste ponto vale: 2.35 m.
------------------------------------------------------------------------------------------------------------
O projétil provavelmente partiu desta edificação, visto que o ponto do edifício fornecido dista 2.35 m da linha de tiro.
Caso infira que o projétil não partiu dessa edificação, remova-a da lista de edificações e reinicie os cálculos.


Considerando um sistema com arrasto, os cálculos terminaram com os seguintes valores:
Foram efetuados 30 iterações para cálculos de trajetória.
Downrange Total = 477.594 m.
Altura de impactação = 89.000 m.
Desvios para direita = 1.440 m.
Ângulo θ (inicial) do disparo = 12.25º.
Ângulo ϕ (ao final da simulação) = 4.10º.
Azimute inicial do disparo = 183.00º.

O efeito de Spin Drift foi responsável por deslocar o projétil aproximadamente, 1.33 m para direta.

O efeito Coriolis foi responsável por deslocar o projétil:
Eixo x: 0.047342 cm;
Eixo y: -0.315902 cm;
Eixo z: -2.078656 cm;


------------------------------------------------------------------------------------------------------------
O projétil partiu, aproximadamente, das coordenadas (N/S, L/O): -8.123597, -34.898155, a uma altura de 15.07 m, partindo de uma edificação.

O projétil levou cerca de 2.7 segundos do momento do disparo à impactação.
Possuía velocidade final de 147.00 m/s e energia cinética de 110.63 J.
O disparo foi efetuado com o cano apontado, aproximadamente, a 15 m acima do alvo, ou 1.70º acima da impactação.
------------------------------------------------------------------------------------------------------------

Energia cinética de alguns projéteis para comparação:
Calibre 	Energia (J)
.25 AUTO	87     
.32 AUTO	175    
.380 AUTO EXPO  259
.38 SPL CHOG    271
9x19mm (124gr)  459
.40 EXPO Gold   568
.357 Mag        724
.454 Casull     2531
```

Após testar as edificações, da mais próxima da impactação à mais longe, os cálculos indicaram ter partido de uma delas.

É possível saber a contribuição para o efeito Coriolis e SpinDrift.


Assim, procurou-se as coordenadas e altura fornecida pelo programa:

![disparo_02](https://github.com/mv-pereira/Forensic-GPM-Trajectories/assets/86118560/0b0927ec-cbac-4c65-8182-f8f96f019f37)


A região que ocorreu o disparo é dada na posição do mouse.

O programa fornece, por fim, a estimativa de duração do tempo de deslocamento do projétil, Velocidade Final e Energia Cinética para comparação com valores nominais de outros calibres.

| Resultado                     | Valores        |
|-------------------------------|----------------|
| Tempo de Trajetória           | 2.7 segundos   |
| Velocidade Final              | 147.00 m/s     |
| Energia Cinética              | 110.63 J       |

> É importante destacar que, caso haja ainda outras edificações próximas à trajetória fornecida, é importante que sejam incluídas de modo ao programa identificar se o projétil passou acima das edificações ou partiu delas.

#### Fica a cargo do **usuário** a indicação precisa da posição e altura das edificações.


> Sugestão de Compilação

> gcc Trajetoria_Balistica.c projetil.c -o Trajetoria -lm
