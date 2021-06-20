# Forensic-GPM-Trajectories
Estimativa de origem de disparos subsônicos baseado em Trajetórias de Massas Pontuais Generalizadas a partir de Projéteis Dinamicamente estáveis em toda a trajetória utilizando a linguagem de programação C.


# Lançamento Horizontal Sem Arrasto

<img src="https://user-images.githubusercontent.com/86118560/122644198-6c21a100-d0ea-11eb-82ff-6556f6ccde0f.png" alt="1" width="650" height="300">
Qual a distância d para que um projétil de velocidade inicial v0 colida a uma altura y0 do solo, formando um ângulo ϕ com a horizontal?
	
Através da segunda Lei de Newton, podemos escrever a equação das forças atuantes na partícula:
![image](https://user-images.githubusercontent.com/86118560/122644182-57dda400-d0ea-11eb-9083-3f68faf2da6a.png)

Escrevendo a força atuante em cada eixo (existe apenas a força Peso):
  ![image](https://user-images.githubusercontent.com/86118560/122644188-61ffa280-d0ea-11eb-8c55-c23276f647e7.png)
, considerando a aceleração da gravidade no sentido.![image](https://user-images.githubusercontent.com/86118560/122644193-662bc000-d0ea-11eb-976e-d1c0994ba691.png)


Considerando as seguintes condições iniciais:
![image](https://user-images.githubusercontent.com/86118560/122644219-80fe3480-d0ea-11eb-8cee-6062bbb4f746.png)
 e , ![image](https://user-images.githubusercontent.com/86118560/122644222-82c7f800-d0ea-11eb-88dc-98ba4cc771b6.png)
onde ![image](https://user-images.githubusercontent.com/86118560/122644224-852a5200-d0ea-11eb-9aa2-683a41e55142.png)
.

A integração e devida aplicação dessas condições iniciais, leva às seguintes equações de movimento, na qual a trajetória dependente do tempo t pode ser descrita em termos das equações paramétricas:
![image](https://user-images.githubusercontent.com/86118560/122644333-126da680-d0eb-11eb-935f-378368d151e3.png)
, onde vo é a velocidade inicial, g, aceleração da gravidade no local e θ, o ângulo inicial do disparo em relação à horizontal.

A tangente do ângulo em relação à horizontal (ϕ) ao longo da trajetória pode ser dado por:
![image](https://user-images.githubusercontent.com/86118560/122673789-7f934180-d1a8-11eb-834c-4b649071b91b.png)


![image](https://user-images.githubusercontent.com/86118560/122673819-92a61180-d1a8-11eb-831b-abb5bb95a49d.png)


Como ϕ varia no tempo, podemos escolher o tempo t0 tal que y0 seja a altura de interesse conhecida.

Seja y(t0) = y0 em [II]:
![image](https://user-images.githubusercontent.com/86118560/122644342-1bf70e80-d0eb-11eb-8973-d703ce558588.png)	
 
<img src="https://user-images.githubusercontent.com/86118560/122644344-1f8a9580-d0eb-11eb-9ca0-372839d045f1.png" alt="2" width="600" height="300">


É possível identificar os termos a = (-1/2)g, b = v0sin(θ) e c = -y0. Com isso, a equação [IV] possui duas raízes dadas por :
![image](https://user-images.githubusercontent.com/86118560/122644353-29ac9400-d0eb-11eb-8aef-df5d8e6de200.png)

![image](https://user-images.githubusercontent.com/86118560/122644354-2b765780-d0eb-11eb-9b58-d3e399331d8b.png)





Substituindo [V] em [III] obteremos o ângulo ϕ na altura y0:
![image](https://user-images.githubusercontent.com/86118560/122644359-2fa27500-d0eb-11eb-97c1-dff93bd538ab.png)

![image](https://user-images.githubusercontent.com/86118560/122644363-30d3a200-d0eb-11eb-9882-291471498bba.png)


Precisamos, então, isolar θ, de modo que outras equações em x, fiquem em função de ϕ e outras constantes conhecidas. De [VI] e:
![image](https://user-images.githubusercontent.com/86118560/122644370-3af5a080-d0eb-11eb-980c-f5b42a0dcb5a.png)

![image](https://user-images.githubusercontent.com/86118560/122644371-3df09100-d0eb-11eb-8b2c-8d25a2a48bb8.png)

![image](https://user-images.githubusercontent.com/86118560/122644373-3f21be00-d0eb-11eb-96f2-3c9f1fb6d9c0.png)


![image](https://user-images.githubusercontent.com/86118560/122644374-41841800-d0eb-11eb-9caa-81a58ad5ed0c.png)
, no entanto, qual o significado dos sinais do argumento do arc sec (x)?

<img src="https://user-images.githubusercontent.com/86118560/122644399-5f517d00-d0eb-11eb-9ca8-a8edcfee33fc.png" alt="1" width="400" height="150">

Pela interpretação física do problema, podemos escolher apenas o argumento positivo da função arcsec para futuras substituições nas equações.



Podemos, agora, manipular [I] para isolar t e substituir em [II] tendo em vista obter y(x).

  ![image](https://user-images.githubusercontent.com/86118560/122644405-64aec780-d0eb-11eb-9d0e-fb4fc2df31ed.png)
[VIII] em [II]:	

![image](https://user-images.githubusercontent.com/86118560/122644409-65dff480-d0eb-11eb-9e21-d0c4bb7678ce.png)

![image](https://user-images.githubusercontent.com/86118560/122644411-67112180-d0eb-11eb-8afe-75c6905c4c62.png)


Calculando as raízes de x em [IX] obtemos:
![image](https://user-images.githubusercontent.com/86118560/122644415-6aa4a880-d0eb-11eb-9baf-a5806600766e.png)


Podemos, por fim, substituir [VII] em [X], utilizando [XI] para obter x em função das variáveis conhecidas v0, ϕ, y0, g.
Utilizando as identidades:
![image](https://user-images.githubusercontent.com/86118560/122644417-6bd5d580-d0eb-11eb-865b-507be14a533f.png)


Obtemos:
![image](https://user-images.githubusercontent.com/86118560/122644419-6e382f80-d0eb-11eb-9851-a8c476a602a1.png)
, assumindo que g, v0, y0 e ϕ são reais.

Onde x+ e x- representam as distâncias à origem do ponto cujo módulo do ângulo vale ϕ. Em outras palavras, o projétil, em sua trajetória parabólica, inicia o lançamento com ângulo θ, que decresce até atingir ϕ no trecho ascendente da trajetória a uma distância x- da origem, chegando a 0 no ponto mais alto (vy = 0 [III]), que decresce a -ϕ no trecho descendente a uma distância x+ da origem.

![image](https://user-images.githubusercontent.com/86118560/122674432-88d1dd80-d1ab-11eb-955a-2cda235eefcf.png)

# Lançamento Horizontal com Arrasto

Com o intuito de resolver o mesmo tipo de problema, ou seja, encontrar a distância d para que um projétil de velocidade inicial v0 colida a uma altura H do solo, formando um ângulo ϕ com a horizontal, vamos introduzir, na segunda Lei de Newton, o termo conhecido como Equação do Arrasto ou, simplesmente, Arrasto (Fd): ![image](https://user-images.githubusercontent.com/86118560/122674456-b4ed5e80-d1ab-11eb-8b7c-b3f8c38bb451.png)

onde:
![image](https://user-images.githubusercontent.com/86118560/122674477-cafb1f00-d1ab-11eb-82bd-2897835b08a4.png)


Fd é a Força de Arrasto;


Cd é o Coeficiente de Arrasto;


ρ é a Densidade do Ar;


A é a Área de Seção transversal do objeto;


V é a velocidade do objeto.


Com isso: ![image](https://user-images.githubusercontent.com/86118560/122674482-cf273c80-d1ab-11eb-81d8-d9edf1ae57c7.png)
, onde k = (1/2)CdρA/m;


A introdução dessa força é o primeiro passo para o estabelecimento do método conhecido como trajetória de massa pontual (Point-mass Trajectory), pois desconsidera diversos efeitos que ocorrem com um projétil tridimensional. Em primeiro momento vamos desconsiderar o efeito do vento e efeito Coriolis, a serem implementados em seguida.


É importante observar, contudo, que a consideração de um coeficiente de arrasto constante é válido, apenas, em velocidade subsônicas.

Para seguir com a solução de [II], podemos, de maneira semelhante à resolução do problema sem arrasto, decompor essa força em duas forças componentes sobre as direções perpendiculares x e y: ![image](https://user-images.githubusercontent.com/86118560/122674484-d4848700-d1ab-11eb-8f66-74a0891c82cf.png) , que, em termos de vx e vy, podem ser escritas como:
![image](https://user-images.githubusercontent.com/86118560/122674491-dbab9500-d1ab-11eb-99e5-3a777ad68ef9.png)
 pois ![image](https://user-images.githubusercontent.com/86118560/122674492-de0def00-d1ab-11eb-8281-3a26cf4d778b.png)
 e ![image](https://user-images.githubusercontent.com/86118560/122674493-dfd7b280-d1ab-11eb-9f13-0af8ecb73c29.png)
.

Considerando ![image](https://user-images.githubusercontent.com/86118560/122674505-ecf4a180-d1ab-11eb-88a7-9c05b96f89f4.png)
 e ![image](https://user-images.githubusercontent.com/86118560/122674510-eebe6500-d1ab-11eb-8c7e-d80fb3b5395c.png)
, vamos reescrever o conjunto de equações diferenciais ordinárias acopladas [IV] como: ![image](https://user-images.githubusercontent.com/86118560/122674513-f251ec00-d1ab-11eb-8e05-86b993dc7862.png)

, com as condições iniciais: ![image](https://user-images.githubusercontent.com/86118560/122674517-f54cdc80-d1ab-11eb-960e-4715a8a6a252.png)


Vamos redefinir as variáveis auxiliares: ![image](https://user-images.githubusercontent.com/86118560/122674535-00a00800-d1ac-11eb-898e-d89afb54d44b.png)
 e substituir em [IV]: ![image](https://user-images.githubusercontent.com/86118560/122674536-0269cb80-d1ac-11eb-85fb-10c1d09aef96.png)

, com condições iniciais [V].

Vamos aplicar o método de Runge-Kutta de 4a Ordem (RK4) no sistema [VII] de modo a encontrar o y de interesse. Note que precisamos, também, que tan ϕ < 0 → vy/vx < 0, ou seja, trajetória descendente. Caso a trajetória de impactação seja ascendente, a condição tan ϕ < 0 é desnecessária.

Esboçando os passos iniciais para melhor entendimento, vamos calcular: ![image](https://user-images.githubusercontent.com/86118560/122674545-085fac80-d1ac-11eb-9115-7b6548e6ce4a.png)


Onde h é um passo arbitrariamente escolhido e kα é a média ponderada das inclinações das retas tangentes nos pontos n, n+h/2 e n+h, conforme segue: ![image](https://user-images.githubusercontent.com/86118560/122674552-0b5a9d00-d1ac-11eb-80e4-4c5028c2d737.png)



com: ![image](https://user-images.githubusercontent.com/86118560/122674571-1f9e9a00-d1ac-11eb-991b-af7518d72ef1.png)

	
Resumindo: Os próximos valores (x(n+1), vx (n+1), ... ) são determinados pelos valores atuais (x(n), vx (n), ... ) somados com o produto do tamanho do intervalo (h) e uma inclinação estimada. A inclinação é uma média ponderada de inclinações: 
    • k1 é a inclinação no início do intervalo; 
    • k2 é a inclinação no ponto médio do intervalo, usando a inclinação k1 para determinar o valor de y no ponto tn + h/2 através do método de Euler; 
    • k3 é novamente a inclinação no ponto médio do intervalo, mas agora usando a inclinação k2 para determinar o valor de y; 
    • k4 é a inclinação no final do intervalo, com seu valor y determinado usando k3.

Vamos realizar um exemplo para as variáveis x e vx, calculando x(1) e vx(1) a partir de [V]:	
kx1 = wx(0,0,V0cosθ0,V0sinθ0), mas como wx = vx (de [VII]), teremos:
kx1 = V0cosθ0


kx2 = wx(h/2, V0cosθ0* h/2, V0cosθ0 + V0cosθ0*h/2, V0sinθ0 + V0cosθ0*h/2);
kx2 = V0cosθ0 + V0cosθ0*h/2 = V0cosθ0(1+ h/2), pois wx = vx;


kx3 = wx(t0 + h/2, x0 + kx2* h/2, vx0+ kx2*h/2, vy0+ kx2*h/2);
kx3 = wx(h/2 , 0 + V0cosθ0(1+ h/2) , V0cosθ0 + V0cosθ0(1+ h/2)(h/2), V0sinθ0 + V0cosθ0(1+ h/2));
kx3 = V0cosθ0 + V0cosθ0(1+ h/2)(h/2) = V0cosθ0( 1 + h/2 + h²/4);
 
 
kx4 = wx(t0 + h, x0 + kx3* h, vx0 + kx3*h, vy0 + kx3*h);
kx4 = wx( h , h* V0cosθ0( 1 + h/2 + h²/4) ,  V0cosθ0 + h* V0cosθ0( 1 + h/2 + h²/4),  V0sinθ0 + h* V0cosθ0( 1 + h/2 + h²/4));
kx4 = V0cosθ0 + h* V0cosθ0( 1 + h/2 + h²/4) = V0cosθ0(1+h+h²/2+h³/4);


Substituindo todos os valores em kx a partir de [IX]:
kx = (1/6)(kx1 + 2*kx2 + 2*kx3 + kx4)
kx = (1/6)V0cosθ0h(6+3h+h²+h³/4), logo, como x(n+1) = xn + kx*h;


x1 = 0 + h*(1/6)V0cosθ0h(6+3h+h²+h³/4).


Vamos calcular somente o primeiro termo de vx para exemplificação, pois o cálculo segue da mesma forma:
	kvx1 = wvx(t0,x0,vx0,vy0);
	kvx1 = wvx(0,0,V0cosθ0,V0sinθ0); de [VII]:
	kvx1 = -(k/m)(V0cosθ0)² sqrt ((V0cosθ0)² + (V0sinθ0)²) = -(k/m)(V0³cos²θ0)


kvx2 = wvx(t0 + h/2, x0 + kvx1* h/2, vx0 + kvx1* h/2, vy0 + kvx1* h/2);
kvx2 = wvx( h/2, -(kh/2m)(V0³cos²θ), V0cosθ0 - (kh/2m)(V0³cos²θ), V0sinθ0 - (kh/2m)(V0³cos²θ));
kvx2 = -(k/m)(V0cosθ0 – (kh/2m)(V0³cos²θ0))² sqrt [(V0cosθ0 – (kh/2m)(V0³cos²θ0))² + (V0sinθ0 - (kh/2m)(V0³cos²θ0))² ] = …
kvx2 = …


kvx3 = …


kvx4 = …


kvx = (1/6)(kvx1 + 2*kvx2 + 2*kvx3 + kvx4), e então: vx(1) = vx(0) + kvx*h;



Portanto, a medida que o tempo t aumenta com o passo h, não somente as equações em vx(n+1) e vy(n+1) modificam-se, como também as de x(n+1) e y(n+1), visto que o valor (n+1) de x e y vão depender dos valores anteriores de vx e vy, que vão se modificando.

Em outras palavras, na equação para x: x(n+1) = xn + kx * h mesmo o primeiro termo xn não dependa diretamente de v, a inclinação da reta k depende, pois wx(t,x,vx,vy) = vx.

Esse cálculo utilizando RK4 iniciou-se de um problema de valor inicial, ou seja, de antemão, temos os valores [V], o que não representa a realidade. A altura H e o ângulo medido ϕ durante a atividade pericial representam o estado final do problema, não o inicial. Na prática, tampouco utilizaremos RK4 para um problema de valor final, pois não sabemos a velocidade final do projétil.


Sugiro, então, utilizarmos as soluções encontradas para o caso sem arrasto como ponto de partida, pois podemos estimar a condição de ângulo inicial do disparo θ0. Como a distância encontrada para o caso sem arrasto será maior que para o caso com arrasto, é possível transladarmos toda a curva até que o valor de y no ponto (xi,yi) coincida com a altura H do prédio, durante a trajetória descendente do projétil.


Transladando ou não a curva, caso o projétil com arrasto chegue no ponto final com uma angulação menor que ϕ, ou seja, ϕf < ϕ (lembrando que, em trajetória descendente, ϕf e ϕ < 0), incrementamos o ângulo incial de disparo θ0 com um fator e recalculamos toda a trajetória. Por outro lado, se ao fim ϕf > ϕ, reduzimos  θ0 com um fato e recalculamos toda a trajetória.


O critério de parada do cálculo é que ![image](https://user-images.githubusercontent.com/86118560/122674610-4e1c7500-d1ac-11eb-90d4-e1afafdb3a1f.png), sendo um δϕ um valor arbitrário aceitável.
![image](https://user-images.githubusercontent.com/86118560/122674757-f3374d80-d1ac-11eb-987e-2abdd020d8b8.png)


Sobre o método RK 4:
https://www.youtube.com/watch?v=l8sOUHnPkgw

Cálculo de coeficientes de arrasto:
https://www.jbmballistics.com/cgi-bin/jbmdrag-5.1.cgi


# Trajetória de Massa Pontual

Para respondermos ao questionamento incial, mas com três dimensões, é necessário medir mais um parâmetro na região do impacto. Nesse caso, precisamos medir não somente o ângulo ϕ, inclinação com o plano horizontal, como também a inclinação com o plano vertical, tomando como ângulo zero, o Norte (Azimute). Chamaremos tal ângulo γ. Em uma vista superior:


![image](https://user-images.githubusercontent.com/86118560/122674995-e7985680-d1ad-11eb-906c-9bb1c2311b55.png)


Adotar um sistema de referencial fixo na superfície terrestre leva ao aparecimento de forças fictícias, dentre elas a força de Coriolis. Vamos introduzir os termos correspondentes ao efeito do vento e efeito Coriolis no sistema já tratado anteriormente: ![image](https://user-images.githubusercontent.com/86118560/122675011-f2eb8200-d1ad-11eb-8e8b-405c9a57c173.png)


onde:
![image](https://user-images.githubusercontent.com/86118560/122675013-f4b54580-d1ad-11eb-9252-630f3a517205.png), representa o somatório de todas as forças aerodinâmicas;


![image](https://user-images.githubusercontent.com/86118560/122675018-f7179f80-d1ad-11eb-9e74-402e8419cbe6.png), a força Peso;


![image](https://user-images.githubusercontent.com/86118560/122675021-f8e16300-d1ad-11eb-8515-77ca37186b7c.png), o vetor aceleração de Coriolis, devido a rotação da terra.



Vamos adicionar as componentes da velocidade do ventona trajetória da massa pontual, de modo que a componente do vento tem valor positivo quando aponta na direção positiva de cada um dos eixos.


![image](https://user-images.githubusercontent.com/86118560/122675061-2e864c00-d1ae-11eb-812f-40d65680daba.png)


A adição da aceleração Coriolis será realizada a adição do vetor ![image](https://user-images.githubusercontent.com/86118560/122675066-334b0000-d1ae-11eb-86a1-5013661ff99a.png).


Para uma massa pontual, o vetor velocidade  na expressão do arrasto precisa ser	substituído por ![image](https://user-images.githubusercontent.com/86118560/122675076-3e059500-d1ae-11eb-8bb5-df2188cf5fa6.png), pois o arrasto aerodinâmico depende da velocidade relativa ao fluxo de ar e não ao solo.


![image](https://user-images.githubusercontent.com/86118560/122675101-54135580-d1ae-11eb-863f-fe73b1c8e319.png)


Com isso, a expresão [I], em um sistema sem o efeito Coriolis, se torna: ![image](https://user-images.githubusercontent.com/86118560/122675113-5d042700-d1ae-11eb-8345-10fcc4d3874f.png), onde k = (1/2)CdρA/m;


Com o termo:![image](https://user-images.githubusercontent.com/86118560/122675118-5fff1780-d1ae-11eb-95be-5a79f57a0bed.png).


É importante entender o significado de de wx, wy e wz. Nas diversas cartas ou mapas de vento, sua direção é dada em graus, medido em sentido horário em relação ao Norte. No nosso caso, desde o problema sem arrasto, denominamos x o eixo responsável pelo alcance, ou seja, vx, vy e vz estão no  sistema de referencial no qual o eixo x é o eixo de deslocamento principal do projétil.


Assim, wx, wy e wz são as componentes da velocidade do vento no nosso referencial de interesse e não em relação ao Norte. Podemos fazer essa transformação através da aplicação de uma matriz mudança de base Norte-Leste para x-z.

Além disso, as direções do vento são dadas pelo lugar de onde o vento vem, não para onde sopram. No nosso caso então, é necessário adicionar 180o (ver Exemplo).


Vamos definir que a direção na qual o disparo foi realizado será chamada de Azimute inicial (AZ0). Com isso, nosso eixo x está rotacionado com AZ0 em relação ao Norte.


![image](https://user-images.githubusercontent.com/86118560/122675177-a05e9580-d1ae-11eb-8aca-9305e6ca7e47.png)


A mudança de base de vN e vE do eixo N-E para os eixos x e z, em termos de AZ0, pode ser definida por: ![image](https://user-images.githubusercontent.com/86118560/122675181-a5234980-d1ae-11eb-8543-c5436bb0e271.png)


O vetor aceleração de Coriolis, sem componentes da velocidade do vento, pode ser escrito como:


![image](https://user-images.githubusercontent.com/86118560/122675186-ac4a5780-d1ae-11eb-8576-0c8628b46288.png), onde L é a latitude do local do disparo (positivo para hemisfério Norte e negativo para Sul e AZ é o azimute do disparo, medido em sentido horário em relação ao Norte.

![image](https://user-images.githubusercontent.com/86118560/122675471-f7b13580-d1af-11eb-983f-03c5e5538b65.png)


Vamos, agora, adicionar as componentes da velocidade do vento em [III]:![image](https://user-images.githubusercontent.com/86118560/122675473-fa138f80-d1af-11eb-9aff-2116177671f1.png)


	
Unindo o resultado obtido em [II] com a aceleração de Coriolis [IV] obtemos, em notação vetorial:![image](https://user-images.githubusercontent.com/86118560/122675475-fc75e980-d1af-11eb-9cc2-240be176c255.png)
![image](https://user-images.githubusercontent.com/86118560/122675478-fe3fad00-d1af-11eb-9658-3f08530d8d79.png)


Definindo as variáveis auxiliares ζ para o cálculo do método RK 4:

![image](https://user-images.githubusercontent.com/86118560/122675484-05ff5180-d1b0-11eb-847d-b60650862899.png)


Com as condições iniciais:

![image](https://user-images.githubusercontent.com/86118560/122675491-0a2b6f00-d1b0-11eb-850a-dc121c88fbb6.png)



Onde vN e vE são as componentes do vetor velocidade do vento nos eixos Norte-Leste, para poder realizar a mudança de base. É importante observar que as condições iniciais podem possuir quaisquer outros valores.


É importante lembrar que o disparo, convencionalmente, ocorre paralelamente a x, e, nesse caso, a variação do eixo z indicará a deriva que ocorrerá por efeito do vento ou pelo efeito Coriolis, considerando que esses efeitos podem inserir variações também nos valores de x e ou y.


Considerações sobre o método empregado no programa


O critério adotado para correção do AZ0 compara a inclinação lateral do projétil na impactação (em relação ao eixo x) com o ângulo medido γ. Como o ângulo γ é medido em relação ao Norte, para comparação precisamos somar o valor da inclinação lateral com AZ0 e asism, essa inclinação lateral da impactação estará medida em relação ao Norte.



O critério de parada fica: ![image](https://user-images.githubusercontent.com/86118560/122675500-16afc780-d1b0-11eb-9d71-5066dcf6b383.png) , onde δγ é um valor arbritariamente pequeno. Que é um critério semelhante a ![image](https://user-images.githubusercontent.com/86118560/122675503-1a434e80-d1b0-11eb-91f1-eecd80900d62.png), utilizado na inclinação do cálculo do arrasto.


![image](https://user-images.githubusercontent.com/86118560/122675514-27603d80-d1b0-11eb-8e6b-5b77ae8fc50e.png)


O foco do programa é variar AZ0 até que o valor final da simulação (IL+AZ0) seja próximo ao valor medido γ.


Podemos aplicar a matriz inversa de mudança de base para obter os pontos em termos do Norte e Leste, e não do eixo x: ![image](https://user-images.githubusercontent.com/86118560/122675519-321ad280-d1b0-11eb-9444-418e0bc70a58.png), que pode ser entendida como a mesma mudança de base anterior através de um ângulo -AZ0.


O foco do programa, no entanto, são os resultados finais. A partir da última mudança de base, trajetória gravada será a baseada nas coordenadas geográficas Norte e Leste, mas o valor de x exibido na tela ao final do programa significa o alcance máximo do projétil (Downrange); já o valor de y significa a altura da impactação; z, o desvio lateral provocado pelo vento e efeito Coriolis (se z<0 significa que o projétil desviou para a esquerda em relação ao eixo x, se z>0, o oposto); θ, o ângulo de disparo em relação ao solo; Azimute: Ângulo azimute inicial do disparo.


![image](https://user-images.githubusercontent.com/86118560/122675536-44950c00-d1b0-11eb-9608-1c12ad443489.png)


Sobre o Efeito Coriolis:
https://www.youtube.com/watch?v=HM_p9gM1MUo

Sobre Aceleração Coriolis:
http://www.coastalwiki.org/wiki/Coriolis_acceleration

Sobre Rotação de Corpos Rígidos:
John Robert Taylor (2004). Classical Mechanics. Sausalito CA: University Science Books. Ch 9.

Sobre Trajetória de Massa Pontual:
Robert L. McCoy. Modern Exterior Ballistics. Schiffer publishing 1999.

Fontes de informação dos Ventos:
pt.windfinder.com
https://earth.nullschool.net/
www.climatempo.com.br


Exemplo:
	Vamos supor que um disparo atingiu o 11o andar de um edifício no final da Avenida Caxangá, bairro da Madalena, Recife/PE.
	Ao medir as informações, a equipe obteve os seguintes valores:
y0 = 36 m; ϕ = 8º; γ = 20º;

Velocidade do Projétil .38 = 305 m/s;

Cd (.38) = 0.2;

Massa (.38) = 8 g;

Diâmetro (.38) = 8.82 mm;

Velocidade do vento, direção e Latitude são dadas pela carta dos ventos.


<img src="https://user-images.githubusercontent.com/86118560/122675558-5676af00-d1b0-11eb-9eeb-a7482cd5cb26.png" alt="1" width="711" height="400">

Figura 14: Valores da velocidade do vento em Recife em 31/03/2020 as 15:00.
Fonte: https://pt.windfinder.com/#16/-8.0570/-34.9108/report

Argumentos no Gnuplot:
set xrange [0:900]  
set yrange [900:0]
set zrange [0:100] 
set parametric 
set urange [0:1000]
set vrange [0:50]
set style line 1 lw 3
splot "data" u 2:4:3 with lines ls 1, u,0,v with lines

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

O projetil é destrogiro ou levogiro? 1 - Destrogiro.    2 - Levogiro. 
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



<img src="https://user-images.githubusercontent.com/86118560/122675566-6098ad80-d1b0-11eb-9cd1-c935540e481e.png" alt="1" width="711" height="400">


<img src="https://user-images.githubusercontent.com/86118560/122675568-62fb0780-d1b0-11eb-942b-951683c8bb39.png" alt="1" width="711" height="400">
