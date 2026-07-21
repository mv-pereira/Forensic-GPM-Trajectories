# Estimador de Curvas de Arrasto Gx

Programa em C para estimar qual curva de arrasto de referência representa melhor o comportamento de um projétil a partir de velocidades conhecidas em diferentes distâncias.

O programa recebe dados normalmente fornecidos pelo fabricante, simula a desaceleração do projétil para diferentes curvas Gx, ajusta um fator multiplicativo de arrasto e seleciona a combinação que produz o menor erro em relação às velocidades informadas.

## Objetivo

Fabricantes de munição frequentemente publicam:

- velocidade inicial;
- velocidade em uma ou mais distâncias;
- massa do projétil;
- diâmetro do projétil;
- coeficiente balístico, geralmente referido à curva G1.

Entretanto, um único coeficiente balístico nem sempre descreve adequadamente o comportamento do projétil em toda a faixa de velocidades. O arrasto varia com o número de Mach, principalmente na região transônica, próxima de Mach 1.

Este projeto utiliza curvas tabuladas de coeficiente de arrasto $C_d$ em função do número de Mach para procurar uma representação compatível com os dados de velocidade observados.

## Funcionamento

O programa executa, em linhas gerais, as seguintes etapas:

1. Solicita os dados físicos do projétil e as velocidades de referência.
2. Calcula a velocidade do som para a temperatura informada.
3. Simula a trajetória e a perda de velocidade para cada curva disponível.
4. Ajusta um fator multiplicativo aplicado ao $C_d$ de cada curva.
5. Calcula o erro entre as velocidades simuladas e as velocidades informadas.
6. Seleciona a curva e o fator que apresentam o menor erro RMS.
7. Exporta uma tabela CSV com $C_d$ em função da velocidade para a melhor curva encontrada.

## Curvas disponíveis

Atualmente, o arquivo `curvas_gx.h` contém as seguintes curvas:

- MCG1;
- MCG2;
- MCG5;
- MCG6;
- MCG7;
- MCG8;
- MCGI;
- RA4.

Cada curva é armazenada como uma tabela de pares:

$$
(M, C_d)
$$

em que:

- $M$ é o número de Mach;
- $C_d$ é o coeficiente de arrasto correspondente.

Quando o número de Mach calculado está entre dois pontos da tabela, o programa determina o valor de $C_d$ por interpolação linear:

$$
C_d(M) = C_{d,1} + \frac{M-M_1}{M_2-M_1}\left(C_{d,2}-C_{d,1}\right)
$$

## Dados de entrada

O programa solicita:

- velocidade inicial em $0\,\text{m}$, em metros por segundo;
- primeira distância de referência, em metros;
- velocidade medida ou informada nessa distância, em metros por segundo;
- segunda distância de referência, em metros;
- velocidade medida ou informada nessa distância, em metros por segundo;
- temperatura ambiente, em graus Celsius;
- massa do projétil, em gramas;
- diâmetro do projétil, em milímetros.

Embora sejam usados dois pontos intermediários de velocidade, as distâncias não precisam ser necessariamente 50 m e 100 m. O usuário pode informar outros valores, desde que sejam coerentes, crescentes e estejam dentro da faixa em que os dados do fabricante sejam confiáveis.

## Número de Mach

O número de Mach é calculado pela relação:

$$
M = \frac{v_{\mathrm{rel}}}{a}
$$

em que:

- $v_{\mathrm{rel}}$ é a velocidade do projétil em relação ao ar;
- $a$ é a velocidade local do som.

A velocidade do som é aproximada no programa por:

$$
a = 331{,}3 + 0{,}606T
$$

em que $T$ é a temperatura em graus Celsius e $a$ é obtida em metros por segundo.

## Ajuste da curva

Para cada curva de referência, o programa aplica um fator multiplicativo ao coeficiente de arrasto tabulado:

$$
C_{d,\mathrm{ajustado}}(M) = i\,C_{d,\mathrm{base}}(M)
$$

em que:

- $C_{d,\mathrm{base}}(M)$ é o valor interpolado da curva de referência;
- $i$ é o fator de ajuste determinado numericamente;
- $C_{d,\mathrm{ajustado}}(M)$ é o coeficiente utilizado na simulação.

No código, esse valor aparece como `fator_ajuste`.

O programa testa fatores dentro de uma faixa predefinida e mantém aquele que produz o menor erro RMS para cada curva.

É importante observar que esse fator é semelhante ao conceito de fator de forma da balística exterior, mas deve ser interpretado como um parâmetro de ajuste do modelo implementado. Ele não deve ser tratado automaticamente como um fator de forma normalizado obtido por ensaio aerodinâmico.

## Cálculo do erro

Para cada curva e fator são calculadas as diferenças entre as velocidades simuladas e as velocidades informadas:

$$
e_1 = v_{1,\mathrm{calculada}} - v_{1,\mathrm{referência}}
$$

$$
e_2 = v_{2,\mathrm{calculada}} - v_{2,\mathrm{referência}}
$$

O erro RMS é:

$$
E_{\mathrm{RMS}} = \sqrt{\frac{e_1^2+e_2^2}{2}}
$$

A curva escolhida é aquela que apresenta o menor valor de $E_{\mathrm{RMS}}$.

O erro é expresso em metros por segundo.

## Cálculo do coeficiente balístico a partir do fator G1

Embora o programa informe a melhor curva a ser utilizada, caso o usuário decida comparar o resultado prático do programa com o Coeficiente Balístico informado pelo fabricante, é possível estimar a partir da curva G1, utiliza-se a densidade seccional do projétil e o fator de forma $i_G1$ fornecido pelo programa:

$$
BC_{G1} = \frac{SD}{i_{G1}}
$$

Na convenção tradicional, a densidade seccional utiliza massa em libras e diâmetro em polegadas. Quando os dados estão em gramas e milímetros, aplicam-se os seguintes fatores de conversão:

$$
m_{\mathrm{lb}} = \frac{m_{\mathrm{g}}}{453{,}59237}
$$

$$
d_{\mathrm{in}} = \frac{d_{\mathrm{mm}}}{25{,}4}
$$

Substituindo essas conversões na fórmula da densidade seccional:

$$
SD =
\frac{
m_{\mathrm{g}} / 453{,}59237
}{
\left(
d_{\mathrm{mm}} / 25{,}4
\right)^2
}
$$

De forma simplificada:

$$
SD =
1{,}42233433
\frac{
m_{\mathrm{g}}
}{
d_{\mathrm{mm}}^2
}
$$

Assim, o coeficiente balístico G1 pode ser estimado diretamente por:

$$
BC_{G1} =
\frac{
1{,}42233433 \, m_{\mathrm{g}}
}{
d_{\mathrm{mm}}^2 \ i_{G1}
}
$$

em que:

- $m_g$ é a massa do projétil em gramas;
- $d_mm$ é o diâmetro do projétil em milímetros;
- $i_{G1}$ é o fator correspondente à curva G1;
- $BC_{G1}$ é o coeficiente balístico na convenção G1.

O resultado segue a convenção normalmente utilizada para valores comerciais de coeficiente balístico G1.
