# Programação Linear

## Alguns conceitos de Programação Linear:

- Simplex
- Pontos Interiores


## To do list

* [x] Transformar em classe
* [x] Adicionar Método das duas Fases
* [x] Adicionar Método Big M
* [] Adicionar Pontos Interiores


## Entrada

* Padronizando a maneira de resolver, o input deverá ser da seguinte forma:
    * Queremos min f(x) = c^T*x sujeita a
    * Ax = b
    * x>=0

* Logo definimos essa a forma padrão do PL que o algoritmo de simplex resolve. Para restrições Ax >= B devemos adicionar variáveis de folga para que fique da forma padrão.

### Exemplo (transformar para forma padrão)

Suponha que queremos:

min x1 + x2, sujeita a :

* x1 + x2 <= 4
* 2x1 + x2 <= 2
* x1,x2 >=0

Adicionamos variáveis de folgas x3,x4 >=0 e transformamos o PL em:

min x1 + x2 + 0x3 + 0x4, sujeita a :

* x1 + x2 + x3 = 4
* 2x1 + x2 + x4 = 8
* x1,x2,x3,x4 >=0


Obs: Caso as restricoes fossem >= subtraímos as variáveis de folgas (que continuam sendo não negativas)

