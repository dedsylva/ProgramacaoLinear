import numpy as np

#Temos o problema

#min f(x) = c^Tx sujeita a 
#Ax = b
#Primeiro passo: Particao basica --> A = [B N]
#Vamos primeiro assumir que ja temos as particoes, e 
#depois com as 2 fases e o big M a gente consegue elas



# exemplo
def f(x):
    return -x[0] -3*x[1]

pare = False
#min f(x) = -x1 - 3x2 sujeito a 
A  = np.array([[-3, 4, 1, 0, 0],
      [1, -1, 0, 1, 0],
      [1, 1, 0, 0, 1]])

N = A[:,:2] #duas primeiras colunas
b = np.array([12, 4, 6])
B = np.identity(3)
c = np.array([-1, -3, 0, 0, 0])

x_basica = np.array([3, 4, 5])
x_n_basica = np.array([1, 2])
c_basica= np.array([0, 0, 0])
c_n_basica = np.array([-1, -3])
contador = 1
while(pare == False):
    print("{} ° iteracao".format(contador))
    contador +=1 

    # Calculo da Solucao Basica
    x_b = np.linalg.solve(B,b)
    res = f(x_b)
    print('Solucao Basica')
    print(x_b)
    print('\n')

    #Vetor Multiplicador Simplex
    lamb = np.linalg.solve(np.transpose(B), c_basica)

    #Calculo dos Custos Relativos
    for i in range(np.shape(c_n_basica)[0]):
        c_rel = c_n_basica - np.dot(np.transpose(lamb), N[:,:i+1])

    #Teste de otimalidade
    aux = 0
    valor_min = c_rel[0]
    indice_min = 0
    for i in range(np.shape(c_rel)[0]):
        #Verificando qual variável entra na base (caso nao estamos na solucao otima)
        if(c_rel[i] <valor_min):
            valor_min = c_rel[i]
            indice_min = i

        if(c_rel[i] >= 0):
            aux += 1

        if(i == np.shape(c_rel)[0]-1 and aux == i+1):
            pare = True

    #A variavel x_n_basica[indice_min] entra na base
    entra = indice_min

    #Calculo da direcao simplex
    y = np.linalg.solve(B, A[:,indice_min])

    aux = 0
    valor_min = 1e4
    indice_min = 0

    #Tamanho do passo
    for i in range(np.shape(y)[0]):
        #Verificando qual variável sai da base
        if(x_b[i]/y[i] <valor_min and y[i]>0):
            valor_min = x_b[i]/y[i]
            indice_min = i

        if(y[i] <= 0):
            aux += 1

        if(i == np.shape(y)[0] and aux == i+1):
            sol_infinita = True
            pare = True 

    #A variavel x_basica[indice_min] sai da base
    sai = indice_min

    #Atualizacao das bases
    B[:,sai] = A[:,x_n_basica[entra]-1]
    N[:,entra] = A[:,x_basica[sai]-1]

    aux = x_basica[sai]
    x_basica[sai] = x_n_basica[entra]
    x_n_basica[entra] = aux

    c_basica = [c[x_basica[i]-1] for i in range(np.shape(x_basica)[0])]
    c_n_basica = [c[x_n_basica[i]-1] for i in range(np.shape(x_n_basica)[0])]

print("Solucao Otima: {}". format(x_b))