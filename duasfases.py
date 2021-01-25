import numpy as np

class DuasFases():
    def __init__(self, A,b,c):

        self.A = A
        self.b = b
        self.c = c

        #Fase 1

        #trocamos as linhas de A se algum b for negativo
        for i in range(len(self.b)):
            if(self.b[i] <0):
                self.b[i] = -self.b[i]
                self.A[i:i+1,:] = -self.A[i:i+1,:]
        
        #Adicionar variaveis auxiliares y
        A2 = np.concatenate((A, np.identity(b.shape[0])), axis=1)
        #c2 = np.ones(A2.shape[1])
        c2 = np.array([0 if i< A.shape[1] else 1 for i in range(A2.shape[1])])

        #funcao a max agr eh y1+y2+...+yn, n=b.shape[0]

        pare = False
        N = self.A #matriz original sem as variaveis y auxiliares
        B = np.identity(b.shape[0])

        x_basica = np.array([i+N.shape[1] for i in range(1,B.shape[0]+1)])
        x_n_basica = np.array([i for i in range(1,N.shape[1]+1)])

        c_basica = np.array([c2[x_basica[i]-1] for i in range(x_basica.shape[0])])
        c_n_basica = np.array([c2[x_n_basica[i]-1] for i in range(x_n_basica.shape[0])])


        contador = 1
        while(pare == False):
            print("{} ° iteracao".format(contador))
            contador +=1 

            # Calculo da Solucao Basica
            x_b = np.linalg.solve(B,self.b)
            res = [c_basica[i]*x_b[i] for i in range(len(c_basica))]
            res = sum(res)
            print('Solucao Basica')
            print(x_b)
            print('Valor da funcao')
            print(res)
            print('\n')

            #Vetor Multiplicador Simplex
            lamb = np.linalg.solve(B.T, c_basica)

            #Calculo dos Custos Relativos
            for i in range(c_n_basica.shape[0]):
                c_rel = c_n_basica[i] - np.dot(lamb.T, N[:,:i+1])

            #Teste de otimalidade
            aux = 0
            valor_min = c_rel[0]
            indice_min = 0
            for i in range(c_rel.shape[0]):
                #Verificando qual variável entra na base (caso nao estamos na solucao otima)
                if(c_rel[i] <valor_min):
                    valor_min = c_rel[i]
                    indice_min = i

                if(c_rel[i] >= 0):
                    aux += 1

                if(i == c_rel.shape[0]-1 and aux == i+1):
                    pare = True

            #A variavel x_n_basica[indice_min] entra na base
            entra = indice_min

            #Calculo da direcao simplex
            y = np.linalg.solve(B, self.A[:,indice_min])

            aux = 0
            valor_min = 1e4
            indice_min = 0

            #Tamanho do passo
            for i in range(y.shape[0]):
                #Verificando qual variável sai da base
                if(x_b[i]/y[i] <valor_min and y[i]>0):
                    valor_min = x_b[i]/y[i]
                    indice_min = i

                if(y[i] <= 0):
                    aux += 1

                if(i == y.shape[0] and aux == i+1):
                    sol_infinita = True
                    pare = True 

            #A variavel x_basica[indice_min] sai da base
            sai = indice_min

            #Atualizacao das bases
            B[:,sai] = A2[:,x_n_basica[entra]-1]
            N[:,entra] = A2[:,x_basica[sai]-1]

            aux = x_basica[sai]
            x_basica[sai] = x_n_basica[entra]
            x_n_basica[entra] = aux

            c_basica = np.array([c2[x_basica[i]-1] for i in range(x_basica.shape[0])])
            c_n_basica = np.array([c2[x_n_basica[i]-1] for i in range(x_n_basica.shape[0])])


        if(res != 0 ):
            print('Problema Infactível')
        print("Solucao Otima: {}". format(x_b))




A = np.array(
    [[1, 1,-1,0],
    [-1, 1,0,-1]])

b = np.array([2,1])
c = np.array([-1,2,0,0])

df = DuasFases(A,b,c)
