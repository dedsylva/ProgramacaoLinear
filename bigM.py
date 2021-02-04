import numpy as np

class BigM():
    def __init__(self, A,b,c):

        self.A = A
        self.b = b
        self.c = c

    def f(self,x, c_basica):
        return np.dot(c_basica.T, x)

    def ImprimeSolucao(self,x_b, x_basica, x_n_basica, res):
        item = [i+1 for i in range(len(self.c))]
        j = 0
        sol = []

        while(j < len(self.c)):
            for i in range(len(x_basica)):
                if j == len(self.c):
                    break
                if x_basica[i] == item[j]:
                    j +=1
                    sol.append(x_b[i])


            for i in range(len(x_n_basica)):
                if j == len(self.c):
                    break                

                if x_n_basica[i] == item[j]:                   
                    j +=1
                    sol.append(0)
              
        print('Solução Otima: {}\n'.format(sol))
        print('Valor Ótimo: {}\n'.format(res))

     
    def resolve(self): 
        M = 1e8
        #trocamos as linhas de A se algum b for negativo
        for i in range(len(self.b)):
            if(self.b[i] <0):
                self.b[i] = -self.b[i]
                self.A[i:i+1,:] = -self.A[i:i+1,:]

        #Adicionar variaveis auxiliares y
        A2 = np.concatenate((self.A, np.identity(self.b.shape[0])), axis=1)
        c2 = np.array([self.c[i] if i< self.A.shape[1] else M for i in range(A2.shape[1])])

        pare = False
        N = self.A #matriz original sem as variaveis y auxiliares
        B = np.identity(self.b.shape[0])

        x_basica = np.array([i+N.shape[1] for i in range(1,B.shape[0]+1)])
        x_n_basica = np.array([i for i in range(1,N.shape[1]+1)])

        c_basica = np.array([c2[x_basica[i]-1] for i in range(x_basica.shape[0])])
        c_n_basica = np.array([c2[x_n_basica[i]-1] for i in range(x_n_basica.shape[0])])


        sol_infinita = False
        while(pare == False):
            # Calculo da Solucao Basica
            x_b = np.linalg.solve(B,self.b)               
            res = self.f(x_b, c_basica)

            #Vetor Multiplicador Simplex
            lamb = np.linalg.solve(B.T, c_basica)

            #Calculo dos Custos Relativos
            for i in range(c_n_basica.shape[0]):
                c_rel = c_n_basica[i] - np.dot(lamb.T, N[:,i:i+1])
            #Teste de otimalidade
                if(i == 0):
                    aux = 0
                    valor_min = c_rel[0]
                    indice_min = 0
                #Verificando qual variável entra na base (caso nao estamos na solucao otima)
                if(c_rel[0] <valor_min):
                    valor_min = c_rel[0]
                    indice_min = i

                if(c_rel[0] >= 0):
                    aux += 1

                if(i == c_n_basica.shape[0]-1 and aux == i+1):
                    pare = True

            #A variavel x_n_basica[indice_min] entra na base
            entra = indice_min


            #Calculo da direcao simplex
            y = np.linalg.solve(B, A2[:,indice_min])
            aux = 0
            valor_min = 1e4
            indice_min = 0

            #Tamanho do passo
            for i in range(y.shape[0]):
                #Verificando qual variável sai da base
                if y[i]>0:
                    if x_b[i]/y[i] <valor_min:
                        valor_min = x_b[i]/y[i]
                        indice_min = i

                else:
                    aux += 1


                if(i == y.shape[0]-1 and aux == i+1):
                    sol_infinita = True
                    pare = True 


            if(pare):
                break

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
        
        if(sol_infinita):
            print('O problema não tem solucao finita')

        else:
            self.ImprimeSolucao(x_b, x_basica, x_n_basica, res)

