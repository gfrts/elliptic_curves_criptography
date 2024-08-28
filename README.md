# elliptic_curves_criptography
	Demonstração de Criptografia de Curvas Elipticas (ECC)

## Descrição

	Realiza de maneira simploria o funcionamento de uma
	criptografia de chave assimetrica, analoga ao ElGamal,
	envolvendo curvas elipticas. Possuindo as funcoes basicas
	acerca do topico: soma de pontos, multiplicacao por escalar,
	verificacao se o ponto esta na curva, troca de chaves, e
	ainda uma versao beta de um algoritmo que tenta quebrar a
	criptografia


## Opções

#### 1 - Soma de pontos

	Esta funcao realiza a soma de dois pontos P e Q pertencentes
	a uma mesma curva eliptica E(p) informada pelo usuario. E
	retorna o resultado em um terceiro ponto R


#### 2 - Multiplicacao por escalar

	Utiliza o algoritmo de exponenciacao left-to-right para realizar
	a multiplicacao de um ponto P de uma curva eliptica E(p)
	por um escalar m. Retornando o valor em um segundo ponto R.
	O algoritmo e mais eficiente que os algoritmos convencionais
	de exponenciacao, e e o algoritmo adotado como padrao para
	trocas de chaves de Diffie-Hellman


#### 3 - Validacao de ponto

	A funcao verifica se um dado ponto P pertence a uma curva
	eliptica E(p), retornando 1 se o ponto esta na curva, e 0
	se o ponto nao esta na curva


#### 4 - Troca de chave Diffie-Hellman

	Inserindo-se um ponto gerador G de uma curva eliptica E(p),
	e um ponto que representa uma mensagem P pertencente a
	E(p) e gerado por G, a funcao realiza uma demonstracao de
	troca de chaves de Diffie-Hellman, analoga ao ElGamal,
	e recebe a informacao decriptada em um terceiro ponto Q
	

#### 5 - Decriptar codigo v.Beta

	Dado dois pontos pertencentes a uma mesma curva eliptica
	E(p), G e P, sendo que P = kG para algum escalar k, inserindo-se
	a ordem do ponto G na curva, a funcao tenta encontrar o escalar k,
	tentando assim, de uma maneira ingenua, quebrar a criptografia
	de curvas elipticas
