/*  Universidade Estadual de Londrina — 2014
    Engenharia Elétrica

    Criptografia de Curvas Elípticas

    Fabrício Leão Rennó Salomon
    Giovani Augusto de Lima Freitas
    Guilherme Brandão da Silva
    Ricardo Fujita
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

#define ORDEM 255

struct _curva_eliptica{
	mpz_t a;
	mpz_t b;
	mpz_t p;
};

struct _ponto{
	mpz_t x;
	mpz_t y;
};

struct _curva_eliptica EC;

void inicializa_ec(){
    printf( "\t\t  +--------------------------------------------+\n"
            "\t\t  | A curva E(p) tem forma y^2 = x^3 + a*x + b |\n"
            "\t\t  +--------------------------------------------+\n\n");

    printf("Informe a: ");
    gmp_scanf("%Zd", EC.a);
    printf("Informe b: ");
    gmp_scanf("%Zd", EC.b);
    printf("Informe p: ");
    gmp_scanf("%Zd", EC.p);
    
}

void dobra_ponto(struct _ponto *P, struct _ponto *R){
	mpz_t curva,temp;
    mpz_init(temp);
	mpz_init(curva);

	if(mpz_cmp_ui(P->y,0)!=0){
		mpz_mul_ui(temp,P->y,2);
		mpz_invert(temp,temp,EC.p);
		mpz_mul(curva,P->x,P->x);
		mpz_mul_ui(curva,curva,3);
		mpz_add(curva,curva,EC.a);
		mpz_mul(curva,curva,temp);
		mpz_mod(curva,curva,EC.p);
		mpz_mul(R->x,curva,curva);
		mpz_sub(R->x,R->x,P->x);
		mpz_sub(R->x,R->x,P->x);
		mpz_mod(R->x,R->x,EC.p);
		mpz_sub(temp,P->x,R->x);
		mpz_mul(R->y,curva,temp);
		mpz_sub(R->y,R->y,P->y);
		mpz_mod(R->y,R->y,EC.p);
	}else{
		mpz_set_ui(R->x,0);
		mpz_set_ui(R->y,0);
	}
}

void soma_ponto(struct _ponto *P, struct _ponto *Q, struct _ponto *R){
    mpz_t temp,curva;
	mpz_mod(P->x,P->x,EC.p);
	mpz_mod(P->y,P->y,EC.p);
	mpz_mod(Q->x,Q->x,EC.p);
	mpz_mod(Q->y,Q->y,EC.p);
	mpz_init(temp);
	mpz_init_set_ui(curva,0);

	if(mpz_cmp_ui(P->x,0)==0 && mpz_cmp_ui(P->y,0)==0){
		mpz_set(R->x,Q->x);
		mpz_set(R->y,Q->y);
		return;
	}

	if(mpz_cmp_ui(Q->x,0)==0 && mpz_cmp_ui(Q->y,0)==0){
		mpz_set(R->x,P->x);
		mpz_set(R->y,P->y);
		return;
	}

	if(mpz_cmp_ui(Q->y,0)!=0){
		mpz_sub(temp,EC.p,Q->y);
		mpz_mod(temp,temp,EC.p);
	}else{
		mpz_set_ui(temp,0);
	}

	if(mpz_cmp(P->y,temp)==0 && mpz_cmp(P->x,Q->x)==0){
		mpz_set_ui(R->x,0);
		mpz_set_ui(R->y,0);
		return;
	}

	if(mpz_cmp(P->x,Q->x)==0 && mpz_cmp(P->y,Q->y)==0){
		dobra_ponto(P,R);
		return;
	}

	else{
		mpz_sub(temp,P->x,Q->x);
		mpz_mod(temp,temp,EC.p);
		mpz_invert(temp,temp,EC.p);
		mpz_sub(curva,P->y,Q->y);
		mpz_mul(curva,curva,temp);
		mpz_mod(curva,curva,EC.p);
		mpz_mul(R->x,curva,curva);
		mpz_sub(R->x,R->x,P->x);
		mpz_sub(R->x,R->x,Q->x);
		mpz_mod(R->x,R->x,EC.p);
		mpz_sub(temp,P->x,R->x);
		mpz_mul(R->y,curva,temp);
		mpz_sub(R->y,R->y,P->y);
		mpz_mod(R->y,R->y,EC.p);
 		return;
	}
}

void mult_escalar(struct _ponto *P, struct _ponto *R, mpz_t m){
	struct _ponto Q,T;
	mpz_init(Q.x);
	mpz_init(Q.y);
	mpz_init(T.x);
	mpz_init(T.y);

	long len,loop;

	len = mpz_sizeinbase(m,2);
	mpz_set_ui(R->x,0);
	mpz_set_ui(R->y,0);

	if(mpz_cmp_ui(m,0)==0){
		return;
	}

    mpz_set(Q.x,P->x);
    mpz_set(Q.y,P->y);

	if(mpz_tstbit(m,0)==1){
		mpz_set(R->x,P->x);
		mpz_set(R->y,P->y);
	}

	for(loop = 1; loop < len; loop++){
		mpz_set_ui(T.x,0);
		mpz_set_ui(T.y,0);
		dobra_ponto(&Q,&T);
		mpz_set(Q.x,T.x);
		mpz_set(Q.y,T.y);
		mpz_set(T.x,R->x);
		mpz_set(T.y,R->y);

		if(mpz_tstbit(m,loop)){
			soma_ponto(&T,&Q,R);
		}
	}
}

void gera_numero(mpz_t *ka, mpz_t *kb, mpz_t *k){
	gmp_randstate_t r_state;
    gmp_randinit_default (r_state);

    gmp_randseed_ui(r_state, clock());

	mpz_urandomb(*ka,r_state, ORDEM);
	mpz_urandomb(*kb,r_state, ORDEM);
	mpz_urandomb(*k,r_state, ORDEM);
}

void inverte_coordenada(struct _ponto *P){
	mpz_mul_si(P->y, P->y, -1);
}

int is_on_curve(struct _ponto *P){
    mpz_t i, j, aux;

    mpz_init(i);
    mpz_init(j);
	mpz_init(aux);

	mpz_mod(j,P->y,EC.p);
	mpz_pow_ui(j,j,2);
	mpz_mod(j,j,EC.p);

    mpz_pow_ui(i,P->x,3);
    mpz_mod(i,i,EC.p);
    mpz_mul(aux,P->x,EC.a);
    mpz_mod(aux,aux,EC.p);
    mpz_add(aux,aux,EC.b);
    mpz_mod(aux,aux,EC.p);
    mpz_add(i,i,aux);
    mpz_mod(i,i,EC.p);

    if(!mpz_cmp(i,j)){
        return 1;
    }
    else{
        return 0;
    }
}

void left_to_right(struct _ponto *P, struct _ponto *R, mpz_t k){
    if(!mpz_cmp_d(k,1)){
        mpz_set(R->x,P->x);
        mpz_set(R->y,P->y);
		return;
    }else if(!mpz_cmp_d(k,0)){
        mpz_set_ui(R->x, 0);
        mpz_set_ui(R->y, 0);
        return;
    }

	int i = 0;
	struct _ponto r_0, r_1, aux, aux1;

	mpz_init_set(aux.x,P->x);
	mpz_init_set(aux.y,P->y);

	mpz_init(aux1.x);
	mpz_init(aux1.y);

	mpz_init_set(r_0.x,P->x);
	mpz_init_set(r_0.y,P->y);
	mpz_init_set(r_1.x,P->x);
	mpz_init_set(r_1.y,P->y);

	for (i=(mpz_sizeinbase(k,2))-2;i>=0;i--){
		dobra_ponto(&r_0, &aux);

		if(mpz_tstbit(k,i)){
			soma_ponto(&aux,&r_1,&aux1);

			mpz_set(aux.x, aux1.x);
			mpz_set(aux.y, aux1.y);
		}
		mpz_set(r_0.x, aux.x);
		mpz_set(r_0.y, aux.y);
	}

	if(mpz_tstbit(k,0)){
        mpz_set(R->x,aux1.x);
        mpz_set(R->y,aux1.y);
		return;
	}
    mpz_set(R->x,aux.x);
    mpz_set(R->y,aux.y);
	return;
}

void diffie_hellman(struct _ponto *G, struct _ponto *P, struct _ponto *Q){
	mpz_t ka, kb, k;
	struct _ponto pa, pb, c1, c2, kb_c1, Pm, aux;

	mpz_init(ka);
	mpz_init(kb);
	mpz_init(k);

	mpz_init(pa.x);
	mpz_init(pb.x);
	mpz_init(c1.x);
	mpz_init(c2.x);
	mpz_init(pa.y);
	mpz_init(pb.y);
	mpz_init(c1.y);
	mpz_init(c2.y);

	mpz_init(kb_c1.x);
	mpz_init(kb_c1.y);

	mpz_init(Pm.x);
	mpz_init(Pm.y);

	mpz_init(aux.x);
	mpz_init(aux.y);

	gera_numero(&ka,&kb,&k);
	gmp_printf("ka = %Zd\n\nkb = %Zd\n\nk = %Zd\n\n", ka, kb, k);

    left_to_right(G, &pa, ka);
	left_to_right(G, &pb, kb);
    gmp_printf("Pa = [%Zd, %Zd]\n\nPb = [%Zd, %Zd]\n\n", pa.x, pa.y, pb.x, pb.y);
    left_to_right(G, &c1, k);
	left_to_right(&pb, &aux, k);
	soma_ponto(P, &aux, &c2);
    left_to_right(&c1, &kb_c1, kb);
	inverte_coordenada(&kb_c1);
	soma_ponto(&c2,&kb_c1, &Pm);

	printf(    "+---------------------------------------------------------+\n");
	gmp_printf("| C1 = [%Zd, %Zd]\n", c1.x, c1.y);
	gmp_printf("| C2 = [%Zd, %Zd]\n", c2.x, c2.y);
	printf(    "+---------------------------------------------------------+\n\n");

	gmp_printf("Pm = [%Zd, %Zd]\n\nP = [%Zd, %Zd]\n", Pm.x, Pm.y, P->x, P->y);

    mpz_set(Q->x, P->x);
    mpz_set(Q->y, P->y);
}

long long int decryption(struct _ponto *G, struct _ponto *P, long long int order){
    struct _ponto R;
    long long int i;

    mpz_t j;
    mpz_init(j);
    mpz_init(R.x);
    mpz_init(R.y);

    for(i=0; i<=order; i++){
        mpz_set_ui(j, i);
        left_to_right(G, &R, j);
        if(!mpz_cmp(R.x,P->x) && !mpz_cmp(R.y,P->y)){
            return i;
        }
    }
    return -1;
}

void help(){
    printf("\n\n");

    FILE *ptr_file;
    char buf[1000];

    ptr_file = fopen("help.txt","r");
    if(!ptr_file){
        return;
    }
    while(fgets(buf,1000, ptr_file)!=NULL){
        printf("%s",buf);
    }
    fclose(ptr_file);
    printf("\n\n");
}

//main----------------------------------------------------------------------------------------------
void main(){
    printf( "\t\t     Universidade Estadual de Londrina - 2014\n"
            "\t\t\t\tEngenharia Eletrica\n\n"
            "\t\t\t Criptografia de Curvas Elipticas\n\n"
            "\t\t\t   Fabricio Leao Renno Salomon\n"
            "\t\t\t Giovani Augusto de Lima Freitas\n"
            "\t\t\t   Guilherme Brandao da Silva\n"
            "\t\t\t          Ricardo Fujita\n"
            "\t\t     ________________________________________\n\n"
            "\t\t\t   /$$$$$$$$  /$$$$$$   /$$$$$$\n"
            "\t\t\t  | $$_____/ /$$__  $$ /$$__  $$\n"
            "\t\t\t  | $$      | $$  |__/| $$  |__/\n"
            "\t\t\t  | $$$$$   | $$      | $$      \n"
            "\t\t\t  | $$__/   | $$      | $$      \n"
            "\t\t\t  | $$      | $$    $$| $$    $$\n"
            "\t\t\t  | $$$$$$$$|  $$$$$$/|  $$$$$$/\n"
            "\t\t\t  |________/ |______/  |______/ \n\n");

    int op;
    mpz_init(EC.a);
    mpz_init(EC.b);
    mpz_init(EC.p);

    inicializa_ec();

    struct _ponto P,R;

    mpz_init(P.x);
    mpz_init(P.y);
    mpz_init_set_ui(R.x,0);
    mpz_init_set_ui(R.y,0);

    while(1){
        printf("\n\t\t\t+---------------Opcao--------------+\n"
                "\t\t\t|1 ->        Soma de pontos        |\n"
                "\t\t\t|2 ->   Multiplicacao por escalar  |\n"
                "\t\t\t|3 ->      Validacao de ponto      |\n"
                "\t\t\t|4 -> Troca de chave Diffie-Hellman|\n"
                "\t\t\t|5 ->    Decriptar codigo v.Beta   |\n"
                "\t\t\t|6 ->           Ajuda              |\n"
                "\t\t\t|0 ->            Sair              |\n"
                "\t\t\t+----------------------------------+\n\n");

        printf("Opcao: ");
        scanf("%d",&op);
       
        if(op==1){
            printf("\t\t\t +--------------------------------+\n"
                    "\t\t\t | Entrou opcao 1: Soma de pontos |\n"
                    "\t\t\t +--------------------------------+\n\n"
                    "+----------------+\n"
                    "| Informe P(x,y) |\n"
                    "+----------------+\n");

            do{
                printf("Informe P.x: ");
                gmp_scanf("%Zd",&P.x);
                printf("Informe P.y: ");
                gmp_scanf("%Zd",&P.y);
            }while(!is_on_curve(&P));

            printf( "+----------------+\n"
                    "| Informe Q(x,y) |\n"
                    "+----------------+\n");

            struct _ponto Q;
            mpz_init(Q.x);
            mpz_init(Q.y);

            do{
                printf("Informe Q.x: ");
                gmp_scanf("%Zd",&Q.x);
                printf("Informe Q.y: ");
                gmp_scanf("%Zd",&Q.y);
            }while(!is_on_curve(&Q));

            soma_ponto(&P,&Q,&R);

            printf(	   "+---------------------------------------------------------+\n");
            gmp_printf("| O ponto P e: [%Zd, %Zd]\n", R.x, R.y);
            printf(	   "+---------------------------------------------------------+\n");

        }else if(op==2){
            printf( "\t\t    +-------------------------------------------+\n"
                    "\t\t    | Entrou opcao 2: Multiplicacao por escalar |\n"
                    "\t\t    +-------------------------------------------+\n\n"
                    "+----------------+\n"
                    "| Informe P(x,y) |\n"
                    "+----------------+\n");

            do{
                printf("Informe P.x: ");
                gmp_scanf("%Zd",&P.x);
                printf("Informe P.y: ");
                gmp_scanf("%Zd",&P.y);
            }while(!is_on_curve(&P));

            struct _ponto Q;
            mpz_init(Q.x);
            mpz_init(Q.y);

            printf( "+-------------------+\n"
                    "| Informe o escalar |\n"
                    "+-------------------+\n");

            mpz_t m;
            mpz_init(m);
            gmp_scanf("%Zd",&m);

            left_to_right(&P, &Q, m);

            printf(	   "+---------------------------------------------------------+\n");
            gmp_printf("| O ponto P e: [%Zd, %Zd]\n",Q.x,Q.y);
            printf(	   "+---------------------------------------------------------+\n");


        }else if(op==3){
            printf( "\t\t       +------------------------------------+\n"
                    "\t\t       | Entrou opcao 3: Validacao de ponto |\n"
                    "\t\t       +------------------------------------+\n\n"
                    "+----------------+\n"
                    "| Informe P(x,y) |\n"
                    "+----------------+\n");

            printf("Informe P.x: ");
            gmp_scanf("%Zd",&P.x);
            printf("Informe P.y: ");
            gmp_scanf("%Zd",&P.y);

            printf(	   "+----------------------------+\n");
            printf(    "|        retornou %d          |\n",is_on_curve(&P));
            printf(	   "+----------------------------+\n");

        }else if(op==4){
            printf( "\t\t\t +--------------------------------+\n"
                    "\t\t\t | Entrou opcao 4: Diffie-Hellman |\n"
                    "\t\t\t +--------------------------------+\n\n"
                    "+-------------------------+\n"
                    "| Informe o ponto gerador |\n"
                    "+-------------------------+\n");

            struct _ponto G;
            mpz_init(G.x);
            mpz_init(G.y);
            do{
                printf("Informe G.x: ");
                gmp_scanf("%Zd",&G.x);
                printf("Informe G.y: ");
                gmp_scanf("%Zd",&G.y);
            }while(!is_on_curve(&G));

            printf( "+------------------------+\n"
                    "| Informe P (informacao) |\n"
                    "+------------------------+\n");

            do{
                printf("Informe P.x: ");
                gmp_scanf("%Zd",&P.x);
                printf("Informe P.y: ");
                gmp_scanf("%Zd",&P.y);
            }while(!is_on_curve(&P));
            printf("\n");

            struct _ponto Q;
            mpz_init(Q.x);
            mpz_init(Q.y);

            diffie_hellman(&G, &P, &Q);

            printf(	 "\n+---------------------------------------------------------+\n");
            gmp_printf("| Mensagem decriptada: [%Zd, %Zd]\n",Q.x,Q.y);
            printf(	   "+---------------------------------------------------------+\n");

         

        }if(op==5){
            printf( "\t\t\t+----------------------------------+\n"
                    "\t\t\t| Entrou opcao 5: Decriptar codigo |\n"
                    "\t\t\t+----------------------------------+\n\n");

            struct _ponto G;
            mpz_init(G.x);
            mpz_init(G.y);

            long long int i;
            long long int order;

            printf( "\t\t\t      +-----------------------+\n"
                    "\t\t\t      | ATENCAO: AREA HACKER! |\n"
                    "\t\t\t      +-----------------------+\n\n"
                    "+-------------------------------------------+\n"
                    "| Algoritmo tenta descobrir k, sendo P = kG |\n"
                    "+-------------------------------------------+\n"
                    "+----------------+\n"
                    "| Informe G(x,y) |\n"
                    "+----------------+\n");

            do{
                printf("Informe G.x: ");
                gmp_scanf("%Zd",&G.x);
                printf("Informe G.y: ");
                gmp_scanf("%Zd",&G.y);
            }while(!is_on_curve(&G));

            printf( "+----------------+\n"
                    "| Informe P(x,y) |\n"
                    "+----------------+\n");

            do{
                printf("Informe P.x: ");
                gmp_scanf("%Zd",&P.x);
                printf("Informe P.y: ");
                gmp_scanf("%Zd",&P.y);
            }while(!is_on_curve(&P));

            printf("Informe ordem: ");
            scanf("%ld",&order);
            i = decryption(&G, &P, order);

            printf("+---------------------------------------------------------+\n");
            printf("|  k = %lld\n", i);
            printf("+---------------------------------------------------------+\n");


        }if(op==6){
            help();

        }if(op==0){
            exit(0);
        }
    }
}

