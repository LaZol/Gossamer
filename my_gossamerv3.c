#include <stdio.h>
#include <mpir.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define MAX_EXEC 10500

int PRINT = 0;
//FILE *fp1, *fp2;


typedef struct tag_data_execution {
    mpz_t IDS;
    mpz_t K1;
    mpz_t K2;
} tag_data_execution;

typedef struct tag {
    mpz_t ID;
    tag_data_execution tag_data_execs[MAX_EXEC];
} tag;



int BITS=96;
//int BITS=10;
int MIXBITS_ROUNDS=32;
//auxiliar global variable
mpz_t BITS_mpz, module, one, due, seed;
//State for random numbers
gmp_randstate_t st;



void print_mpz(const char* desc, const mpz_t n, int base) {
	printf("%s:",desc);
	mpz_out_str(stdout,base,n);
	printf("\n");

}

void mpz_add_mod(mpz_t res, const mpz_t op1, const mpz_t op2) {
	mpz_add(res, op1, op2);
	mpz_mod(res, res, module);
}


void left_circular_shift(mpz_t res, const mpz_t n, const mp_bitcnt_t k){
int j, index;
//init result to zero
mpz_init(res);
if (k <0){
    printf("Number of shift can't be negative.");
    return;
}
    for(j=BITS-1; j>=0; j--) {
        index = (j-k);
        if (index < 0) index = BITS + index;
        //if n[index] ==1 => res[j] = 1;
        if (mpz_tstbit(n, index)) mpz_setbit(res, j);
        //printf("n[%ld]: %d\n", index, mpz_tstbit(n, index));
        //printf("res[%ld]: %d\n", j, mpz_tstbit(res, j));
    }
}

void mixbits(mpz_t result, const mpz_t x, const mpz_t y) {
	mpz_t z,t1,t2;
	mpz_init_set(z,x);
	mpz_inits(t1,t2, NULL);
	for (int i=0; i<MIXBITS_ROUNDS; i++) {
        //t1 = z>>1
		mpz_tdiv_q(t1,z,due);
        //t2=z+t1
		mpz_add(t2,z,t1);
		//t1= t2+z
		mpz_add(t1,t2,z);
		//t2=t1+y
		mpz_add(t2,t1,y);
		//z=t2mod96
		mpz_mod(z,t2,module);
	}
	mpz_set(result,z);
	mpz_clears(z,t1,t2,NULL);
}

int main(int argc, char* argv[]) {
    struct timeval stop, start;
    int fprint= 0;
    long long int r;
    long long int microsec, sec, tot_sec;
    stop.tv_sec=0;
    start.tv_sec = 0;
    stop.tv_usec =0;
    start.tv_usec =0;
    tot_sec = 0;

    

	//Check optional parameters
	if (argv[1]) {
		if (!strncmp(argv[1],"-p",2)) {
			PRINT = 1;
            }
		}


    int tags_num = 0;
    int exec_num = 0;

    mpz_init(seed);
    //print_mpz("seed", seed, 16);

    /* Intializes random number generator */

    mpz_init_set_str(seed,"40407349974688373120161327220",16);

    //Init random generators

    srand(time(NULL));
    r = rand();
    mpz_set_si(seed, r);

	gmp_randinit_lc_2exp_size(st, 128);
	gmp_randseed (st, seed);



	tags_num = 1;
	exec_num = 1000;

    //printf("Number of tags:");
    //scanf("%d", &tags_num);
    //printf("Number of execution:");
    //scanf("%d", &exec_num);

    if (exec_num>MAX_EXEC){
    	printf("\nNumber of execution can't be more than %d.\n", MAX_EXEC);
    }


    tag tags[tags_num];
    for (int i=0; i<tags_num ;i++){
        mpz_inits(tags[i].ID, NULL);
    }

    for (int i=0; i<tags_num ;i++){
        mpz_urandomb(tags[i].ID, st, BITS);
        print_mpz("ID", tags[i].ID, 10);
    }


    //INITIALIZING VARIABLE
    mpz_t nk1,nk2,pi,n1,n2,n3,n4,tmp1,tmp2,tmp3, num_shift, a, b, c, d;

    mpz_init(BITS_mpz);
	mpz_set_si(BITS_mpz,BITS);
	//print_mpz("BITS_mpz",BITS_mpz, 10);

	//one and two mpz_t
	mpz_init_set_str(one,"1",10);
	mpz_init_set_str(due,"2",10);

	mpz_init(module);
	mpz_mul_2exp(module,one,BITS);

    mpz_init_set_str(pi,"3243F6A8885A308D313198A2",16);

    /*
    if (tags_num ==2){
            printf("2 TAG\n");
            printf("Number of execution: %d", exec_num);
        if ( (fp1 = fopen("neural_network/training_set.csv", "w") ) == NULL )  {
            printf("Error: can't create file.");
            exit(0);
           }
        if ( (fp2 = fopen("neural_network/test_set.csv", "w") ) == NULL ){
            printf("Error: can't create file.");
            exit(0);
        }
    }*/
    

    //start_t=time(NULL);

    for (int i=0; i<tags_num ;i++){

       for(int k=0; k<exec_num;k++ ){

            mpz_inits(nk1,nk2,n1,n2,n3,n4,tmp1,tmp2,tmp3, num_shift, a, b, c, d, NULL);

            gettimeofday(&start, NULL);
            mpz_urandomb(n1,st,BITS);
            mpz_urandomb(n2,st,BITS);

            if (PRINT) printf("\n");
            if (PRINT) printf("Tag %d execution %d :", i, k);
            if (PRINT) printf("\n");
            if (PRINT) print_mpz("n1",n1, 10);
            if (PRINT) print_mpz("n2",n2, 10);

            if(k==0){
                mpz_urandomb(tags[i].tag_data_execs[k].IDS, st, BITS);
                mpz_urandomb(tags[i].tag_data_execs[k].K1, st, BITS);
                mpz_urandomb(tags[i].tag_data_execs[k].K2, st, BITS);

                if (PRINT) printf("\n");
                if (PRINT) print_mpz("k1",tags[i].tag_data_execs[k].K1, 10);
                if (PRINT) print_mpz("k2",tags[i].tag_data_execs[k].K2, 10);
                if (PRINT) print_mpz("IDS",tags[i].tag_data_execs[k].IDS, 10);
                if (PRINT) print_mpz("ID",tags[i].ID, 10);
                if (PRINT) printf("\n");

           }

           ////EXCUTION

            if (PRINT) printf("\n");
            if (PRINT) printf("Computing of A,B,C,D:\n");
            //A=ROT(ROT(IDS+k1+pi+n1, k2)+k1, k1)
            //	//IDS+k1+pi+n1
            mpz_add_mod(tmp1,tags[i].tag_data_execs[k].IDS,tags[i].tag_data_execs[k].K1);
            mpz_add_mod(tmp2,pi,tmp1);
            mpz_add_mod(tmp1,n1,tmp2);
            //	//k2mod96
            mpz_mod(num_shift,tags[i].tag_data_execs[k].K2,BITS_mpz);
            //	//ROT(IDS+k1+pi+n1, k2mod96)
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            //+ k1
            mpz_add_mod(tmp1,tags[i].tag_data_execs[k].K1,tmp2);
            //ROT(..., k1mod96)
            mpz_mod(num_shift,tags[i].tag_data_execs[k].K1,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            mpz_set(a,tmp2);
            if (PRINT) print_mpz("A",a, 10);
              

            //B=ROT(ROT(IDS+k2+pi+n2, k1)+k2, k2)
            mpz_add_mod(tmp1,tags[i].tag_data_execs[k].IDS,tags[i].tag_data_execs[k].K2);
            mpz_add_mod(tmp2,pi,tmp1);
            mpz_add_mod(tmp1,n2,tmp2);
            mpz_mod(num_shift,tags[i].tag_data_execs[k].K1,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            mpz_add_mod(tmp2,tags[i].tag_data_execs[k].K2,tmp1);
            mpz_mod(num_shift,tags[i].tag_data_execs[k].K2,BITS_mpz);
            left_circular_shift(tmp1, tmp2, mpz_get_si(num_shift));
            mpz_set(b,tmp1);
            if (PRINT) print_mpz("B",b, 10);


            //n3
            mixbits(n3, n1, n2);
            if (PRINT) print_mpz("n3",n3, 10);

            //nk1 and nk2
            mpz_add_mod(tmp1,n2,tags[i].tag_data_execs[k].K1);
            mpz_add_mod(tmp2,pi,tmp1);
            mpz_add_mod(tmp1,n3,tmp2);
            mpz_mod(num_shift,n2,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            // k2 XOR n3
            mpz_xor(tmp3, tags[i].tag_data_execs[k].K2, n3);
            //adding XOR to shifted value
            mpz_add_mod(tmp1,tmp2,tmp3);
            mpz_mod(num_shift,n1,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            mpz_xor(tmp1, tmp2, n3);
            mpz_set(nk1,tmp1);
            if (PRINT) print_mpz("nk1",nk1, 10);



            mpz_add_mod(tmp1,n1,tags[i].tag_data_execs[k].K2);
            mpz_add_mod(tmp2,pi,tmp1);
            mpz_add_mod(tmp1,n3,tmp2);
            mpz_mod(num_shift,n1,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            mpz_add_mod(tmp3, tags[i].tag_data_execs[k].K1, n3);
            mpz_add_mod(tmp1,tmp2,tmp3);
            mpz_mod(num_shift,n2,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            mpz_add_mod(tmp1, tmp2, n3);
            mpz_set(nk2,tmp1);
            if (PRINT) print_mpz("nk2",nk2, 10);

            //n4
            mixbits(n4, n3, n2);
            if (PRINT) print_mpz("n4",n4, 10);

            //C
            mpz_add_mod(tmp1,n3,nk1);
            mpz_add_mod(tmp2,pi,tmp1);
            mpz_add_mod(tmp1,n4,tmp2);
            mpz_mod(num_shift,n3,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            // nk2 XOR n4
            mpz_xor(tmp3, nk2, n4);
            //adding XOR to shifted value
            mpz_add_mod(tmp1,tmp2,tmp3);
            mpz_mod(num_shift,n2,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            mpz_xor(tmp1, tmp2, n4);
            mpz_set(c,tmp1);
            if (PRINT) print_mpz("C",c, 10);


            //D
            mpz_add_mod(tmp1,n2,nk2);
            mpz_add_mod(tmp2,tags[i].ID,tmp1);
            mpz_add_mod(tmp1,n4,tmp2);
            mpz_mod(num_shift,n2,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            mpz_add_mod(tmp3, nk1, n4);
            mpz_add_mod(tmp1,tmp2,tmp3);
            mpz_mod(num_shift,n3,BITS_mpz);
            left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
            mpz_add_mod(tmp1, tmp2, n4);
            mpz_set(d,tmp1);
            if (PRINT) print_mpz("D",d, 10);
            
            gettimeofday(&stop, NULL);
            sec = stop.tv_sec - start.tv_sec;
            microsec = stop.tv_usec - start.tv_usec;
            tot_sec =tot_sec+ microsec;
        
            printf("This esecution took %llu sec %llu sec^-6\n", sec, microsec);
        
/*
            if (tags_num ==2){
                if ((i==0)||(i==1)){
                        if (k<9998){
                            fprintf(fp1, "%d,%d,", i, k);
                            mpz_out_str (fp1, 2, tags[i].tag_data_execs[k].IDS);
                            fprintf(fp1, ",");
                            mpz_out_str (fp1, 2, a);
                            fprintf(fp1, ",");
                            mpz_out_str (fp1, 2, b);
                            fprintf(fp1, ",");
                            mpz_out_str (fp1, 2, c);
                            fprintf(fp1, ",");
                            mpz_out_str (fp1, 2, d);
                            fprintf(fp1, "\n");
                            fflush(fp1);
                            }else {
                            fprintf(fp2, "%d, %d,", i, k);
                            mpz_out_str (fp2, 2, tags[i].tag_data_execs[k].IDS);
                            fprintf(fp2, ",");
                            mpz_out_str (fp2, 2, a);
                            fprintf(fp2, ",");
                            mpz_out_str (fp2, 2, b);
                            fprintf(fp2, ",");
                            mpz_out_str (fp2, 2, c);
                            fprintf(fp2, ",");
                            mpz_out_str (fp2, 2, d);
                            fprintf(fp2, "\n");
                            fflush(fp2);

                    }
                }

                    /*if (i==0){
                        //fprint = fprintf(fp1, "[");
                        //printf("fprintf: %d\n", fprint);
                        fprintf(fp1, "%d, %d, ", i, k);
                        mpz_out_str (fp1, 2, tags[i].execs[k].D);
                        fprintf(fp1, "\n");
                        fflush(fp1);
                        }
                    if (i==1){
                        //fprint = fprintf(fp1, "[");
                        //printf("fprintf: %d\n", fprint);
                        fprintf(fp2, "%d, %d, ", i, k);
                        mpz_out_str (fp2, 2, tags[i].execs[k].D);
                        fprintf(fp2, "\n");
                        fflush(fp2);
                    }*/

                //}


                /////////////////////////////UPDATING
                if (PRINT) printf("\n");
                if (PRINT) printf("Updating!\n");
                if (PRINT) print_mpz("ID",tags[i].ID, 10);

                ////LOOSING OLD N1!
                mixbits(n1,n4,n3);

                ///////////////IDSnext
                mpz_add_mod(tmp1,n4,nk1);
                mpz_add_mod(tmp2,tags[i].tag_data_execs[k].IDS,tmp1);
                mpz_add_mod(tmp1,n1,tmp2);
                mpz_mod(num_shift,n4,BITS_mpz);
                left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
                // nk2 XOR n4
                mpz_xor(tmp3, nk2, n4);
                //adding XOR to shifted value
                mpz_add_mod(tmp1,tmp2,tmp3);
                mpz_mod(num_shift,n3,BITS_mpz);
                left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
                mpz_xor(tmp1, tmp2, n1);
                mpz_set(tags[i].tag_data_execs[k+1].IDS,tmp1);
                if (PRINT) print_mpz("newIDS",tags[i].tag_data_execs[k+1].IDS, 10);

                ///////////////////////k1next
                mpz_add_mod(tmp1,n3,nk2);
                mpz_add_mod(tmp2,pi,tmp1);
                mpz_add_mod(tmp1,n1,tmp2);//new n1
                mpz_mod(num_shift,n3,BITS_mpz);
                left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
                mpz_add_mod(tmp3, nk1, n1);//new n1
                mpz_add_mod(tmp1,tmp2,tmp3);
                mpz_mod(num_shift,n4,BITS_mpz);
                left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
                mpz_add_mod(tmp1, tmp2, n1);//new n1
                mpz_set(tags[i].tag_data_execs[k+1].K1,tmp1);
                if (PRINT) print_mpz("new k1",tags[i].tag_data_execs[k+1].K1, 10);

                //k2next
                mpz_add_mod(tmp1,tags[i].tag_data_execs[k+1].IDS,nk2); //nextIDS
                mpz_add_mod(tmp2,pi,tmp1);
                mpz_add_mod(tmp1,tags[i].tag_data_execs[k+1].K1,tmp2);//k1 next
                mpz_mod(num_shift,tags[i].tag_data_execs[k+1].IDS,BITS_mpz);//IDS next
                left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
                mpz_add_mod(tmp3, nk1, tags[i].tag_data_execs[k+1].K1);//k1 next
                mpz_add_mod(tmp1,tmp2,tmp3);
                mpz_mod(num_shift,n1,BITS_mpz);//n1 next
                left_circular_shift(tmp2, tmp1, mpz_get_si(num_shift));
                mpz_add_mod(tmp1, tmp2, tags[i].tag_data_execs[k+1].K1); //k1 next
                mpz_set(tags[i].tag_data_execs[k+1].K2,tmp1);
                if (PRINT) print_mpz("new k2",tags[i].tag_data_execs[k+1].K2, 10);
                if (PRINT) printf("\n");





        }

    }
   
    printf("Total time: %lu microsec\n", tot_sec);
    tot_sec = tot_sec/tags_num;
    tot_sec = tot_sec/exec_num;

    printf("Average time for one execution: %lu microsec\n", tot_sec);
    
    //fclose(fp1);
    //fclose(fp2);
    mpz_clears(nk1,nk2,pi,n1,n2,n3,n4,tmp1,tmp2,tmp3,num_shift,a, b, c, d,NULL);

    exit(0);
}


