/* -------------------------------------------------------------- */
/* C code for Matsumoto & Nishimura - Mersene Twister             */
/*                                                                */
/* This version is compact with a limited API and 3 new functions */
/* 2 functions to save and restore statuses                       */
/* 1 function mtRand optimized for speed (MT is already one if    */
/*   the fastest high quality generator, if not the fastest with  */
/*   MLFG from Michael Mascagni et al                             */
/*                                                                */
/* Updates by D. Hill - March 2016                                */
/* -------------------------------------------------------------- */

extern "C" {

#include <stdio.h>
#include <stdlib.h>

#include <stdio.h>

/* -------------- Prototypes --------------------------------- */
void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
double mtRand(void);

}
