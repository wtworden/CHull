

#include <gmp.h>
#include <stdio.h>
#include <assert.h>

int main(){

  char inputStr[1024];
  /*
     mpz_t is the type defined for GMP integers.
     It is a pointer to the internals of the GMP integer data structure
   */
  mpz_t a;
  mpz_t b;
  mpz_t d;

  int flag_a;
  int flag_b;

  printf ("Enter the first number: ");
  scanf("%1023s" , inputStr); /* NOTE: never every write a call scanf ("%s", inputStr);
                                  You are leaving a security hole in your code. */
  /* 1. Initialize the number */

  mpz_init(d);
  mpz_set_ui(d,0);

  /* 1. Initialize the number */
  mpz_init(a);
  mpz_set_ui(a,0);

  /* 2. Parse the input string as a base 10 number */
  flag_a = mpz_set_str(a,inputStr, 10);
  assert (flag_a == 0); /* If flag is not 0 then the operation failed */

  printf ("Enter the second number: ");
  scanf("%1023s" , inputStr);

/* 1. Initialize the number */
  mpz_init(b);
  mpz_set_ui(b,0);

  /* 2. Parse the input string as a base 10 number */
  flag_b = mpz_set_str(b,inputStr, 10);
  assert (flag_b == 0); /* If flag is not 0 then the operation failed */


  /* Print n */
  printf ("a = ");
  mpz_out_str(stdout,10,a);
  printf ("\n");

    /* Print n */
  printf ("b = ");
  mpz_out_str(stdout,10,b);
  printf ("\n");

  /* 3. Add one to the number */
  mpz_gcd_ui(d, a, b);

  /* 4. Print the result */

  printf (" gcd(a,b) = ");
  mpz_out_str(stdout,10,d);
  printf ("\n");





  /* 6. Clean up the mpz_t handles or else we will leak memory */
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(d);

}

