#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N (1)
#define NEQ (10) //number of equations

#define a (1.5) 
#define g (0.025) 
#define m (0.1) 
#define d (0.5) 
#define phi (100.0) 


#define AA 471
#define B 1586
#define CC 6988
#define DD 9689
#define M 16383
#define RIMAX 2147483648.0        /* = 2^31 */
#define RandomInteger (ra[nd & M] = ra[(nd-AA) & M] ^ ra[(nd-B) & M] ^ ra[(nd-CC) & M] ^ ra[(nd-DD) & M])

static long ra[M+1] = {}, nd = 0;

void seed(long seed)
{
    int  i;

    if(seed <= 0) { fprintf(stderr, "SEED error!"); exit(1); }
    ra[0] = (long) fmod(16807.0*(double)seed, 2147483647.0);
    for(i = 1; i <= M; i++)
        ra[i] = (long)fmod( 16807.0 * (double) ra[i-1], 2147483647.0);
}

long randl(long num)    // random long integer between 0 and n-1; n must be larger than 0
{
    ++nd;
    long ri = RandomInteger;
    return(ri % num);
}

double randd(void)		// random double between 0 and 1
{
    ++nd;
    return((double) RandomInteger / RIMAX);
}

typedef struct _IND
{
  int X, Y, P, R; // molecule types within the compartments
}COMP;

COMP comp[N];

// filling up the compartments
void init()
{
 int i;

 for(i = 0; i < N; i++)
 {
   comp[i].X = randl(10)+1;
   comp[i].Y = randl(10)+1;
   comp[i].P = 0;
   comp[i].R = randl(10)+1;
 }
  
}

int kk;
double totprop;
double mu[NEQ];
double b;

void Gillespie_step(int i)
{
  int l;
  double r, prop = 0.0;
  
  totprop = 0.0;
  mu[0] = a * comp[i].R * comp[i].X; // R+X --a--> 2X (R--, X++)
  mu[1] = g * comp[i].X * comp[i].Y;  // X+Y --g--> P (X--, Y--, P++)
  mu[2] = m * comp[i].X;  // X --m--> Y (X--, Y++)
  mu[3] = b * comp[i].R * comp[i].Y;  // R+Y --b--> 2Y (R--, Y++)
  mu[4] = m * comp[i].Y;  // Y --m--> X (Y--, X++)
  mu[5] = d * comp[i].R;  // R --d--> 0 (R--)
  mu[6] = phi;  // --phi--> R (R++)
  mu[7] = d * comp[i].X;  // X --d--> 0 (X--)
  mu[8] = d * comp[i].Y;  // Y --d--> 0 (Y--)
  mu[9] = d * comp[i].P;  // P --d--> 0 (P--)
  
  for(l = 0; l < NEQ; l++)
    totprop += mu[l];
  
  r = totprop * randd();

  for(kk = 0; kk < NEQ; kk++)
  {
    prop += mu[kk];
    if(prop > r)
      break;
  }

  if(kk == 0)
  {
    comp[i].R--;
    comp[i].X++;
  }
  else if(kk == 1)
  {
    comp[i].X--;
    comp[i].Y--;
    comp[i].P++;
  }
  else if(kk == 2)
  {
    comp[i].X--;
    comp[i].Y++;
  }
  else if(kk == 3)
  {
    comp[i].R--;
    comp[i].Y++;
  }
  else if(kk == 4)
  {
    comp[i].X++;    
    comp[i].Y--;
  }
  else if(kk == 5)
  {
    comp[i].R--;    
  }
  else if(kk == 6)
  {
    comp[i].R++;
  }
  else if(kk == 7)
  {
    comp[i].X--;
  }
  else if(kk == 8)
  {
    comp[i].Y--;
  }
  else if(kk == 9)
  {
    comp[i].P--;
  }  
  
  else
  {
    printf("ERROR - 1\n");
    printf("totprop=%lf\tprop=%lf\tr=%lf\tkk=%d\n", totprop, prop, r, kk);
    exit(1);    
  }
  
}


int main(void)
{
  int sim, numb;
  double avex, avey, avep, aver, num;

  seed(76231);

  init();
  
  for(sim = 1; sim < 10000; sim++)
    Gillespie_step(0);

  b = 0.0;

  for(numb = 0; numb < 120; numb++)
  {
    avex = 0.0;
    avey = 0.0;
    avep = 0.0;
    aver = 0.0;
    num = 0.0;

    if(numb < 60)
      b += 0.05;
    else
      b -= 0.05;

    for(sim = 1; sim < 1000000; sim++)
    {
      Gillespie_step(0);
      avex += 1.0 * comp[0].X;
      avey += 1.0 * comp[0].Y;
      avep += 1.0 * comp[0].P;
      aver += 1.0 * comp[0].R;
      num += 1.0;
    }
    printf("%lf\t%lf\t%lf\t%lf\t%lf\n", b, avex/num, avey/num, avep/num, aver/num);
    fflush(stdout);
  }

  return 0;
}
