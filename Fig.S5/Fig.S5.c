#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N (1000) // number of compartments
#define D (2) // number of daughter compartments after division
#define NEQ (10) // number of reactions
#define OWR (1) // besides one daughter replacing the parent, OWR daughter compartments overwrite OWR randomly selected compartments in the population (OWR <= D - 1)
#define DIVSIZE (100) // total molecule count threshold for compartment division

#define a1 (2.0) 
#define b1 (2.4)
#define aa2 (1.5) 
#define bb2 (1.2)
#define gx (0.2) 
#define gy (0.1)
#define hx (0.8)
#define hy (0.6)
#define mx (0.1) 
#define my (0.1)
# define R0 (2.0)
#define c (4)

#define AA 471
#define B 1586
#define CC 6988
#define DD 9689
#define M 16383
#define RIMAX 2147483648.0        /* = 2^31 */
#define RandomInteger (ra[nd & M] = ra[(nd-AA) & M] ^ ra[(nd-B) & M] ^ ra[(nd-CC) & M] ^ ra[(nd-DD) & M])

typedef struct _IND // structure of the compartments
{
  int X1, X2, Y1, Y2, V, W, VP, WP; // molecule types within the compartments
  double a2;
  double b2;
}COMP;
COMP comp[N];

typedef struct _IND2 // structure of the daughter compartments
{
  int X1, X2, Y1, Y2, V, W, VP, WP;
}DAU;
DAU daughter[D];


static long ra[M + 1] = {}, nd = 0;
double time = 0.0, theta;
int numdiv = 0, sim;

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

// the environment affects the autocatalyitic rates
double env(int k)
{
  // CASE 1: environment E_X
  if(theta==1.0){
    if(comp[k].X1 + comp[k].X2 + comp[k].Y1 +comp[k].Y2==0){
      double multi = 1;
      return(multi);
    }
    double multi = 1+c*(((double)comp[k].X1 + (double)comp[k].X2)/((double)comp[k].X1+(double)comp[k].X2 + (double)comp[k].Y1 + (double)comp[k].Y2));
    return(multi);
  }
  // CASE 2: environment E_Y
  if(theta==-1.0){
    if(comp[k].X1 + comp[k].X2 + comp[k].Y1 +comp[k].Y2==0){
      double multi = 1;
      return(multi);
    }
    double multi = 1+c*(((double)comp[k].Y1 + (double)comp[k].Y2)/((double)comp[k].X1+(double)comp[k].X2 + (double)comp[k].Y1 + (double)comp[k].Y2));
    return(multi);
  }
}

void shuffle(int *array, size_t n)  
{
  if(n > 1)
  {
    size_t i;
    for(i = 0; i < n - 1; i++)
    {
      size_t j = i + randl(n - i);
      if((j > n))
          exit(4);
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}

// filling the initial compartments with a random number of X and Y molecules, and a constant number of R molecules
void init() 
{
 int i;

 for(i = 0; i < N; i++)
 {
   comp[i].X1 = randl(10)+1;
   comp[i].Y1 = randl(10)+1;
   comp[i].X2 = 0;
   comp[i].Y2 = 0;
   comp[i].R = R0;
   comp[i].V = 0;
   comp[i].W = 0;
   comp[i].VP = 0;
   comp[i].WP = 0;
   comp[i].a2 = aa2 * env(i);
   comp[i].b2 = bb2 * env(i);
 }
}

int Gillespie_step()
{
  int i, l;
  double totprop = 0.0;
  double mu[N][NEQ];
  double r, prop = 0.0;
  
   for(i = 0; i < N; i++)
  {
    mu[i][0] = a1 * comp[i].R * comp[i].X1; // R+X1 --a1--> X2 (X1--,X2++)
    mu[i][1] = hx * comp[i].X1;  // X1 --hx--> V (X1--, V++)
    mu[i][2] = gx * comp[i].X1 * comp[i].W;  // X1+W --gx--> WP (X1--, W--, WP++)
    mu[i][3] = mx * comp[i].X1;  // X1 --mx--> Y1 (X1--, Y1++)
    mu[i][4] = my * comp[i].Y1;  // Y1 --my--> X1 (Y1--, X1++)
    mu[i][5] = gy * comp[i].Y1 * comp[i].V;  // Y1+V --gy--> VP (Y1--, V--, VP++)
    mu[i][6] = b1 * comp[i].R * comp[i].Y1; //R+Y1 --b1--> Y2 (Y1--,Y2++)
    mu[i][7] = hy * comp[i].Y1;  // Y1 --hy--> W (Y1--, W++)
    mu[i][8] = comp[i].a2 * comp[i].X2; // X2 --a2 --> 2X1 (X2--,X1+2)
    mu[i][9] = comp[i].b2 * comp[i].Y2; // Y2 --b2 --> 2Y1 (Y2--,Y1+2)
  }
  
  for(i = 0; i < N; i++)
    for(l = 0; l < NEQ; l++)
      totprop += mu[i][l];
  
  r = totprop * randd();

  for(i = 0; i < N; i++)
  {
    for(l = 0; l < NEQ; l++)
    {
      prop += mu[i][l];
      if(prop > r)
        goto jump_here;
    }
  }

  jump_here:
  
  if(l == 0)
  {
    comp[i].X1--;
    comp[i].X2++;
  }
  else if(l == 1)
  {
    comp[i].X1--;
    comp[i].V++;
  }
  else if(l == 2)
  {
    comp[i].X1--;
    comp[i].W--;
    comp[i].WP++;
  }
  else if(l == 3)
  {
    comp[i].X1--;
    comp[i].Y1++;
  }
  else if(l == 4)
  {
    comp[i].X1++;    
    comp[i].Y1--;
  }
  else if(l == 5)
  {
    comp[i].Y1--;
    comp[i].V--;
    comp[i].VP++;
  }
  else if(l == 6)
  {
    comp[i].Y1--;
    comp[i].Y2++;
  }
  else if(l == 7)
  {
    comp[i].Y1--;
    comp[i].W++;
  }
  else if(l == 8)
  {
    comp[i].X2--;
    comp[i].X1 = comp[i].X1+2;
  }
  else if(l == 9)
  {
    comp[i].Y2--;
    comp[i].Y1 = comp[i].Y1+2 ;
  }
  else
  {
    printf("ERROR - 1\n");
    printf("totprop=%lf\tprop=%lf\tr=%lf\tkk=%dreac_count=%d\n", totprop, prop, r, i, sim);
    exit(1);    
  }

  //updating the autocatalytic rates in the compartments based on the new within cell molecular composition
  comp[i].a2 = aa2*env(i);
  comp[i].b2 = bb2*env(i);
  
  time += 1.0/totprop * log(1.0/randd());
  
  return(i);
  
}

void divide(int k) 
{
  int i, order[N-1]={}, nn;
  
  for(i = 0; i < D; i++)
  {
    daughter[i].X1 = 0;
    daughter[i].X2 = 0;
    daughter[i].Y1 = 0;
    daughter[i].Y2 = 0;
    daughter[i].V = 0;
    daughter[i].W = 0;
    daughter[i].VP = 0;
    daughter[i].WP = 0;
  }

  // randomly distributing the molecules between the daughter compartments   
  for(i = 0; i < comp[k].X1; i++)
    daughter[randl(D)].X1++;
  for(i = 0; i < comp[k].X2; i++)
    daughter[randl(D)].X2++;
  for(i = 0; i < comp[k].Y1; i++)
    daughter[randl(D)].Y1++;
  for(i = 0; i < comp[k].Y2; i++)
    daughter[randl(D)].Y2++;
  for(i = 0; i < comp[k].V; i++)
    daughter[randl(D)].V++;
  for(i = 0; i < comp[k].W; i++)
    daughter[randl(D)].W++;
  for(i = 0; i < comp[k].VP; i++)
    daughter[randl(D)].VP++;
  for(i = 0; i < comp[k].WP; i++)
    daughter[randl(D)].WP++;
  
  // overwriting the parent
  comp[k].X1 = daughter[0].X1;
  comp[k].X2 = daughter[0].X2;
  comp[k].Y1 = daughter[0].Y1;
   comp[k].Y2 = daughter[0].Y2;
  comp[k].V = daughter[0].V;
  comp[k].W = daughter[0].W;
  comp[k].VP = daughter[0].VP;
  comp[k].WP = daughter[0].WP;
  comp[k].R = R0;
  comp[k].a2 = aa2*env(k);
  comp[k].b2 = bb2*env(k);
  
  // shuffling the order of the compartments
  nn = 0;
  for(i = 0; i < N; i++)
    if(i != k)
      order[nn++] = i;

  shuffle(order, N-1);
  
  // overwriting OWR randomly selected compartments
  for(i = 0; i < OWR; i++)
  {
    comp[order[i]].X1 = daughter[i + 1].X1;
    comp[order[i]].X2 = daughter[i + 1].X2;
    comp[order[i]].Y1 = daughter[i + 1].Y1;
    comp[order[i]].Y2 = daughter[i + 1].Y2;
    comp[order[i]].V = daughter[i + 1].V;
    comp[order[i]].W = daughter[i + 1].W;
    comp[order[i]].VP = daughter[i + 1].VP;
    comp[order[i]].WP = daughter[i + 1].WP;
    comp[order[i]].R = R0;
    comp[order[i]].a2 = aa2 * env(order[i]);
    comp[order[i]].b2 = bb2 *   env(order[i]);
  }
  numdiv++;

}

// data for evalution and plotting
void statistics()
{
  int k, xdom = 0, ydom = 0, totx1 = 0, totx2 = 0, toty1 = 0, toty2 = 0, totv = 0, totw = 0, totvp = 0, totwp = 0;
  double avea2 = 0.0, aveb2 = 0.0;
  
  for(k = 0; k < N; k++)
  {
    avea2 += comp[k].a2;
    aveb2 += comp[k].b2;
    totx1 += comp[k].X1;
    totx2 += comp[k].X2;
    toty1 += comp[k].Y1;
    toty2 += comp[k].Y2;
    totv += comp[k].V;
    totw += comp[k].W;
    totvp += comp[k].VP;
    totwp += comp[k].WP;
    if(comp[k].X1+comp[k].X2 > comp[k].Y1+comp[k].Y2) xdom++;
    if(comp[k].Y1+comp[k].Y2 > comp[k].X1+comp[k].X2) ydom++;
  }
  printf("%lf\t%d\t%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\n",time, sim, theta, totx1, totx2, toty1, toty2, totv, totw, totvp, totwp, xdom, ydom, numdiv, avea2/N, aveb2/N);
  fflush(stdout);
}



int main(void)
{
  int i, ind=0;
  
  seed(13541);
  theta = 1.0;
  init();

  
  for(sim = 0; sim < 1.2E10; sim++) 
  {
    i = Gillespie_step();
    if(comp[i].X1 + comp[i].X2 + comp[i].Y1 + comp[i].Y2 + comp[i].V + comp[i].W + comp[i].VP + comp[i].WP >= DIVSIZE) // checking if the total number of molecules reaches the division threshold
      divide(i);
    if(!(sim%2000))
       statistics();

    // changing the environment periodically
    if(time > 100 && ind == 0)
    {
       ind = 1;
       theta *= -1.0;
    }
    
    if(time > 200 && ind == 1)
    {
       ind = 2;
       theta *= -1.0;
    }

    if(time > 300 && ind == 2)
    {
       ind = 3;
       theta *= -1.0;
    }

    if(time > 400 && ind == 3)
    {
       ind = 4;
       theta *= -1.0;
    }

    if(time > 500 && ind == 4)
    {
       ind = 5;
       theta *= -1.0; 
    }
    if(time > 600 && ind == 5)
    {
       ind = 6;
       theta *= -1.0; 
    }
    if(time > 700 && ind == 6)
    {
       ind = 7;
       theta *= -1.0; 
    }
    if(time > 800 && ind == 7)
    {
       ind = 8;
       theta *= -1.0; 
    }
    if(time > 900 && ind == 8)
    {
      exit(1);
    }
  }
  return 0;
}
