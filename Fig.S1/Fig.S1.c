#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N (1000) // number of compartments
#define D (2) // number of daughter compartments after division
#define NEQ (14) // number of reactions
#define OWR (1) // besides one daughter replacing the parent, OWR daughter compartments overwrite OWR randomly selected compartments in the population (OWR <= D - 1)
#define DIVSIZE (100) // total molecule count threshold for compartment division

#define aa (1.5) 
#define bb (1.2)
#define gx (0.2) 
#define gy (0.1)
#define gxr (0.02)
#define gyr (0.01)
#define hx (0.8)
#define hy (0.6)
#define hxr (0.08)
#define hyr (0.06)
#define mx (0.1) 
#define my (0.1)
#define R0 (70.0) 
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
  int X, Y, V, W, VP, WP; // molecule types within the compartments
  double R;
  double a;
  double ar;
  double b;
  double br;
}COMP;
COMP comp[N];

typedef struct _IND2 // structure of the daughter compartments
{
  int X, Y, V, W, VP, WP;
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
    if(comp[k].X+comp[k].Y==0){
      double multi = 1;
      return(multi);
    }
    double multi = 1+c*((double)comp[k].X/((double)comp[k].X+(double)comp[k].Y));
    return(multi);
  }
  // CASE 2: environment E_Y
  else if(theta==-1.0){
    if(comp[k].X+comp[k].Y==0){
      double multi = 1;
      return(multi);
    }
    double multi = 1+c*((double)comp[k].Y/((double)comp[k].X+(double)comp[k].Y));
    return(multi);
  }
  else
    exit(1);
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
   comp[i].X = randl(10)+1;
   comp[i].Y = randl(10)+1;
   comp[i].R = R0;
   comp[i].V = 0;
   comp[i].W = 0;
   comp[i].VP = 0;
   comp[i].WP = 0;
   comp[i].a = aa * env(i);
   comp[i].ar = comp[i].a/10;
   comp[i].b = bb * env(i);
   comp[i].br = comp[i].b/10;
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
    mu[i][0] = comp[i].a * comp[i].R * comp[i].X; // R+X --a--> 2X (X++)
    mu[i][1] = hx * comp[i].X;  // X --hx--> V (X--, V++)
    mu[i][2] = gx * comp[i].X * comp[i].W;  // X+W --gx--> WP (X--, W--, WP++)
    mu[i][3] = mx * comp[i].X;  // X --mx--> Y (X--, Y++)
    mu[i][4] = my * comp[i].Y;  // Y --my--> X (Y--, X++)
    mu[i][5] = gy * comp[i].Y * comp[i].V;  // Y+V --gy--> VP (Y--, V--, VP++)
    mu[i][6] = comp[i].b * comp[i].R * comp[i].Y; //R+Y --b--> 2Y (Y++)
    mu[i][7] = hy * comp[i].Y;  // Y --hy--> W (Y--, W++)
    mu[i][8] = comp[i].ar * comp[i].X * comp[i].X;// X+X --ar--> X+R (X--)
    mu[i][9] = hxr * comp[i].V;  // V --hxr--> X (V--, X++)
    mu[i][10] = gxr * comp[i].WP;  // WP --gxr--> X + W (X++, W++, WP--)
    mu[i][11] = gyr * comp[i].VP;  // CVP --gyr--> Y + V (Y++, V++, VP--)
    mu[i][12] = comp[i].br * comp[i].Y * comp[i].Y; //2Y --br--> Y + R (Y--)
    mu[i][13] = hyr * comp[i].W;  // W --hyr--> Y (Y++, W--)
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
    comp[i].X++;
  }
  else if(l == 1)
  {
    comp[i].X--;
    comp[i].V++;
  }
  else if(l == 2)
  {
    comp[i].X--;
    comp[i].W--;
    comp[i].WP++;
  }
  else if(l == 3)
  {
    comp[i].X--;
    comp[i].Y++;
  }
  else if(l == 4)
  {
    comp[i].X++;    
    comp[i].Y--;
  }
  else if(l == 5)
  {
    comp[i].Y--;
    comp[i].V--;
    comp[i].VP++;
  }
  else if(l == 6)
  {
    comp[i].Y++;
  }
  else if(l == 7)
  {
    comp[i].Y--;
    comp[i].W++;
  }
  else if(l == 8)
  {
    comp[i].X--;
  }
  else if(l == 9)
  {
    comp[i].X++;
    comp[i].V--;
  }
  else if(l == 10)
  {
    comp[i].X++;
    comp[i].W++;
    comp[i].WP--;
  }
  else if(l == 11)
  {
    comp[i].Y++;
    comp[i].V++;
    comp[i].VP--;
  }
  else if(l == 12)
  {
    comp[i].Y--;
  }
  else if(l == 13)
  {
    comp[i].Y++;
    comp[i].W--;
  }
  else
  {
    printf("ERROR - 1\n");
    printf("totprop=%lf\tprop=%lf\tr=%lf\tkk=%dreac_count=%d\n", totprop, prop, r, i, sim);
    exit(1);    
  }

  //updating the autocatalytic rates in the compartments based on the new within cell molecular composition
  comp[i].a = aa*env(i);
  comp[i].b = bb*env(i);
  comp[i].ar = comp[i].a/10;
  comp[i].br = comp[i].b/10;
  
  time += 1.0/totprop * log(1.0/randd());
  
  return(i);
  
}

void divide(int k)
{
  int i, order[N-1]={}, nn;

  for(i = 0; i < D; i++)
  {
    daughter[i].X = 0;
    daughter[i].Y = 0;
    daughter[i].V = 0;
    daughter[i].W = 0;
    daughter[i].VP = 0;
    daughter[i].WP = 0;
  }

  // randomly distributing the molecules between the daughter compartments
  for(i = 0; i < comp[k].X; i++)
    daughter[randl(D)].X++;
  for(i = 0; i < comp[k].Y; i++)
    daughter[randl(D)].Y++;
  for(i = 0; i < comp[k].V; i++)
    daughter[randl(D)].V++;
  for(i = 0; i < comp[k].W; i++)
    daughter[randl(D)].W++;
  for(i = 0; i < comp[k].VP; i++)
    daughter[randl(D)].VP++;
  for(i = 0; i < comp[k].WP; i++)
    daughter[randl(D)].WP++;
  
  // overwriting the parent
  comp[k].X = daughter[0].X;
  comp[k].Y = daughter[0].Y;
  comp[k].V = daughter[0].V;
  comp[k].W = daughter[0].W;
  comp[k].VP = daughter[0].VP;
  comp[k].WP = daughter[0].WP;
  comp[k].R = R0;
  comp[k].a = aa*env(k);
  comp[k].b = bb*env(k);
  comp[k].ar = comp[k].a/10;
  comp[k].br = comp[k].b/10;
  
  // shuffling the order of the compartments
  nn = 0;
  for(i = 0; i < N; i++)
    if(i != k)
      order[nn++] = i;

  shuffle(order, N-1);
  
  // overwriting OWR randomly selected compartments 
  for(i = 0; i < OWR; i++)
  {
    comp[order[i]].X = daughter[i + 1].X;
    comp[order[i]].Y = daughter[i + 1].Y;
    comp[order[i]].V = daughter[i + 1].V;
    comp[order[i]].W = daughter[i + 1].W;
    comp[order[i]].VP = daughter[i + 1].VP;
    comp[order[i]].WP = daughter[i + 1].WP;
    comp[order[i]].R = R0;
    comp[order[i]].a = aa * env(order[i]);
    comp[order[i]].b = bb *   env(order[i]);
    comp[order[i]].ar = comp[order[i]].a/10;
    comp[order[i]].br = comp[order[i]].b/10;      
  }
  numdiv++;

}

// data for evalution and plotting
void statistics()
{
  int k, xdom = 0, ydom = 0, totx = 0, toty = 0, totv = 0, totw = 0, totvp = 0, totwp = 0;
  double avea = 0.0, aveb = 0.0, avear = 0.0, avebr = 0.0;
  
  for(k = 0; k < N; k++)
  {
    avea += comp[k].a;
    aveb += comp[k].b;
    avear += comp[k].ar;
    avebr += comp[k].br;
    totx += comp[k].X;
    toty += comp[k].Y;
    totv += comp[k].V;
    totw += comp[k].W;
    totvp += comp[k].VP;
    totwp += comp[k].WP;
    if(comp[k].X > comp[k].Y) xdom++;
    if(comp[k].Y > comp[k].X) ydom++;
  }
  printf("%lf\t%d\t%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n",time, sim, theta, totx, toty, totv, totw, totvp, totwp, xdom, ydom, numdiv, avea/N, aveb/N, avear/N, avebr/N);
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
    if(comp[i].X + comp[i].Y + comp[i].V + comp[i].W + comp[i].VP + comp[i].WP >= DIVSIZE) // checking if the total number of molecules reaches the division threshold
      divide(i);
    if(!(sim%2000))
       statistics();

    // changing the environment periodically
    if(time > 1 && ind == 0)
    {
       ind = 1;
       theta *= -1.0;
    }
    
    if(time > 2 && ind == 1)
    {
       ind = 2;
       theta *= -1.0;
    }

    if(time > 3 && ind == 2)
    {
       ind = 3;
       theta *= -1.0;
    }

    if(time > 4 && ind == 3)
    {
       ind = 4;
       theta *= -1.0;
    }

    if(time > 5 && ind == 4)
    {
       ind = 5;
       theta *= -1.0; 
    }

    if(time > 6 && ind == 5)
    {
       ind = 6;
       theta *= -1.0;
    }

    if(time > 7 && ind == 6)
    {
       ind = 7;
       theta *= -1.0;
    }

    if(time > 8 && ind == 7)
    {
       ind = 8;
       theta *= -1.0; 
    }
    
    if(time > 9)
    {
      exit(1);
    }
  }
  return 0;
}
