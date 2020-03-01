/*
   Gray Scott Algorithm for fluid diffusion simulation

   References
   -----------

https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system
https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/
https://mrob.com/pub/comp/xmorphia/

*/

#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include <pthread.h>

//For graphics
#ifdef WITH_X11
#include "flame.h"
#endif

//
#include "rdtsc.h"

//
#define max(a, b) ((a) > (b)) ? (a) : (b)
#define min(a, b) ((a) < (b)) ? (a) : (b)

//
#define GRID_H 201   //Y axis
#define GRID_W 201   //X axis
#define GRID_C 2     //Chemicals: A & B

#define index(A,i,j,c) A[(c)*GRID_H*GRID_W + (i)*GRID_W + (j)]

int modifA = 0;
int modifB = 0;

//For pixel RGB channels
typedef unsigned char byte;

//Laplacian - Stencil
/*float laplacian(float g[GRID_H][GRID_W][GRID_C], int i, int j, int c)
  {
  float sum = 0.0;

/*sum += (g[i - 1][j - 1][c] *   0.05);
sum += (g[i - 1][j][c]     *   0.2);
sum += (g[i - 1][j + 1][c] *   0.05);

sum += (g[i][j - 1][c]     *   0.2);
sum -= (g[i][j][c]);
sum += (g[i][j + 1][c]     *   0.2);

sum += (g[i + 1][j - 1][c] *   0.05);
sum += (g[i + 1][j][c]     *   0.2);
sum += (g[i + 1][j + 1][c] *   0.05);

sum += (g[i - 1][j - 1][c] + g[i - 1][j + 1][c] + g[i + 1][j - 1][c] + g[i + 1][j + 1][c]) * 0.05;
sum += (g[i - 1][j][c] + g[i][j-1][c] + g[i][j + 1][c] + g[i+1][j][c]) * 0.2;
sum -= (g[i][j][c]);

return sum;
}*/

//Laplacian - Stencil
float laplacian(float * restrict g, int i, int j, int c)
{
  float sum = 0.0;

  /*sum += (index(g,i - 1,j - 1,c) *   0.05);
    sum += (index(g,i - 1,j,c)     *   0.2);
    sum += (index(g,i - 1,j + 1,c) *   0.05);

    sum += (index(g,i,j - 1,c)     *   0.2);
    sum -= (index(g,i,j,c));
    sum += (index(g,i,j + 1,c)     *   0.2);

    sum += (index(g,i + 1,j - 1,c) *   0.05);
    sum += (index(g,i + 1,j,c)     *   0.2);
    sum += (index(g,i + 1,j + 1,c) *   0.05);
    */

  sum += (index(g,i - 1,j - 1,c) + index(g,i - 1,j + 1,c) + index(g,i + 1,j - 1,c) + index(g,i + 1,j + 1,c)) * 0.05;
  sum += (index(g,i - 1,j,c) + index(g,i,j-1,c) + index(g,i,j + 1,c) + index(g,i+1,j,c)) * 0.2;
  sum -= (index(g,i,j,c));

  return sum;
}

//Laplacian - Stencil
void laplacian2(float * g, float * g1, int i, int j, int c, int max_j, int chem)
{
  const float diff_rate_A = 1.0, diff_rate_B = 0.5;

  //Coral Growth, f = 0.0545
  //Mitosis, f = 0.0545
  const float f = 0.0321, k = 0.062, dt = 1.1; //f: feed, k: kill, dt: delta time

  float * gtemp = &index(g, i, j-1, c);
  float * res = &index(g1, i, j, c);

  int l;
  memset(res, 0, 14 * sizeof(float));



  res[0] += (gtemp[-1] * 0.05);
  res[0] += (gtemp[0]  * 0.2);
  res[1] += (gtemp[0]  * 0.05);

  for(l = 1; l < 13 && l + j < max_j; l++) {
    res[l-1] += (gtemp[l] * 0.05);
    res[l]   += (gtemp[l] * 0.2);
    res[l+1] += (gtemp[l] * 0.05);
  }

  res[l-1] += (gtemp[l]   * 0.05);
  res[l]   += (gtemp[l]   * 0.2);
  res[l]   += (gtemp[l+1] * 0.05);


  gtemp += GRID_W;


  res[0] += (gtemp[-1] * 0.2);
  res[0] -= gtemp[0];
  res[1] += (gtemp[0]  * 0.2);

  for(l = 1; l < 13 && l + j < max_j; l++) {
    res[l-1] += (gtemp[l] * 0.2);
    res[l]   -= gtemp[l];
    res[l+1] += (gtemp[l] * 0.2);
  }

  res[l-1] += (gtemp[l]   * 0.2);
  res[l]   -= gtemp[l];
  res[l]   += (gtemp[l+1] * 0.2);


  gtemp += GRID_W;


  res[0] += (gtemp[-1] * 0.05);
  res[0] += (gtemp[0]  * 0.2);
  res[1] += (gtemp[0]  * 0.05);

  for(l = 1; l < 13 && l + j < max_j; l++) {
    res[l-1] += (gtemp[l] * 0.05);
    res[l]   += (gtemp[l] * 0.2);
    res[l+1] += (gtemp[l] * 0.05);
  }

  res[l-1] += (gtemp[l]   * 0.05);
  res[l]   += (gtemp[l]   * 0.2);
  res[l]   += (gtemp[l+1] * 0.05);


  float * dat0 = &index(g, i, j, 0);
  float * dat1 = &index(g, i, j, 1);
  float ta, tb;

  switch(chem) {
    case 0: // Chem A
      {
        for(l = 0; l < 14 && l + j < max_j; l++) {
          ta = dat0[l];
          tb = dat1[l];

          //index(G1,i,j,0) = ta + (diff_rate_A * laplacian(G, i, j, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;
          res[l] = ta + (diff_rate_A * res[l] - (ta * tb * tb) + (f * (1 - ta))) * dt;
        }
        break;
      }
    case 1: // Chem B
      {
        for(l = 0; l < 14 && l + j < max_j; l++) {
          ta = dat0[l];
          tb = dat1[l];

          //index(G1,i,j,1) = tb + (diff_rate_B * laplacian(G, i, j, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;
          res[l] = tb + (diff_rate_B * res[l] + (ta * tb * tb) - ((k + f) * tb)) * dt;
        }
        break;
      }
  }
}

/*void convol(float G[GRID_H][GRID_W][GRID_C], float G1[GRID_H][GRID_W][GRID_C], int iter) {
  float ta, tb;
  const float diff_rate_A = 1.0, diff_rate_B = 0.5;

//Coral Growth, f = 0.0545
//Mitosis, f = 0.0545
const float f = 0.0321, k = 0.062, dt = 1.1; //f: feed, k: kill, dt: delta time

int min_i = (GRID_H >> 1) - 10 - iter, max_i = min(GRID_H - 1, (GRID_H >> 1) + 10 + iter);
if(min_i < 1)
min_i = 1;

int min_j = (GRID_W >> 1) - 10 - iter, max_j = min(GRID_W - 1, (GRID_W >> 1) + 10 + iter);
if(min_j < 1)
min_j = 1;

for (int i = min_i; i < max_i; i++) {
for (int j = min_j; j < max_j; j++) {
ta = G[i][j][0];
tb = G[i][j][1];

//Chemical A
G1[i][j][0] = ta + (diff_rate_A * laplacian(G, i, j, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;
}
}

for (int i = min_i; i < max_i; i++) {
for (int j = min_j; j < max_j; j++) {
ta = G[i][j][0];
tb = G[i][j][1];

//Chemical B
G1[i][j][1] = tb + (diff_rate_B * laplacian(G, i, j, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;
}
}
}*/


void convol(float * restrict G, float * restrict G1, int iter) {
  float ta, tb;
  const float diff_rate_A = 1.0, diff_rate_B = 0.5;

  //Coral Growth, f = 0.0545
  //Mitosis, f = 0.0545
  const float f = 0.0321, k = 0.062, dt = 1.1; //f: feed, k: kill, dt: delta time

  int min_i = (GRID_H >> 1) - 10 - iter, max_i = min(GRID_H - 1, (GRID_H >> 1) + 10 + iter);
  if(min_i < 1)
    min_i = 1;

  int min_j = (GRID_W >> 1) - 10 - iter, max_j = min(GRID_W - 1, (GRID_W >> 1) + 10 + iter);
  if(min_j < 1)
    min_j = 1;

  float * g = __builtin_assume_aligned (G, 32);

#pragma omp parallel num_threads(2)
  {
#pragma omp single nowait
    {
      for (int i = min_i; i < max_i; i++) {
        for (int j = min_j; j < max_j; j++) {
          ta = index(g,i,j,0);
          tb = index(g,i,j,1);

          //Chemical A
          index(G1,i,j,0) = ta + (diff_rate_A * laplacian(g, i, j, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;
        }
      }
    }

#pragma omp single nowait
    {
      for (int i = min_i; i < max_i; i++) {
        for (int j = min_j; j < max_j; j++) {
          ta = index(g,i,j,0);
          tb = index(g,i,j,1);

          //Chemical B
          index(G1,i,j,1) = tb + (diff_rate_B * laplacian(g, i, j, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;
        }
      }
    }
  }
}

typedef struct {
  float *G, *G1;
  int iter;
  int chem;
} Args;

void * convol4thread(void * a) {
  Args *args = (Args*)a;

  float ta, tb;
  const float diff_rate_A = 1.0, diff_rate_B = 0.5;

  //Coral Growth, f = 0.0545
  //Mitosis, f = 0.0545
  const float f = 0.0321, k = 0.062, dt = 1.1; //f: feed, k: kill, dt: delta time

  int min_i = (GRID_H >> 1) - 10 - args->iter, max_i = min(GRID_H - 1, (GRID_H >> 1) + 10 + args->iter);
  if(min_i < 1)
    min_i = 1;

  int min_j = (GRID_W >> 1) - 10 - args->iter, max_j = min(GRID_W - 1, (GRID_W >> 1) + 10 + args->iter);
  if(min_j < 1)
    min_j = 1;

  float * g = __builtin_assume_aligned (args->G, 32);

  if(args->chem == 0) {

    for (int i = min_i; i < max_i; i++) {
      for (int j = min_j; j < max_j; j++) {
        ta = index(g,i,j,0);
        tb = index(g,i,j,1);

        //Chemical A
        index(args->G1,i,j,0) = ta + (diff_rate_A * laplacian(g, i, j, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;
      }
    }
  } else {

    for (int i = min_i; i < max_i; i++) {
      for (int j = min_j; j < max_j; j++) {
        ta = index(g,i,j,0);
        tb = index(g,i,j,1);

        //Chemical B
        index(args->G1,i,j,1) = tb + (diff_rate_B * laplacian(g, i, j, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;
      }
    }
  }

  return NULL;
}

void convol4(float * restrict G, float * restrict G1, int iter, int chem) {

  float ta, tb;
  const float diff_rate_A = 1.0, diff_rate_B = 0.5;

  //Coral Growth, f = 0.0545
  //Mitosis, f = 0.0545
  const float f = 0.0321, k = 0.062, dt = 1.1; //f: feed, k: kill, dt: delta time

  int min_i = (GRID_H >> 1) - 10 - iter, max_i = min(GRID_H - 1, (GRID_H >> 1) + 10 + iter);
  if(min_i < 1)
    min_i = 1;

  int min_j = (GRID_W >> 1) - 10 - iter, max_j = min(GRID_W - 1, (GRID_W >> 1) + 10 + iter);
  if(min_j < 1)
    min_j = 1;

  float * g = __builtin_assume_aligned (G, 32);

  if(chem == 0) {

    for (int i = min_i; i < max_i; i++) {
      for (int j = min_j; j < max_j; j++) {
        ta = index(g,i,j,0);
        tb = index(g,i,j,1);

        //Chemical A
        index(G1,i,j,0) = ta + (diff_rate_A * laplacian(g, i, j, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;
      }
    }
  } else {

    for (int i = min_i; i < max_i; i++) {
      for (int j = min_j; j < max_j; j++) {
        ta = index(g,i,j,0);
        tb = index(g,i,j,1);

        //Chemical B
        index(G1,i,j,1) = tb + (diff_rate_B * laplacian(g, i, j, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;
      }
    }
  }
}







/*
   void convol2(float * G, float * G1, int iter) {
   float ta, tb;
   const float diff_rate_A = 1.0, diff_rate_B = 0.5;

//Coral Growth, f = 0.0545
//Mitosis, f = 0.0545
const float f = 0.0321, k = 0.062, dt = 1.1; //f: feed, k: kill, dt: delta time

int min_i = (GRID_H >> 1) - 10 - iter, max_i = min(GRID_H - 1, (GRID_H >> 1) + 10 + iter);
if(min_i < 1)
min_i = 1;

int min_j = (GRID_W >> 1) - 10 - iter, max_j = min(GRID_W - 1, (GRID_W >> 1) + 10 + iter);
if(min_j < 1)
min_j = 1;

float add;

for (int i = min_i; i < max_i; i++) {
for (int j = min_j; j < max_j; j += 14) {
//ta = index(G,i,j,0);
//tb = index(G,i,j,1);

//Chemical A
//index(G1,i,j,0) = ta + (diff_rate_A * laplacian(G, i, j, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;

laplacian2(G, G1, i, j, 0, max_j, 0);
}
}

for (int i = min_i; i < max_i; i++) {
for (int j = min_j; j < max_j; j += 14) {
//ta = index(G,i,j,0);
//tb = index(G,i,j,1);

//Chemical B
//index(G1,i,j,1) = tb + (diff_rate_B * laplacian(G, i, j, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;

laplacian2(G, G1, i, j, 1, max_j, 1);
}
}
}

void convol3(float * G, float * G1, int iter) {
float ta, tb;
const float diff_rate_A = 1.0, diff_rate_B = 0.5;

//Coral Growth, f = 0.0545
//Mitosis, f = 0.0545
const float f = 0.0321, k = 0.062, dt = 1.1; //f: feed, k: kill, dt: delta time

int min_i = (GRID_H >> 1) - 10 - iter, max_i = min(GRID_H - 1, (GRID_H >> 1) + 10 + iter);
if(min_i < 1) {
min_i = 1;
}

int min_j = (GRID_W >> 1) - 10 - iter, max_j = min(GRID_W - 1, (GRID_W >> 1) + 10 + iter);
if(min_j < 1) {
min_j = 1;
}

for (int j = min_j; j < max_j; j += 16) {
for (int i = min_i; i < max_i; i++) {
for(int l = 0; l < 16 && l + j < max_j; l++) {
ta = index(G,i,j+l,0);
tb = index(G,i,j+l,1);

//Chemical A
index(G1,i,j+l,0) = ta + (diff_rate_A * laplacian(G, i, j+l, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;
}
}
}

for (int j = min_j; j < max_j; j++) {
  for (int i = min_i; i < max_i; i++) {
    for(int l = 0; l < 16 && l+j < max_j; l++) {
      ta = index(G,i,j+l,0);
      tb = index(G,i,j+l,1);

      //Chemical B
      index(G1,i,j+l,1) = tb + (diff_rate_B * laplacian(G, i, j+l, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;
    }
  }
}
}*/

/*//
  void init(float G1[GRID_H][GRID_W][GRID_C])
  {
//Initialize grid with chemical A only
for (int i = 1; i < GRID_H - 1; i++)
for (int j = 1; j < GRID_W - 1; j++)
{
G1[i][j][0] = 1.0; //Chemical A
G1[i][j][1] = 0.0;
}

unsigned long long n = (GRID_H >> 1);

//Add some of chemical B
for (int i = n - 5; i < n + 5; i++)
for (int j = n - 5; j < n + 5; j++)
G1[i][j][1] = 1.0; //Chemical B
}*/

void init(float* restrict G1)
{
  //Initialize grid with chemical A only
  for (int i = 1; i < GRID_H - 1; i++)
    for (int j = 1; j < GRID_W - 1; j++)
    {
      index(G1,i,j,0) = 1.0; //Chemical A
      index(G1,i,j,1) = 0.0;
    }

  unsigned long long n = (GRID_H >> 1);

  //Add some of chemical B
  for (int i = n - 5; i < n + 5; i++)
    for (int j = n - 5; j < n + 5; j++)
      index(G1,i,j,1) = 1.0; //Chemical B
}



//
void diffuse(
#ifdef WITH_X11
    flame_obj_t *fo,
#endif
    int bx, int by,
    //float G[GRID_H][GRID_W][GRID_C],
    float * restrict G,
    //float G1[GRID_H][GRID_W][GRID_C],
    float * restrict G1,
#ifdef WITH_X11
    unsigned frame[GRID_H][GRID_W],
#endif
    unsigned long long iter)
{
  int color;
  float ta, tb;
  const float diff_rate_A = 1.0, diff_rate_B = 0.5;

  //Coral Growth, f = 0.0545
  //Mitosis, f = 0.0545
  const float f = 0.0321, k = 0.062, dt = 1.1; //f: feed, k: kill, dt: delta time

  //
  for (int i = 1; i < GRID_H - 1; i++)
  {
    for (int j = 1; j < GRID_W - 1; j++)
    {
      //G[i][j][0] = G1[i][j][0];
      //G[i][j][1] = G1[i][j][1];
      index(G,i,j,0) = index(G1,i,j,0);
      index(G,i,j,1) = index(G1,i,j,1);
#ifdef WITH_X11
      if ((iter % 50) == 0)
      {
        //ta = G[i][j][0] - G1[i][j][1];
        ta = index(G,i,j,0) - index(G1,i,j,1);

        //Set r, g, b according to temp
        if (ta > 0.8)
          color = 0xFFFFFF00; //White
        else
          if (ta > 0.5)
            color = 0xBEBEBE00; //Gray
          else
            if (ta > 0.1)
              color = 0xA9A9A900; //Dark Gray
            else
              color = 0; //Black

        frame[i][j] = color;
      }
#endif
    }
  }
  //printf("A : %d\n", modifA * 100 / (GRID_H * GRID_W));
  //printf("B : %d\n", modifB * 100 / (GRID_H * GRID_W));

  /*int min_i = (GRID_H >> 1) - 10 - iter, max_i = min(GRID_H - 1, (GRID_H >> 1) + 10 + iter);
    if(min_i < 1)
    min_i = 1;

    int min_j = (GRID_W >> 1) - 10 - iter, max_j = min(GRID_W - 1, (GRID_W >> 1) + 10 + iter);
    if(min_j < 1)
    min_j = 1;

    float* G_align = __builtin_assume_aligned (G, 32);

  //
  for (int j = min_j; j < max_j; j += 16) {
  for (int i = min_i; i < max_i; i++) {
  for (int l = 0; l+j  < max_j && l < 16; l++) {
  //ta = G[i][j+l][0];
  //tb = G[i][j+l][1];

  ta = index(G,i,j+l,0);
  tb = index(G,i,j+l,1);

  //Chemical A
  index(G1,i,j+l,0) = ta + (diff_rate_A * laplacian(G_align, i, j+l, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;

  //Chemical B
  index(G1,i,j+l,1) = tb + (diff_rate_B * laplacian(G_align, i, j+l, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;
  }
  }
  }*/


  Args args;
  args.G = G;
  args.G1 = G1;
  args.iter = iter;
  args.chem = 0;

  pthread_t pid;
  pthread_create(&pid, NULL, convol4thread, &args);
  convol4(G, G1, iter, 1);
  pthread_join(pid, NULL);




  //convol(G, G1, iter);




  /*for (int i = min_i; i < max_i; i++) {
    for (int j = min_j; j < max_j; j++)
    {
  //ta = G[i][j][0];
  //tb = G[i][j][1];

  ta = index(G,i,j,0);
  tb = index(G,i,j,1);
  //Chemical A
  index(G1,i,j,0) = ta + (diff_rate_A * laplacian(G, i, j, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;

  //Chemical B
  index(G1,i,j,1) = tb + (diff_rate_B * laplacian(G, i, j, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;
  }
  }*/
}

//
int main(int argc, char **argv)
{
  if (argc < 2)
    return printf("Usage: %s [number of iterations]\n", argv[0]), 2;

  //
  char c;
  int x_min = 50, x_max = 1900;
  int y_min = 50, y_max = 1000;
#ifdef WITH_X11
  flame_obj_t *fo = flame_open("Diffusion", x_max, y_max);
#endif

  double after, before;

  int bx = 5, by = y_max >> 1;

#ifdef WITH_X11
  XEvent event;
#endif

  //
  unsigned frame[GRID_H][GRID_W]; 

  //
  //float G[GRID_H][GRID_W][GRID_C];
  //float G1[GRID_H][GRID_W][GRID_C];

  void *G, *G1;
  if(!posix_memalign (&G, 32, GRID_H * GRID_W * GRID_C * sizeof(float))) perror("posix_memalign");
  if(!posix_memalign (&G1, 32, GRID_H * GRID_W * GRID_C * sizeof(float))) perror("posix_memalign");

  //G = malloc(GRID_H * GRID_W * GRID_C * sizeof(float));
  //G1 = malloc(GRID_H * GRID_W * GRID_C * sizeof(float));
  //updated = malloc(GRID_H * GRID_W * GRID_C * sizeof(float));


  //Memory sizes
  unsigned s      = GRID_H * GRID_W;
  unsigned sG     = (sizeof(float) * s * GRID_C) >> 10;
  unsigned sG1    = (sizeof(float) * s * GRID_C) >> 10;
  unsigned sframe = (sizeof(unsigned) * s) >> 10;

  //In KiB & MiB
  unsigned stotal_k = (sG + sG1 + sframe);
  unsigned stotal_m = (sG + sG1 + sframe) >> 10;

  //
  unsigned long max_iter = atol(argv[1]);

  fprintf(stderr,
      "Size G     : %ukiB\n"
      "Size G1    : %ukiB\n"
      "Size frame : %ukiB\n"
      "Total size : %ukiB ~ %uMiB\n",
      sG, sG1, sframe, stotal_k, stotal_m);

  //Magic starts here
#ifdef WITH_X11
  flame_clear_display(fo);
#endif

  //
  init(G1);

  //
  for (unsigned long iter = 0; iter < max_iter; iter++)
  {
    //Start probe
    before = rdtsc();

    //Call the solver
#ifdef WITH_X11
    diffuse(fo, bx, by, G, G1, frame, iter);
#else
    diffuse(bx, by, G, G1, iter);
#endif

    //Stop probe
    after = rdtsc();

    //
    printf("Iteration %lu: %.2f ref-cycles / cell\n",
        iter, (double)(after - before) / ((GRID_H-2) * (GRID_W-2)));

#ifdef WITH_X11
    //Draw frame (graphics)
    for (int i = 1; i < GRID_H; i++)
      for (int j = 1; j < GRID_W; j++)
      {
        byte r = (frame[i][j] >> 24) & 0x000000FF;
        byte g = (frame[i][j] >> 16) & 0x000000FF;
        byte b = (frame[i][j] >>  8) & 0x000000FF;

        flame_set_color(fo, r, g, b);
        flame_draw_point(fo, bx + (i * 1), by - (j * 1));
      }
#endif
  }

  /*for (int i = 1; i < GRID_H - 1; i++)
  {
    for (int j = 1; j < GRID_W - 1; j++)
    {
      printf("%g ", index(((float*)G1), i, j, 0));
    }
    printf("\n");
  }


  printf("\n");
  printf("\n");


  for (int i = 1; i < GRID_H - 1; i++)
  {
    for (int j = 1; j < GRID_W - 1; j++)
    {
      printf("%g ", index(((float*)G1), i, j, 1));
    }
    printf("\n");
  }*/





#ifdef WITH_X11
  //
  flame_close(fo);
#endif

  //
  return 0;
}
