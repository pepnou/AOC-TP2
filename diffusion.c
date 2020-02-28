/*
  Gray Scott Algorithm for fluid diffusion simulation

  References
  -----------

  https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system
  https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/
  https://mrob.com/pub/comp/xmorphia/

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

//For graphics
#ifdef WITH_X11
#include "flame.h"
#endif

//
#include "rdtsc.h"

//
#define max(a, b) ((a) > (b)) ? (a) : (b)

//
#define GRID_H 201   //Y axis
#define GRID_W 201   //X axis
#define GRID_C 2     //Chemicals: A & B

//For pixel RGB channels
typedef unsigned char byte;

//Laplacian - Stencil
float laplacian(float g[GRID_H][GRID_W][GRID_C], int i, int j, int c)
{
   float sum = -g[i][j][c];

   sum += (g[i - 1][j][c]     *   0.2);
   sum += (g[i + 1][j][c]     *   0.2);
   sum += (g[i][j - 1][c]     *   0.2);
   sum += (g[i][j + 1][c]     *   0.2);
   sum += (g[i - 1][j - 1][c] *   0.05);
   sum += (g[i - 1][j + 1][c] *   0.05);
   sum += (g[i + 1][j - 1][c] *   0.05);
   sum += (g[i + 1][j + 1][c] *   0.05);

   return sum;
}

//
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
}

//
void diffuse(
#ifdef WITH_X11
             flame_obj_t *fo,
#endif
             int bx, int by,
             float G[GRID_H][GRID_W][GRID_C],
             float G1[GRID_H][GRID_W][GRID_C],
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
               G[i][j][0] = G1[i][j][0];
               G[i][j][1] = G1[i][j][1];

#ifdef WITH_X11
               if ((iter % 50) == 0)
                  {
                     ta = G[i][j][0] - G1[i][j][1];

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

   //
   for (int i = 1; i < GRID_H - 1; i++) {
     for (int j = 1; j < GRID_W - 1; j++)
     {
       ta = G[i][j][0];
       tb = G[i][j][1];

       //Chemical A
       G1[i][j][0] = ta + (diff_rate_A * laplacian(G, i, j, 0) - (ta * tb * tb) + (f * (1 - ta))) * dt;

       //Chemical B
       G1[i][j][1] = tb + (diff_rate_B * laplacian(G, i, j, 1) + (ta * tb * tb) - ((k + f) * tb)) * dt;
     }
   }
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
   float G[GRID_H][GRID_W][GRID_C], G1[GRID_H][GRID_W][GRID_C];

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

#ifdef WITH_X11
   //
   flame_close(fo);
#endif

   //
   return 0;
}
