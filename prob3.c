#include <stdio.h>
#include <math.h>
#include "zekelib2.c"

int main(void)
{
      int Nobj = 2, n = 2048;
      int i, j, k;
      struct planet sol[Nobj];
      struct planet oldSol[Nobj];
      double T = 10000, dt = 1/(double)(n);
      double dx[3], d, t, E[Nobj], E0[Nobj];
      FILE *poss, *energy1, *angMom;
      poss = fopen("possition4", "wt");
      energy1 = fopen("energy", "wt");
      angMom = fopen("momentum", "wt");
      printf("dt = %f\n", dt);
      //initiate my two objects
      initSolSystem(&sol[0], Nobj);
      sol[0].v[0] = - sol[1].v[0]*sol[1].mass/sol[0].mass;
      sol[0].v[1] = - sol[1].v[1]*sol[1].mass/sol[0].mass;
      sol[0].v[2] = 0;
      sol[1].v[2] = 0;

      //this will convert the coordinates to a
      //center of mass coordinate system
      center(&sol[0], Nobj);

      //get the initial energy and angular momentum
      energy(sol, &E0[0], Nobj);      
      momentum(sol, Nobj, angMom); 

      //now start solving for the acceleration
      for(t=0.0; t<=T; t+=dt)
      {
            //save to a file every 16 time steps
            if((t<1) || ((t>2999) && (t<3000)) || ((t>5999) && (t<6000)) || (t>9999))
            {
                  //printf("do i get here? %.8e \n", sol[1].r[1]);
                  // print the energies
                  energy(sol, &E[0], Nobj);
                  momentum(sol, Nobj, angMom);
                  for(j=0; j<Nobj; j++)
                  {
                        fprintf(energy1, "%.8e  %.8e\t", E[j], E[j]-E0[j]);
                  }
                  fprintf(energy1, "\n");

                  for(j=0; j<Nobj; j++)
                  {
                        for(i=0; i<3; i++)
                        {
                              fprintf(poss, "%.8e  ", sol[j].r[i]);
                        }
                        fprintf(poss, "\t");
                  }
                  fprintf(poss, "\n");
            }
            
            //compute the acceleration
            accelEin(&sol[0], Nobj);

            //now do time steps for every planet
            if(t==0.0)
            {
                  printf("I do get here yes?\n");
                  for(j=0; j<Nobj; j++)
                  {
                        for(i=0; i<3; i++)
                        {
                              oldSol[j].a[i] = sol[j].a[i];
                              oldSol[j].v[i] = sol[j].v[i];
                              oldSol[j].r[i] = sol[j].r[i];
                        }
                  }
                  //euler(&sol[0], dt, Nobj);
            }
            else
            {
                  //leapFrog(&sol[0], &oldSol[0], dt, Nobj);
            }
            //euler(&sol[0], dt, Nobj);
            RK4(&sol[0], dt, Nobj);
      }


      //close all my files
      fclose(poss);
      fclose(energy1);
      fclose(angMom);
      
      return 0;
}
