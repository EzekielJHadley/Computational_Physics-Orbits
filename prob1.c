#include <stdio.h>
#include <math.h>
#include "zekelib2.c"

int main(void)
{
      int Nobj = 2, n = 1024;
      int i, j, k;
      struct planet sol[Nobj];
      struct planet oldSol[Nobj];
      double T = 30, dt = 1/(double)(n);
      double dx[3], d, t, E[Nobj], E0[Nobj];
      FILE *file, *energy1, *angMom;
      file = fopen("possition", "wt");
      energy1 = fopen("energy", "wt");
      angMom = fopen("momentum", "wt");
      printf("dt = %f\n", dt);
      //initiate my two objects
      for(j=0; j<Nobj; j++)
      {
            for(i=0; i<3; i++)
            {
                  sol[j].r[i]=0;
                  sol[j].v[i]=0;
                  sol[j].a[i]=0;
            }
      }
      sol[0].mass = 10; //39.4883261338;
      sol[0].r[0] = -1;
      sol[0].v[1] = -1;
      sol[1].mass = 10; //pow(1.1859004826, -4);
      sol[1].r[0] = 1;
      sol[1].v[1] = 1; //6.32834939021;
      
      //sol[0].v[1] = - sol[1].v[1]*sol[1].mass/sol[0].mass;

      //this will convert the coordinates to a
      //center of mass coordinate system
      center(&sol[0], Nobj);

      //get the initial energy and angular momentum
      energy(sol, &E0[0], Nobj);      
      momentum(sol, Nobj, angMom); 

      //now start solving for the acceleration
      for(t=0; t<=T; t+=dt)
      {
            //save to a file every 16 time steps
            if((int)(t/dt)%32 == 0)
            {
                  //printf("do i get here? %.8e \n", sol[1].r[1]);
                  // print the energies
                  energy(sol, &E[0], Nobj);
                  momentum(sol, Nobj, angMom);
                  fprintf(energy1, "%.8e   %.8e       %.8e   %.8e\n", E[0], E[0]-E0[0], E[1], E[1]-E0[1]);
                  fprintf(file, "%.8e   %.8e         %.8e   %.8e\n", sol[0].r[0], sol[0].r[1], sol[1].r[0], sol[1].r[1]);
            }
            
            //compute the acceleration
            accel(&sol[0], Nobj);

            //now do time steps for every planet
            if(t==0)
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
            euler(&sol[0], dt, Nobj);
            //RK4(&sol[0], dt, Nobj);
      }


      //close all my files
      fclose(file);
      fclose(energy1);
      fclose(angMom);
      
      return 0;
}
