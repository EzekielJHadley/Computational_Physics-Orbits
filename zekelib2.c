#include <stdio.h>
#include <math.h>

struct planet
{
      double r[3];
      double v[3];
      double a[3];
      double mass;
};

//now calculate the distance between two planets
double dist(struct planet p1, struct planet p2)
{
      double length=0;
      int i;
      for (i=0; i<3; i++)
      {
            length += pow(p1.r[i] - p2.r[i], 2);
      }
      length = sqrt(length);
      return length;
}

double accel(struct planet *p, int Nobj)
{
      double d, dx[3];
      int i, j, k;
      //calculate the acceleration of each planet
      for(k=0; k<3; k++)
      {
            for(i=0; i<Nobj; i++)
            {
                  (*(p+i)).a[k] = 0;
            }
      }
      for(j=0; j<Nobj-1; j++)
      {
            for(i=j+1; i<=Nobj; i++)
            {
                  for(k=0; k<3; k++)
                  {
                        dx[k] = (*(p+i)).r[k] - (*(p+j)).r[k];
                  }
                  d = dist(*(p+i), *(p+j));
                    
                  for(k=0; k<3; k++)
                  {
                        (*(p+j)).a[k] += ((*(p+i)).mass)*dx[k]/pow(d, 3);
                        (*(p+i)).a[k] -= ((*(p+j)).mass)*dx[k]/pow(d, 3);
                        //printf("change? %f %f\n", sol[0].mass, sol[1].mass);
                  }
            }
      }
}

double accelEin(struct planet *p, int Nobj)
{
      double d, dx[3];
      int i, j, k;
      //calculate the schwartzchild radius of the sun
      double c2 = pow(63239.7263,2);
      double Rs = (((*(p)).mass)/c2 *2.0) * 1000.0;
      //calculate the acceleration of each planet
      for(k=0; k<3; k++)
      {
            for(i=0; i<Nobj; i++)
            {
                  (*(p+i)).a[k] = 0;
            }
      }
      for(j=0; j<Nobj-1; j++)
      {
            for(i=j+1; i<=Nobj; i++)
            {
                  for(k=0; k<3; k++)
                  {
                        dx[k] = (*(p+i)).r[k] - (*(p+j)).r[k];
                  }
                  d = dist(*(p+i), *(p+j));
                    
                  for(k=0; k<3; k++)
                  {
                        //find the accleration on the planet with the correction from relativity
                        (*(p+j)).a[k] += ((*(p+i)).mass)*dx[k]/pow(d, 3)*(1 + (3.0*Rs)/(2.0*d));
                        (*(p+i)).a[k] -= ((*(p+j)).mass)*dx[k]/pow(d, 3)*(1 + (3.0*Rs)/(2.0*d));
                        //printf("change? %f %f\n", sol[0].mass, sol[1].mass);
                  }
            }
      }
}

//do a single euler step
double euler(struct planet *p, double dt, int Nobj)
{
      int i, j;
      for(j=0; j<Nobj; j++)
      {
            for(i = 0; i<3; i++)
            {
                  //and new velocity at the same time
                  (*(p+j)).v[i] += ((*(p+j)).a[i])*dt;
                  //compute new possition
                  (*(p+j)).r[i] += ((*(p+j)).v[i])*dt;
            }
      }
}

//do a single Leap frog step
//an euler step needs to be done first
double leapFrog(struct planet *p, struct planet *old, double dt, int Nobj)
{
      int i, j;
      struct planet current[Nobj];
      //set aside the current values to be saved as the old
      for(j=0; j<Nobj; j++)
      {
            for(i=0; i<3; i++)
            {
                  current[j].a[i] = (*(p+j)).a[i];
                  current[j].v[i] = (*(p+j)).v[i];
                  current[j].r[i] = (*(p+j)).r[i];
            }
      }
      for(j=0; j<Nobj; j++)
      {
            for(i=0; i<3; i++)
            {
                  //now update possition
                  (*(p+j)).r[i] = ((*(old+j)).r[i]) + ((*(p+j)).v[i])*2.0*dt;
                  //now update velocity
                  (*(p+j)).v[i] = ((*(old+j)).v[i]) + ((*(p+j)).a[i])*2.0*dt;
            }
      }
      //now save the current to the old
      for(j=0; j<Nobj; j++)
      {
            for(i=0; i<3; i++)
            {
                  (*(old+j)).a[i] = current[j].a[i];
                  (*(old+j)).v[i] = current[j].v[i];
                  (*(old+j)).r[i] = current[j].r[i];
            }
      }
}

//do a single Runge-Kutta step aka RK4
double RK4(struct planet *p, double dt, int Nobj)
{
      int i, j;
      double v1[3][Nobj], v2[3][Nobj], v3[3][Nobj], v4[3][Nobj];
      double x1[3][Nobj], x2[3][Nobj], x3[3][Nobj], x4[3][Nobj];
      double oldR[3][Nobj], oldV[3][Nobj];
      for(j=0; j<Nobj; j++)
      {
            //save the current values for later use
            for(i = 0; i<3; i++)
            {
                  oldV[i][j] = (*(p+j)).v[i];
                  oldR[i][j] = (*(p+j)).r[i];
            }
            //get the values for K1
            for(i = 0; i<3; i++)
            {
                  v1[i][j] = (*(p+j)).a[i]*dt;
                  x1[i][j] = (*(p+j)).v[i]*dt;
                  (*(p+j)).v[i] += .5*v1[i][j];
                  (*(p+j)).r[i] += .5*x1[i][j];
            }
      }
      //update the acceleration
      accel(p, Nobj);
      //get values for K2
      for(j=0; j<Nobj; j++)
      {
            for(i = 0; i<3; i++)
            {
                  v2[i][j] = (*(p+j)).a[i]*dt;
                  x2[i][j] = (*(p+j)).v[i]*dt;
                  (*(p+j)).v[i] = oldV[i][j] + .5*v2[i][j];
                  (*(p+j)).r[i] = oldR[i][j] + .5*x2[i][j]; 
            }
      }
      //up date accleration
      accel(p, Nobj);
      //get values for K3
      for(j=0; j<Nobj; j++)
      {
            for(i = 0; i<3; i++)
            {
                  v3[i][j] = (*(p+j)).a[i]*dt;
                  x3[i][j] = (*(p+j)).v[i]*dt;
                  (*(p+j)).v[i] = oldV[i][j] + v3[i][j];
                  (*(p+j)).r[i] = oldR[i][j] + x3[i][j];
            }
      }
      //update acceleration
      accel(p, Nobj);
      //get values for K4
      for(j=0; j<Nobj; j++)
      {
            for(i = 0; i<3; i++)
            {
                  v4[i][j] = (*(p+j)).a[i]*dt;
                  x4[i][j] = (*(p+j)).v[i]*dt;
            }
      
            //now update the whole thing
            for(i=0; i<3; i++)
            {
                  (*(p+j)).v[i] = oldV[i][j] + (v1[i][j] + 2.0*v2[i][j] + 2.0*v3[i][j] + v4[i][j])/6.0;
                  (*(p+j)).r[i] = oldR[i][j] + (x1[i][j] + 2.0*x2[i][j] + 2.0*x3[i][j] + x4[i][j])/6.0;
            }
      }
}


//I need to find the cneter of mass of the system
//and turn it into the origin
double center(struct planet *p1, int Nobj)
{
      double cmX=0, cmY=0, cmZ=0;
      double TMass = 0;
      int i;
      //find the total mass
      for(i = 0; i<Nobj; i++)
      {
            TMass += (*(p1+i)).mass;
      }
      
      for(i = 0; i<Nobj; i++)
      {
            cmX += (*(p1+i)).mass * (*(p1+i)).r[0];
            cmY += (*(p1+i)).mass * (*(p1+i)).r[1];
            cmZ += (*(p1+i)).mass * (*(p1+i)).r[2];
      }

      cmX = cmX/TMass;
      cmY = cmY/TMass;
      cmZ = cmZ/TMass;

      for(i = 0; i<Nobj; i++)
      {
            (*(p1+i)).r[0] -= cmX;
            (*(p1+i)).r[1] -= cmY;
            (*(p1+i)).r[2] -= cmZ;
      }

}

//calculate the energies, both KE and U
//Then sum them up to total energy
double energy(struct planet p1[], double *E, int Nobj)
{
      double KE[Nobj], U[Nobj], V2=0, d;
      int i, j;
      //set the arrays to zero
      for(i=0; i<Nobj; i++)
      {
            KE[i] = 0;
            U[i] = 0;
      }

      //start calculating the potentials
      for(i=0; i<Nobj-1; i++)
      {
            for(j=i+1; j<Nobj; j++)
            {
                  d = dist(p1[i], p1[j]);
                  U[i] += -(p1[j].mass * p1[i].mass)/d;
                  U[j] += U[i];
            }
      }
      
      for(i=0; i<Nobj; i++)
      {      
            //calculate the KE
            V2 = 0;
            for(j = 0; j<3; j++)
            {
                  V2 += pow(p1[i].v[j], 2);
            }
            KE[i] = .5*p1[i].mass*V2;
            *(E+i) = KE[i] + U[i];
      }
}

//calculate the angular momentum
double momentum(struct planet p[], int Nobj, FILE *angMom)
{
      double L[Nobj][3];
      int i, j;
      for(i=0; i<Nobj; i++)
      {
            L[i][0] = p[i].mass*(p[i].r[1]*p[i].v[2] - p[i].r[2]*p[i].v[1]);
            L[i][1] = p[i].mass*(p[i].r[2]*p[i].v[0] - p[i].r[0]*p[i].v[2]);
            L[i][2] = p[i].mass*(p[i].r[0]*p[i].v[1] - p[i].r[1]*p[i].v[0]);
      }
      for(i=0; i<Nobj; i++)
      {
            fprintf(angMom, "%.8e \t", L[i][2]);
      }
      fprintf(angMom, "\n");
}

double initSolSystem(struct planet *p, int Nobj)
{
      int i, j;
      double unitMass = 1.98523584203e-29;
      double unitTime = 365.242198781;
      for(j=0; j<Nobj; j++)
      {
            for(i=0; i<3; i++)
            {
                  (*(p+j)).a[i] = 0;
                  (*(p+j)).v[i] = 0;
                  (*(p+j)).r[i] = 0;
                  (*(p+j)).mass = 0;
            }
      }
      //data from 1012 oct 03

      //for the sun
      (*(p+0)).mass = 1.9891e30*unitMass;
      (*(p+0)).v[0] = 6.33989263356e-6*unitTime;//5.905956508651875e-6*unitTime;//6.33989263356e-6*unitTime;
      (*(p+0)).v[1] = 0;//-2.301650373103100e-6*unitTime;
      (*(p+0)).v[2] = 0;//-1.277571249480375e-7*unitTime;
      (*(p+0)).r[0] = 0;//-1.829354013542788e-3;
      (*(p+0)).r[1] = -2.293188147356324e-3;
      (*(p+0)).r[2] = 0;//-3.036829152435623e-5;
      
      if(Nobj > 1)
      {
            //for mercury
            (*(p+1)).mass = 3.302e23*unitMass;
            (*(p+1)).v[0] = 1.62965737377e-2*unitTime;//1.943196031387812e-2*unitTime; //1.62965737377e-2*unitTime;
            (*(p+1)).v[1] = 0;//-1.136288032813267e-2*unitTime;
            (*(p+1)).v[2] = 0;//-2.710759466619381e-3*unitTime;
            (*(p+1)).r[0] = 0;//-2.1113003530905613e-1;
            (*(p+1)).r[1] = 4.65880520737e-1;//-4.150386252451500e-1; //4.65880520737e-1;
            (*(p+1)).r[2] = 0;//-1.455016657250516e-2;
      }
      else if(Nobj > 2)
      {
            //Venus
            (*(p+2)).mass = 48.685e23*unitMass;
            (*(p+2)).v[0] = -2.020996777812236e-2*unitTime;//2.03147262302e-2*unitTime;
            (*(p+2)).v[1] = 1.682250970681684e-3*unitTime;
            (*(p+2)).v[2] = 1.189678900463839e-3*unitTime;
            (*(p+2)).r[0] = 6.18478317557035e-2;
            (*(p+2)).r[1] = 7.148117232340789e-1;//7.17508477832e-1;
            (*(p+2)).r[2] = 6.120604188366068e-3;
      }
      else if(Nobj > 3)
      {
            //Earth
            (*(p+3)).mass = 5.9736e24*unitMass;
            (*(p+3)).v[0] = -3.260735654464636e-3*unitTime;//1.71829104585e-2*unitTime;
            (*(p+3)).v[1] = 1.687068507248210e-2*unitTime;
            (*(p+3)).v[2] = -3.824432953572980e-8*unitTime;
            (*(p+3)).r[0] = 9.834884170947437e-1;
            (*(p+3)).r[1] = 1.716247109979325e-1;//9.98350894941e-1;
            (*(p+3)).r[2] = -3.800852830374635e-5;
      }
      else if(Nobj > 4)
      {
            //Mars
            (*(p+4)).mass = 6.4185e23*unitMass;
            (*(p+4)).v[0] = 1.450847567676623e-2*unitTime;
            (*(p+4)).v[1] = 4.4222951313010323e-4*unitTime;
            (*(p+4)).v[2] = -3.468979065292308e-4*unitTime;
            (*(p+4)).r[0] = -8.082076946558051e-2;
            (*(p+4)).r[1] = -1.461715498863398e0;
            (*(p+4)).r[2] = -2.866806304677543e-2;
      }
      else if(Nobj > 5)
      {
            //Jupiter
            (*(p+5)).mass = 1898.13e24*unitMass;
            (*(p+5)).v[0] = -6.970166176636945e-3*unitTime;
            (*(p+5)).v[1] = 3.463008784094558e-3*unitTime;
            (*(p+5)).v[2] = 1.416166917469149e-4*unitTime;
            (*(p+5)).r[0] = 2.068672788690857e00;
            (*(p+5)).r[1] = 4.588206257139396e00;
            (*(p+5)).r[2] = -6.542759750160609e-2;
      }
      else if(Nobj > 6)
      {
            //Saturn
            (*(p+6)).mass = 5.68319e26*unitMass;
            (*(p+6)).v[0] = 2.608878065516174e-3*unitTime;
            (*(p+6)).v[1] = -4.771272541042325e-3*unitTime;
            (*(p+6)).v[2] = -2.116235786925556e-5*unitTime;
            (*(p+6)).r[0] = -8.327612684710576e0;
            (*(p+6)).r[1] = -5.094385183775453e0;
            (*(p+6)).r[2] = 4.199942267734638e-1;
      }
      else if(Nobj > 7)
      {
            //Uranus
            (*(p+7)).mass = 86.8103e24*unitMass;
            (*(p+7)).v[0] = -4.670902041438952e-4*unitTime;
            (*(p+7)).v[1] = 3.725219577753643e-3*unitTime;
            (*(p+7)).v[2] = 1.991137473477052e-5*unitTime;
            (*(p+7)).r[0] = 1.993405528582205e1;
            (*(p+7)).r[1] = 2.235992140925989e0;
            (*(p+7)).r[2] = -2.499490813599057e-1;
      }
      else if(Nobj > 8)
      {
            //Neptune
            (*(p+8)).mass = 102.41e24*unitMass;
            (*(p+8)).v[0] = 1.463438294703336e-3*unitTime;
            (*(p+8)).v[1] = 2.784855619935944e-3*unitTime;
            (*(p+8)).v[2] = -9.139710872903500e-5*unitTime;
            (*(p+8)).r[0] = 2.642907673505876e1;
            (*(p+8)).r[1] = -1.417691162096594e1;
            (*(p+8)).r[2] = -3.171300241886222e-1;
      }
      else if(Nobj > 9)
      {
            //Pluto
            (*(p+9)).mass = 1.314e22*unitMass;
            (*(p+9)).v[0] = 3.154161855758214e-3*unitTime;
            (*(p+9)).v[1] = -1.661644887600797e-4*unitTime;
            (*(p+9)).v[2] = -8.992173908334235e-4*unitTime;
            (*(p+9)).r[0] = 4.823262525837811e0;
            (*(p+9)).r[1] = -3.188501285727357e1;
            (*(p+9)).r[2] = 2.016697629154323e0;
      }
}

double initBirthday(struct planet *p, int Nobj)
{
      int i, j;
      double unitMass = 1.98523584203e-29;
      double unitTime = 365.242198781;
      for(j=0; j<Nobj; j++)
      {
            for(i=0; i<3; i++)
            {
                  (*(p+j)).a[i] = 0;
                  (*(p+j)).v[i] = 0;
                  (*(p+j)).r[i] = 0;
                  (*(p+j)).mass = 0;
            }
      }
      //data from 1012 oct 03
      //for the sun
      (*(p+0)).mass = 1.9891e30*unitMass;
      (*(p+0)).v[0] = 4.011510646586643e-8*unitTime;//6.33989263356e-6*unitTime;
      (*(p+0)).v[1] = -7.022384395418655e-6*unitTime;
      (*(p+0)).v[2] = 4.578644757741857e-8*unitTime;
      (*(p+0)).r[0] = -4.117181893959415e-3;
      (*(p+0)).r[1] = 4.030104134802236e-3;
      (*(p+0)).r[2] = 1.020933277550742e-5;

      //for mercury
      (*(p+1)).mass = 3.302e23*unitMass;
      (*(p+1)).v[0] = 2.1328554738055118e-2*unitTime; //1.62965737377e-2*unitTime;
      (*(p+1)).v[1] = -6.584182302397356e-3*unitTime;
      (*(p+1)).v[2] = -2.495465766678154e-3*unitTime;
      (*(p+1)).r[0] = -1.363577554080288e-1;
      (*(p+1)).r[1] = -4.426029900363971e-1; //4.65880520737e-1;
      (*(p+1)).r[2] = -2.432601184705706e-2;

      //Venus
      (*(p+2)).mass = 48.685e23*unitMass;
      (*(p+2)).v[0] = -1.606860544558468e-2*unitTime;//2.03147262302e-2*unitTime;
      (*(p+2)).v[1] = 1.227147869048631e-2*unitTime;
      (*(p+2)).v[2] = 1.095129024912683e-3*unitTime;
      (*(p+2)).r[0] = 4.373922337744308e-1;
      (*(p+2)).r[1] = 5.7755228177878194e-1;//7.17508477832e-1;
      (*(p+2)).r[2] = -1.768781991530748e-2;

      //Earth
      (*(p+3)).mass = 5.9736e24*unitMass;
      (*(p+3)).v[0] = 1.691789197514915e-2*unitTime;//1.71829104585e-2*unitTime;
      (*(p+3)).v[1] = 5.743007779890157e-4*unitTime;
      (*(p+3)).v[2] = -2.097885393462414e-7*unitTime;
      (*(p+3)).r[0] = 3.4422105598773606e-2;
      (*(p+3)).r[1] = -1.011705332533769e0;//9.98350894941e-1;
      (*(p+3)).r[2] = -1.954009115263795e-5;

      //Mars
      (*(p+4)).mass = 6.4185e23*unitMass;
      (*(p+4)).v[0] = -1.103647285307430e-2*unitTime;
      (*(p+4)).v[1] = -6.690469371173629e-3*unitTime;
      (*(p+4)).v[2] = 1.316311517607374e-4*unitTime;
      (*(p+4)).r[0] = -9.267313961087821e-1;
      (*(p+4)).r[1] = 1.359841114139311e0;
      (*(p+4)).r[2] = 5.110978358055619e-2;

      //Jupiter
      (*(p+5)).mass = 1898.13e24*unitMass;
      (*(p+5)).v[0] = -1.945102751124479e-3*unitTime;
      (*(p+5)).v[1] = 7.667536090825971e-3*unitTime;
      (*(p+5)).v[2] = 1.190361013582632e-5*unitTime;
      (*(p+5)).r[0] = 4.796572892023593e00;
      (*(p+5)).r[1] = 1.215618777931425e00;
      (*(p+5)).r[2] = -1.125204014582976e-1;
      
      //Saturn
      (*(p+6)).mass = 5.68319e26*unitMass;
      (*(p+6)).v[0] = 5.165553653401751e-3*unitTime;
      (*(p+6)).v[1] = -1.131997900150658e-3*unitTime;
      (*(p+6)).v[2] = -1.859530096219814e-4*unitTime;
      (*(p+6)).r[0] = -2.004697654019709e00;
      (*(p+6)).r[1] = -9.820464838975864e00;
      (*(p+6)).r[2] = 2.509024139848868e-1;
      
      //Uranus
      (*(p+7)).mass = 86.8103e24*unitMass;
      (*(p+7)).v[0] = 3.888494147837374e-3*unitTime;
      (*(p+7)).v[1] = -5.272462117463120e-4*unitTime;
      (*(p+7)).v[2] = -5.238588203435185e-5*unitTime;
      (*(p+7)).r[0] = -1.680275759605606e00;
      (*(p+7)).r[1] = -1.913757623852844e1;
      (*(p+7)).r[2] = -4.927765004886140e-2;
      
      //Neptune
      (*(p+8)).mass = 102.41e24*unitMass;
      (*(p+8)).v[0] = 3.097617494068595e-3*unitTime;
      (*(p+8)).v[1] = 3.868551921065853e-4*unitTime;
      (*(p+8)).v[2] = -7.946112808072177e-5*unitTime;
      (*(p+8)).r[0] = 3.557174521315203e00;
      (*(p+8)).r[1] = -3.001533022382608e1;
      (*(p+8)).r[2] = 5.360910435892920e-1;

      //Pluto
      (*(p+9)).mass = 1.314e22*unitMass;
      (*(p+9)).v[0] = 2.204445637382696e-3*unitTime;
      (*(p+9)).v[1] = -2.730599258465571e-3*unitTime;
      (*(p+9)).v[2] = -3.378824667413577e-4*unitTime;
      (*(p+9)).r[0] = -2.210222035328136e1;
      (*(p+9)).r[1] = -1.798753682760293e1;
      (*(p+9)).r[2] = 8.317965397401359e0;
}
