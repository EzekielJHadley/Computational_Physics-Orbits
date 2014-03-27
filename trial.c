#include <stdio.h>
#include <math.h>

struct stuff
{
      double x, v, a;
};

double chang(struct stuff p[], struct stuff *new, int j)
{
      (*new).x = j--;
      (*new).v = j--;
      (*new).a = j;
      printf("change? %f   %f    %f\n", p[0].x, p[0].v, p[0].a);
}

int main(void)
{
      struct stuff man[2];
      man[0].x=0;
      man[1].x=0;
      man[0].v=0;
      man[1].v=0;
      man[0].a=0;
      man[1].a=0;
      chang(man, &man[0], 5);
      printf("here are the values: x=%f, v=%f, a=%f\n", man[0].x, man[0].v, man[0].a);
      return 0;
}
