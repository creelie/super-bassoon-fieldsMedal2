/* For every discriminant D<0 (D = 0,1 mod 4) with |D| <= X, enumerate primitive reduced forms (a,b,c),
   |b|<=a<=c, b>=0 if |b|==a or a==c, gcd(a,b,c)=1.  Print D, h(D), and whether every reduced form is ambiguous
   (b==0 || b==a || a==c), which holds iff the form class group has exponent <= 2.  Prints only D with exponent <= 2. */
#include <stdio.h>
#include <stdlib.h>
static long g(long a,long b){ if(a<0)a=-a; if(b<0)b=-b; while(b){long t=a%b;a=b;b=t;} return a; }
int main(int argc,char**argv){
  long X=atol(argv[1]);
  for(long N=3; N<=X; N++){ if(N%4==1||N%4==2) continue; long D=-N;
    long h=0; int amb=1;
    for(long a=1; 3*a*a<=N; a++){
      for(long b=-a+1; b<=a; b++){
        if(((b%2)+2)%2 != N%2) continue;
        long num=b*b+N; if(num%(4*a)) continue; long c=num/(4*a); if(c<a) continue;
        if(c==a && b<0) continue;
        if(g(g(a,b),c)!=1) continue;
        h++; if(!(b==0||b==a||a==c)) amb=0;
      }
    }
    if(amb) printf("%ld %ld\n",D,h);
  }
  return 0;
}
