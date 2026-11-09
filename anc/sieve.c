/* Fundamental discriminants d<0 with Q0<|d|<=2^40 such that no prime p with p^2<|d|/4 splits in Q(sqrt d).
   Every such field with class group of exponent <=2 must appear in the output (see the paper, Lemma on small split primes).
   case 0: d=-q, q = 3 mod 8 (d = 5 mod 8; d = 1 mod 8 would make 2 split).  Condition: (-q/p) != 1 for odd p.
   case 1: d=-4m, m = 1,2 mod 4.  Condition: (-m/p) != 1 for odd p.
   Output: candidate n (q or m) passing all odd primes p with p^2 < |d|/4 and p < 1000.  gcc -O2 -fopenmp */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static int primes[200], np=0;
static unsigned char bad[1000][1000]; /* bad[p][r]=1 iff r!=0 and -r is a nonzero square mod p */
int main(int argc,char**argv){
  int cas=atoi(argv[1]);
  const uint64_t QMAX=1ULL<<40, Q0=(argc>2)?strtoull(argv[2],0,10):1000000ULL;
  for(int p=3;p<1000;p++){int ip=1;for(int k=2;k*k<=p;k++)if(p%k==0)ip=0; if(ip) primes[np++]=p;}
  for(int i=0;i<np;i++){int p=primes[i]; unsigned char sq[1000]={0}; for(int x=1;x<p;x++) sq[(x*x)%p]=1;
     for(int r=1;r<p;r++) bad[p][r]=sq[(p-r)%p];}
  /* wheel over primes 3..23 plus 2-part */
  uint64_t M=(cas==0)?8:4; uint64_t *res=malloc(sizeof(uint64_t)*4); size_t nr=0;
  if(cas==0){res[nr++]=3;} else {res[nr++]=1; res[nr++]=2;}
  int nw=0; for(int i=0;i<np && primes[i]<=23;i++){ int p=primes[i]; nw=i+1;
     size_t cnt=0; uint64_t *nres=malloc(sizeof(uint64_t)*nr*p);
     for(size_t j=0;j<nr;j++) for(int t=0;t<p;t++){ uint64_t x=res[j]+t*M; if(!bad[p][x%p]) nres[cnt++]=x; }
     free(res); res=nres; nr=cnt; M*=p; }
  fprintf(stderr,"case %d: modulus %llu, %zu residues\n",cas,(unsigned long long)M,nr);
  uint64_t nmax = (cas==0)?QMAX:QMAX/4, nmin=(cas==0)?Q0:Q0/4;
  long long total=0, out=0;
  #pragma omp parallel for schedule(dynamic,64) reduction(+:total,out)
  for(size_t j=0;j<nr;j++){
    for(uint64_t n=res[j]; n<=nmax; n+=M){
      if(n<=nmin) continue;
      total++;
      uint64_t q=(cas==0)?n:4*n; int okc=1;
      for(int i=nw;i<np;i++){ uint64_t p=primes[i]; if(4*p*p>=q) break; if(bad[p][n%p]){okc=0;break;} }
      /* the wheel primes must also be within range: they are automatically fine (admissible for all p<=23) */
      if(okc){
        #pragma omp critical
        { printf("%d %llu\n",cas,(unsigned long long)n); }
        out++;
      }
    }
  }
  fprintf(stderr,"case %d: tested %lld, survivors %lld\n",cas,total,out);
  return 0;
}
