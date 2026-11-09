"""Complete list of negative discriminants D whose form class group has exponent <= 2, given the fundamental list.
Input: exp2_1e6.txt (output of forms.c for |D|<=10^6). Usage: python3 -I exp2.py exp2_1e6.txt"""
import sys
from sympy import factorint, jacobi_symbol
def kron(d,p):
    if p==2: return 0 if d%2==0 else (1 if d%8==1 else -1)
    return jacobi_symbol(d%p,p)
def fundamental(D):
    N=-D
    if D%4==1: return all(e==1 for e in factorint(N).values())
    if D%4==0:
        m=N//4
        return m%4 in (1,2) and all(e==1 for e in factorint(m).values())
    return False
def mu(D):
    N=-D; r=sum(1 for p in factorint(N) if p>2)
    if D%4==1: return r
    n=N//4
    if n%4==3: return r
    if n%4 in (1,2) or n%8==4: return r+1
    return r+2   # n = 0 mod 8
lst=[tuple(map(int,l.split())) for l in open(sys.argv[1])]
small={D:h for D,h in lst}
fund=sorted([D for D in small if fundamental(D)],reverse=True)
print('fundamental with exponent<=2:',len(fund),'largest |d| =',max(-d for d in fund))
found={}
# phi(f)/2^omega(f) is multiplicative with value p^(a-1)(p-1)/2 on p^a, >=1 except at p^a=2 (value 1/2).
# phi(f)<=6*2^omega(f) forces each such value <=12, so f divides F below.
from itertools import product
F={2:5,3:3,5:2,7:1,11:1,13:1,17:1,19:1,23:1}
for p,a in F.items():
    assert p**a//p*(p-1)/2<=12 and (p**(a+1)//p*(p-1)/2>12)
assert all(p*(p-1)//p/2>12 or p in F for p in [29,31,37])  # p>=29: (p-1)/2>=14
DIVS=[]
for es in product(*[range(a+1) for a in F.values()]):
    f=1
    for p,e in zip(F,es): f*=p**e
    DIVS.append(f)
print('divisors of F:',len(DIVS))
for d in fund:
    hd=small[d]; u0={-3:3,-4:2}.get(d,1)
    assert hd==2**(len(factorint(-d))-1)
    # f with phi(f)<= u*2^(omega(f)+1): enumerate f up to a safe bound
    for f in DIVS:
        fac=factorint(f); om=len(fac); u=1 if f==1 else u0
        val=f
        num=1;den=1
        from fractions import Fraction
        prod=Fraction(f)
        for p in fac: prod*=Fraction(p-kron(d,p),p)
        if prod > u*2**(om+1): continue
        D=f*f*d; hD=prod*hd/u
        assert hD.denominator==1
        if hD==2**(mu(D)-1): found[D]=int(hD)
# safety: phi(f) <= 6*2^omega(f) forces f < 5000
print('discriminants found:',len(found),'largest |D|=',max(-D for D in found))
print('agree with forms list:', found=={D:h for D,h in small.items()})
idoneal=sorted(-D//4 for D in found if D%4==0); print('idoneal:',len(idoneal),idoneal)
