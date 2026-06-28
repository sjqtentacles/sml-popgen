# sml-popgen

Zero-dependency Standard ML library for population genetics calculations.

## API

```sml
signature POPGEN =
sig
  val alleleFreq         : {aa:int, ab:int, bb:int} -> {p:real, q:real}
  val genotypeExpected   : {p:real, n:int} -> {aa:real, ab:real, bb:real}
  val hardyWeinbergChiSq : {obs: int * int * int} -> {chiSq:real, df:int}

  val fst : {hT:real, hS:real} -> real

  val wrightFisher : {p:real, pop2N:int, seed:int} -> {p':real, seed':int}

  val ld : {pAB:real, pA:real, pB:real} -> {d:real, dPrime:real, r2:real}
end
```

## Worked example

```sml
(* Allele frequencies from observed genotype counts *)
val {p, q} = Popgen.alleleFreq {aa=25, ab=50, bb=25}
(* p = 0.5, q = 0.5 *)

(* Hardy-Weinberg chi-square test *)
val {chiSq, df} = Popgen.hardyWeinbergChiSq {obs=(25, 50, 25)}
(* chiSq ≈ 0.0, df = 1 → fits HW equilibrium *)

(* Wright-Fisher drift simulation *)
val {p'=newP, seed'=s} = Popgen.wrightFisher {p=0.5, pop2N=200, seed=12345}
```

## Scope and limitations

- Covers allele/genotype frequencies, Hardy-Weinberg test, Fst, Wright-Fisher drift, and LD (D, D', r²).
- Wright-Fisher uses a simple linear congruential generator; results are deterministic but not cryptographically random.
- Does not model selection, mutation, migration, or multiple alleles.
- `hardyWeinbergChiSq` always reports `df=1` (appropriate for one locus, two alleles).

## Build and test

Requires [MLton](http://mlton.org/) and Poly/ML in PATH.

```
make all-tests
```
