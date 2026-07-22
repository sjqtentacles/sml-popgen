# sml-popgen

[![CI](https://github.com/sjqtentacles/sml-popgen/actions/workflows/ci.yml/badge.svg)](https://github.com/sjqtentacles/sml-popgen/actions/workflows/ci.yml)

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

## Example

`make example` builds and runs [`examples/demo.sml`](examples/demo.sml), which
computes allele frequencies, expected Hardy-Weinberg genotype counts, a
chi-square goodness-of-fit test, Fst, a five-generation Wright-Fisher drift
simulation from a fixed seed, and linkage disequilibrium for fixed inputs
(output is byte-identical under MLton and Poly/ML):

```
Population genetics worked example

Allele frequencies from genotype counts (AA=30, AB=40, BB=30):
  p = 0.5000   q = 0.5000

Expected Hardy-Weinberg genotype counts (n=100):
  E[AA] = 25.0000   E[AB] = 50.0000   E[BB] = 25.0000

Hardy-Weinberg chi-square goodness of fit (obs = 30,40,30):
  chiSq = 4.0000   df = 1

Fst (hT=0.5, hS=0.3):
  Fst = 0.4000

Wright-Fisher drift, p0=0.5, 2N=20, seed=12345, 5 generations:
  p' = 0.4500   seed' = 24224
  p' = 0.4000   seed' = 29426
  p' = 0.5000   seed' = 19571
  p' = 0.6000   seed' = 32745
  p' = 0.5000   seed' = 27958

Linkage disequilibrium (pAB=0.3, pA=0.5, pB=0.5):
  D = 0.0500   D' = 0.2000   r2 = 0.0400
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
