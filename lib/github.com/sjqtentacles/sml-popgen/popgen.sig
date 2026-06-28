signature POPGEN =
sig
  (* Allele frequencies from genotype counts: AA, AB, BB counts *)
  val alleleFreq : {aa:int, ab:int, bb:int} -> {p:real, q:real}

  (* Expected HW genotype counts for given p, total count n *)
  val genotypeExpected : {p:real, n:int} -> {aa:real, ab:real, bb:real}

  (* Hardy-Weinberg chi-square goodness of fit (df=1).
     obs: observed AA, AB, BB counts. *)
  val hardyWeinbergChiSq : {obs: int * int * int} -> {chiSq:real, df:int}

  (* Fst = (hT - hS) / hT *)
  val fst : {hT:real, hS:real} -> real

  (* Wright-Fisher drift step: p=allele freq, pop2N=2N (diploid individuals),
     seed=LCG seed. Returns (p', seed'). *)
  val wrightFisher : {p:real, pop2N:int, seed:int} -> {p':real, seed':int}

  (* Linkage disequilibrium: pAB=freq(AB haplotype), pA=freq(A allele), pB=freq(B allele) *)
  val ld : {pAB:real, pA:real, pB:real} -> {d:real, dPrime:real, r2:real}
end
