structure Tests =
struct
  open Harness
  fun close name (e, a, eps) = check name (Real.abs (e - a) <= eps)

  fun run () =
  let
    val () = section "allele frequency"
    val {p, q} = Popgen.alleleFreq {aa=25, ab=50, bb=25}
    val () = close "p=0.5" (0.5, p, 1e~12)
    val () = close "q=0.5" (0.5, q, 1e~12)
    val () = close "p+q=1" (1.0, p+q, 1e~12)

    val () = section "Hardy-Weinberg"
    val e = Popgen.genotypeExpected {p=0.5, n=100}
    val () = close "expected AA" (25.0, #aa e, 1e~12)
    val () = close "expected AB" (50.0, #ab e, 1e~12)
    val () = close "expected BB" (25.0, #bb e, 1e~12)
    val {chiSq, df} = Popgen.hardyWeinbergChiSq {obs=(25,50,25)}
    val () = close "chi-sq=0 for HW" (0.0, chiSq, 1e~9)
    val () = checkInt "df=1" (1, df)
    val p2 = 0.5*0.5 + 2.0*0.5*0.5 + 0.5*0.5
    val () = close "HW identity" (1.0, p2, 1e~12)

    val () = section "Fst"
    val () = close "Fst hT=hS" (0.0, Popgen.fst {hT=0.4, hS=0.4}, 1e~12)
    val () = close "Fst" (0.25, Popgen.fst {hT=0.4, hS=0.3}, 1e~9)

    val () = section "Wright-Fisher"
    val {p'=p1, seed'=s1} = Popgen.wrightFisher {p=0.5, pop2N=100, seed=42}
    val {p'=p2, seed'=_} = Popgen.wrightFisher {p=0.5, pop2N=100, seed=42}
    val () = close "WF reproducible" (p1, p2, 1e~12)
    val () = check "WF p' in [0,1]" (p1 >= 0.0 andalso p1 <= 1.0)
    val _ = s1

    val () = section "Linkage disequilibrium"
    val {d=d0, dPrime=dp0, r2=r20} = Popgen.ld {pAB=0.25, pA=0.5, pB=0.5}
    val () = close "D=0 at equilibrium" (0.0, d0, 1e~12)
    val () = close "D'=0 at equilibrium" (0.0, dp0, 1e~12)
    val () = close "r2=0 at equilibrium" (0.0, r20, 1e~12)
    val {d=d1, dPrime=dp1, r2=r21} = Popgen.ld {pAB=0.5, pA=0.5, pB=0.5}
    val () = close "D max" (0.25, d1, 1e~12)
    val () = close "D' max LD" (1.0, dp1, 1e~9)
    val () = close "r2 max LD" (1.0, r21, 1e~9)

  in Harness.run () end
end
