structure Popgen :> POPGEN =
struct
  fun alleleFreq {aa, ab, bb} =
    let
      val total = 2 * (aa + ab + bb)
      val pCount = 2*aa + ab
      val p = Real.fromInt pCount / Real.fromInt total
    in {p=p, q=1.0-p} end

  fun genotypeExpected {p, n} =
    let
      val q = 1.0 - p
      val nr = Real.fromInt n
    in {aa = p*p*nr, ab = 2.0*p*q*nr, bb = q*q*nr} end

  fun hardyWeinbergChiSq {obs=(aa, ab, bb)} =
    let
      val n = aa + ab + bb
      val {p, q=_} = alleleFreq {aa=aa, ab=ab, bb=bb}
      val {aa=eAA, ab=eAB, bb=eBB} = genotypeExpected {p=p, n=n}
      fun chi2term obs e =
        if e < 1e~10 then 0.0
        else let val d = Real.fromInt obs - e in d*d/e end
    in {chiSq = chi2term aa eAA + chi2term ab eAB + chi2term bb eBB, df=1} end

  fun fst {hT, hS} =
    if Real.abs hT < 1e~10 then 0.0
    else (hT - hS) / hT

  (* LCG pseudo-random generator with small modulus to avoid overflow *)
  val lcgM = 32749  (* prime near 2^15; max product 3141*32748 < 2^27, fits any int *)
  fun lcgNext seed = (seed * 3141 + 2718) mod lcgM

  fun wrightFisher {p, pop2N, seed} =
    let
      fun loop 0 s acc = (acc, s)
        | loop n s acc =
            let
              val s' = lcgNext s
              val u = Real.fromInt s' / Real.fromInt lcgM
              val hit = if u < p then 1 else 0
            in loop (n-1) s' (acc + hit) end
      val (successes, seed') = loop pop2N seed 0
      val p' = Real.fromInt successes / Real.fromInt pop2N
    in {p'=p', seed'=seed'} end

  fun ld {pAB, pA, pB} =
    let
      val d = pAB - pA * pB
      val dMax =
        if d >= 0.0
        then Real.min (pA*(1.0-pB), (1.0-pA)*pB)
        else Real.min (pA*pB, (1.0-pA)*(1.0-pB))
      val dPrime =
        if Real.abs dMax < 1e~10 then 0.0
        else d / dMax
      val denom = pA * (1.0 - pA) * pB * (1.0 - pB)
      val r2 =
        if denom < 1e~10 then 0.0
        else d*d / denom
    in {d=d, dPrime=dPrime, r2=r2} end
end
