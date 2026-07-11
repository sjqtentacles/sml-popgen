(* demo.sml - population genetics calculations: allele frequencies, Hardy-
   Weinberg equilibrium, Fst, Wright-Fisher drift, and linkage disequilibrium.
   Deterministic: fixed literal inputs and a fixed LCG seed for the drift
   simulation, so output is identical on every run and both compilers. *)

structure G = Popgen

fun fmt n x =
  let val x = if Real.== (x, 0.0) then 0.0 else x
  in Real.fmt (StringCvt.FIX (SOME n)) x
  end

val () = print "Population genetics worked example\n"

val () = print "\nAllele frequencies from genotype counts (AA=30, AB=40, BB=30):\n"
val {p, q} = G.alleleFreq {aa=30, ab=40, bb=30}
val () = print ("  p = " ^ fmt 4 p ^ "   q = " ^ fmt 4 q ^ "\n")

val () = print "\nExpected Hardy-Weinberg genotype counts (n=100):\n"
val {aa=eAA, ab=eAB, bb=eBB} = G.genotypeExpected {p=p, n=100}
val () = print ("  E[AA] = " ^ fmt 4 eAA ^ "   E[AB] = " ^ fmt 4 eAB
                ^ "   E[BB] = " ^ fmt 4 eBB ^ "\n")

val () = print "\nHardy-Weinberg chi-square goodness of fit (obs = 30,40,30):\n"
val {chiSq, df} = G.hardyWeinbergChiSq {obs=(30,40,30)}
val () = print ("  chiSq = " ^ fmt 4 chiSq ^ "   df = " ^ Int.toString df ^ "\n")

val () = print "\nFst (hT=0.5, hS=0.3):\n"
val () = print ("  Fst = " ^ fmt 4 (G.fst {hT=0.5, hS=0.3}) ^ "\n")

val () = print "\nWright-Fisher drift, p0=0.5, 2N=20, seed=12345, 5 generations:\n"
fun drift 0 (_:real) (_:int) = ([] : (real * int) list)
  | drift k p seed =
      let val {p', seed'} = G.wrightFisher {p=p, pop2N=20, seed=seed}
      in (p', seed') :: drift (k-1) p' seed' end
val trajectory = drift 5 0.5 12345
val () =
  List.app
    (fn (p', seed') => print ("  p' = " ^ fmt 4 p' ^ "   seed' = " ^ Int.toString seed' ^ "\n"))
    trajectory

val () = print "\nLinkage disequilibrium (pAB=0.3, pA=0.5, pB=0.5):\n"
val {d, dPrime, r2} = G.ld {pAB=0.3, pA=0.5, pB=0.5}
val () = print ("  D = " ^ fmt 4 d ^ "   D' = " ^ fmt 4 dPrime ^ "   r2 = " ^ fmt 4 r2 ^ "\n")
