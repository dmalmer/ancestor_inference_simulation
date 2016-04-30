
(* run with: math -run "<<run_rabbit.m" *)

packageDir = "/Users/daniel/Downloads/RABBIT/RABBIT_Packages/"
wkDir = "/Users/daniel/Dropbox/AncestorInference/simulation/rabbit/"

SetDirectory[packageDir]
Needs["MagicReconstruct`"]
SetDirectory[wkDir]

(* parameters: *)
(*  40 gens of random mating, followed by 30 gens of sibling inbreeding *)
popScheme = Join[Table["RM1-E-196", {40}], Table["Sibling", {30}]];
epsF = eps = 0.005;
model = "jointModel";
estfun = "origViterbiDecoding";

inputFile = "./rabbit_markers.csv"
mouseSNPs = Import[inputFile]

outputFile = "rabbit_output_chr1.txt";

Print["starting magicReconstruct"]

magicReconstruct[inputFile, model, epsF, eps, popScheme, outputFile, HMMMethod -> estfun, PrintTimeElapsed -> True];

Print["magicReconstruct finished"]

summaryFile = StringDrop[outputFile, -4] <> "_summary.csv";
saveAsSummaryMR[outputFile, summaryFile]

res = getSummaryMR[summaryFile];
Print[First[res] // TableForm]

Exit[]
