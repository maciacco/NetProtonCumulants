root -b -l <<EOF
.L ReadTree.cxx+
//ReadTree("newTree_mc_TOF_extendedPt_06122024", "LHC21pp", 364, 365, true)
//.L PrEff.cxx
//PrEff("LHC21pp0_var", "prEff", 364, 365)
ReadTree("newTree_TOF_extendedPt_04122024", "LHC18pp", 364, 365, false)
//.L ProcessTuple.cxx+
.q
EOF

#cat cmdz | xargs -P 10 -I CMD bash -c CMD
