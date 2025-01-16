root -b -l <<EOF
//.L ReadTree.cxx+
//ReadTree("newTree_LHC22f3_WD", "LHC22f3", 364, 365, true)
//.L ReadTree.cxx+
//ReadTree("newTree_LHC22f3_WD", "LHC22f3_wd", 364, 365, true, 1)
//.L PrEffwFD.cxx
//PrEffwFD("LHC22f30_var", "LHC22f3_wd0_var", "prEff", 364, 365)
//PrEff("LHC22f30_var", "prEff", 364, 365)
//PrEff("LHC22f3_wd0_var", "prEffWD", 364, 365)
//.L ReadTree.cxx+
//ReadTree("newTree_LHC18v2_20241230", "LHC18ppTrig_HM", 364, 365, false)
//.L ProcessTuple.cxx+
.L ProcessTupleSingleParticleHighOrder.cxx+
.q
EOF

cat cmdz | xargs -P 20 -I CMD bash -c CMD
