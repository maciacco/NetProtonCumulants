min=$1
max=$2

root -b -l <<EOF
///.L La2PrPt.cxx
///La2PrPt()

//.L ReadTree.cxx
//ReadTree("newTree_LHC22f3_20250221", "LHC22f3_01", ${min}, ${max}, true)

///.L ReadTree.cxx
///ReadTree("newTree_LHC22f3_20250221", "LHC22f3_wd", ${min}, ${max}, true, 1)

//.L PrEffwFD.cxx
//PrEffwFD("LHC22f3_010_var", "", "prEff", ${min}, ${max})


//.L ReadTree.cxx
//ReadTree("newTree_LHC18pp_20250611", "LHC18ppTrig_HM", ${min}, ${max}, false)


.L ReadTreeMix.cxx
ReadTreeMix("newTree_LHC18pp_20250611", "LHC18ppTrig_HM", ${min}, ${max}, false)

//.L ProcessTuple.cxx+
//.L ProcessTupleSingleParticleHighOrder.cxx+
.q
EOF

#cat cmdz | xargs -P 20 -I CMD bash -c CMD
