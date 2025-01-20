dir=$1

root -b -l <<EOF
.L AnalysisSingleParticleHighOrder.cxx
AnalysisSingleParticleHighOrder("${1}", "k2k1")
AnalysisSingleParticleHighOrder("${1}", "k1")
AnalysisSingleParticleHighOrder("${1}", "k2")
AnalysisSingleParticleHighOrder("${1}", "k3")
AnalysisSingleParticleHighOrder("${1}", "k4")
AnalysisSingleParticleHighOrder("${1}", "k5")
AnalysisSingleParticleHighOrder("${1}", "k6")
.q
EOF


