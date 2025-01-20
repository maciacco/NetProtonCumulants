dir=$1

root -b -l <<EOF
.L Analysis.cxx
Analysis("${1}", "k2k2sk")
Analysis("${1}", "k6k2")
Analysis("${1}", "k4k2")
Analysis("${1}", "k1")
Analysis("${1}", "k2")
Analysis("${1}", "k3")
Analysis("${1}", "k4")
Analysis("${1}", "k5")
Analysis("${1}", "k6")
Analysis("${1}", "k2sk")
.q
EOF


