iS="$1"

root -b -l <<EOF
.L ProcessTupleSingleParticle.cxx+
ProcessTupleSingleParticle($iS, 364, 365, false)
.q
EOF