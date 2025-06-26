iS="$1"

root -b -l <<EOF
.L ProcessTupleSingleParticleHighOrder.cxx+
ProcessTupleSingleParticleHighOrder($iS, 81, 82, false)
//.L ProcessTuple.cxx+
//ProcessTuple($iS, 81, 82, false)
.q
EOF
