iS="$1"

root -b -l <<EOF
//.L ProcessTupleSingleParticleHighOrder.cxx+
//ProcessTupleSingleParticleHighOrder($iS, 0, 750, false)
.L ProcessTuple.cxx+
ProcessTuple($iS, 0, 750, false)
.q
EOF
