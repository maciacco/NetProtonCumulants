iS="$1"

root -b -l <<EOF
.L ProcessTupleSingleParticleHighOrder.cxx+
ProcessTupleSingleParticleHighOrder($iS, 364, 365, false)
//.L ProcessTuple.cxx+
//ProcessTuple($iS, 364, 365, false)
.q
EOF
