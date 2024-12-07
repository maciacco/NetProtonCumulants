iS="$1"

root -b -l <<EOF
.L ProcessTuple.cxx+
ProcessTuple($iS, 364, 365, false)
.q
EOF