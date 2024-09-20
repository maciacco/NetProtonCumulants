iS="$1"

root -b -l <<EOF
.x ProcessTuple.cxx($iS, 0, 3)
.q
EOF