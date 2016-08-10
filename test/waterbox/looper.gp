set ylabel "Average potential energy"
set ytics nomirror
set y2label "Weighted histogram"
set y2range [0:]
set y2tics nomirror
plot "<head -n61 wbst.rst | tail -n +2" u (1/$1):(($4>0)?$2:1/0) w lp t "Pot. Energy", \
     "" u (1/$1):($4*$8) w lp axes x1y2 t "Flat"
reset
pause 10
reread
