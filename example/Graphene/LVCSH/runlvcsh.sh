#!/bin/bash
ensampledir='ensample'
for i in $(seq 6 6 36)
  do
  mkdir epwkq${i}
  cp ../epw/epw${i}.out epwkq${i}
  cp LVCSH.in epwkq${i}
  cp lvcsh.bsub epwkq${i}
  sed -i "2s:lvcsh:lvcsh-kq${i}-s0:g" epwkq${i}/lvcsh.bsub
  sed -i "s:epw.out:epw${i}.out:g" epwkq${i}/LVCSH.in
  
  cd epwkq${i}
    cat > runlvcsh.sh <<-EOF
		#!/bin/bash
		ensampledir='ensample'
		for j in \$(seq 0 1 0)
		  do
		  mkdir ${ensampledir}\${j}
		  cp LVCSH.in ${ensampledir}\${j}
		  cp lvcsh.bsub ${ensampledir}\${j}
		  cd ${ensampledir}\${j}
		  sed -i "2s:lvcsh-kq${i}-s0:lvcsh-kq${i}-s\${j}:g" lvcsh.bsub
		  bsub < lvcsh.bsub
		  cd ..
		  done
		EOF
		chmod +x runlvcsh.sh
		
    for j in $(seq 0 1 0)
      do
        mkdir ${ensampledir}${j}
        cp LVCSH.in ${ensampledir}${j}
        cp lvcsh.bsub ${ensampledir}${j}
        cd ${ensampledir}${j}
        sed -i "2s:lvcsh-kq${i}-s0:lvcsh-kq${i}-s${j}:g" lvcsh.bsub#
#        bsub < lvcsh.bsub
        cd ..
      done
  cd ..
  done
