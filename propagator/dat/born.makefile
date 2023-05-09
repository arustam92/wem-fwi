
dslow.H:
	Spike n1=30 n2=200 n3=1 o1=0 o2=0 o3=0 d1=10 d2=10 d3=1 k1=15 k2=100 mag=0.00001 | \
	Pad beg2=50 end2=50 extend=1 | Transp > $@


born.H: slow.H wave.H ${B}/checkBorn.x dslow.H
	${B}/checkBorn.x < slow.H dvel=dslow.H wave=wave.H par=geometry.par wfld=born_wfld.H > $@
adjBorn.H: bg.H ${B}/check_adjBorn.x wave.H
	${B}/check_adjBorn.x < slow.H data=ph1.H wave=wave.H par=geometry.par > $@

dotBorn: slow.H wave.H ${B}/dotBorn.x dslow.H
	${B}/dotBorn.x < slow.H dvel=dslow.H wave=wave.H par=geometry.par > /dev/null
