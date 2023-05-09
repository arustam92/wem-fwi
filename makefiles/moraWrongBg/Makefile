Wav/moraWav.H:
	Spike n1=750 d1=0.008 o1=0 k1=100 | \
	Bandpass flo=4 fhi=12 nplo=20 nphi=20 > $@

##########################################################################│
Vel/moraSlowTrue.H:
	Gauss n1=170 n2=500 n3=1 d1=10 d2=10 d3=1 o1=0.0 o2=0.0 o3=0.0 velback=0.0 maxvel=-150 var=1e6 max1=800 max2=2500 zero_at=400 | \
	Smooth rect1=7 rect2=7 repeat=3 | \
	Transp > t.H
	Vel n1=170 n2=500 n3=1 o1=0 o2=0 o3=0 d1=10 d2=10 d3=1 vc=2200 vr1=2100 z1=1400 const1=1 |  Transp > tt.H
	Add tt.H t.H | Pow pow=-1 | Rtoc > $@
	Rm t.H tt.H

Dat/moraData.H:
	python scripts/run.py vel=Vel/moraSlowTrue.H geom=Par/geom.json wave=Wav/moraWav.H out=$@

Vel/moraSlowBg$(BG).H:
	Vel n1=170 n2=500 n3=1 o1=0 o2=0 o3=0 d1=10 d2=10 d3=1 vc=$(BG) | \
	Transp | Pow pow=-1 | Rtoc > $@

Vel/moraSlowBg%_ext.H: Vel/moraSlowBg%.H
	Real < $< | Pad extend=1 n3out=73 end2=50 | Rtoc > $@
