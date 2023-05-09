B:=../bin

wave.H:
	Spike n1=1500 d1=0.002 o1=0 k1=100 | Bandpass flo=10 fhi=20 > $@
	#Wavelet wavelet=ricker2 phase=minimum tdelay=0 o1=0 d1=0.002 n1=1500 fund=5 | \
	#Window3d n1=750 | Pad end1=750 > $@

slow.H:
	Vel n1=40 n2=200 n3=1 o1=0 o2=0 o3=0 d1=10 d2=10 d3=1 vc=1000 vr1=2000 z1=250 const1=1 | \
        Pow pow=-1 | \
        Pad beg2=50 end2=50 extend=1 | Transp > $@


bg2.H: slow.H
	Smooth < $< rect1=5 rect2=5 repeat=3 > $@

slow2.H: slow.H
	Gauss n1=30 n2=200 n3=1 o1=0 o2=0 o3=0 d1=10 d2=10 d3=1 max1=150 max2=100 \
	maxvel=-1000 velback=0 var=2000 ngauss=1 | Pad beg2=50 end2=50 extend=1 | \
	Transp > t.H
	Vel n1=30 n2=200 n3=1 o1=0 o2=0 o3=0 d1=10 d2=10 d3=1 vc=2000 vr1=2500 z1=250 const1=1 |\
	Pad beg2=50 end2=50 extend=1 | Transp > tt.H
	Add tt.H t.H | Pow pow=-1 > $@
	rm t.H tt.H
bg.H:
	Vel n1=40 n2=200 n3=1 o1=0 o2=0 o3=0 d1=10 d2=10 d3=1 vc=1000 const1=1 | \
	Pow pow=-1 | Pad beg2=50 end2=50 extend=1 | Transp > $@


horizontal.H:
	Vel n1=200 n2=200 n3=1 o1=0 o2=0 o3=0 d1=0.5 d2=0.5 d3=1 vc=500 vr1=1000 z1=50 const1=1 | \
	Pow pow=-1 | Pad beg1=50 end1=50 extend=1 > $@


include ph.makefile
include born.makefile
