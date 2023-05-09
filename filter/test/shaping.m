


$(T)/test.H: $(B)/Shaping.x
	@mkdir test -p 
	Gauss n1=50 o1=0 n2=2 > $(T)/g1.H
	Gauss n1=50 o1=-20 n2=2 > $(T)/g2.H
	$< < $(T)/g1.H wave=$(T)/g2.H neglag=1 poslag=2 > $@

$(T)/test2.H: $(B)/Shaping.x
	@mkdir test -p 
	Wavelet wavelet=ricker2 phase=minimum n1=200 o1=0 d1=0.002 > $(T)/w1.H
	Wavelet wavelet=ricker2 phase=180 n1=200 o1=-0.1 tdelay=0.1 d1=0.002 > $(T)/w2.H
	Cat3d $(T)/w1.H $(T)/w1.H axis=2 > $(T)/w11.H
	Cat3d $(T)/w2.H $(T)/w2.H axis=2 > $(T)/w22.H
	$< < $(T)/w11.H wave=$(T)/w22.H neglag=70 poslag=50 eps=0.1 > $@

plot: $(T)/test2.H $(B)/Conv.x
	$(B)/Conv.x < $(T)/w11.H filter=$< | Graph grid=1 | Tube & 

$(T)/test3.H: $(B)/Conv.x
	Spike n1=100 k1=50 n2=2 > $(T)/sp1.H
	Spike n1=100 k1=80 n2=2 > $(T)/sp2.H
	$< < $(T)/sp1.H filter=$(T)/sp2.H > $@
	