

$(T)/kol.s.H: $(B)/Kolmogorov.x
	Pad end1=0 < $(T)/r2.H > $(T)/rp2.H
	$< < $(T)/r1.H eps=0.0 > $(T)/kol.f.H
	Filter < $(T)/rp2.H filter=$(T)/kol.f.H > $@
	Rm $(T)/rp2.H 
