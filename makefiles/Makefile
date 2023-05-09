
Fig/const_fwi_2000.v:
	Window3d n3=1 < Inv/const_fwi_2000_model.H | Cabs | Pow pow=-1 | Transp | Scale dscale=1e-3 | \
	Grey gainpanel=a color=j newclip=1 wantscalebar=1 out=$@ title="FWI"> /dev/null

Fig/const_fwi_2100.v:
	Window3d n3=1 < Inv/const_fwi_model.H | Cabs | Pow pow=-1 | Transp | Scale dscale=1e-3 | \
	Grey gainpanel=a color=j newclip=1 wantscalebar=1 out=$@ title="FWI"> /dev/null

Fig/%.v: Inv/%.H
	Window3d n3=1 n1=1 min1=5000 < $< | Cabs | Pow pow=-1 | Scale dscale=1e-3 | \
	Grey gainpanel=a color=j newclip=1 wantscalebar=1 out=$@ title="TFWI at x=5000"> /dev/null
