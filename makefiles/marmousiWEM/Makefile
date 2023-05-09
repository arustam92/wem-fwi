Fig/right_obj.v:
	Scale < Inv/marm_3-10_obj.H | \
	Graph min2=0 max2=1 grid=1 max1=50 min1=0 \
	title=" " label1="Iteration" label2="Scaled obj. function" \
	plotfat=10 labelfat=5 axisfat=5 gridfat=5 plotcol=6 \
	out=$@ > /dev/null

Fig/wrong_obj.v:
		Scale < Inv/wrong_marm_3-10_obj.H | \
		Graph min2=0 max2=1 grid=1 max1=50 min1=0 \
		title=" " label1=" " label2=" " \
		plotfat=10 labelfat=5 axisfat=5 gridfat=5 plotcol=5 \
		out=$@ > /dev/null

Fig/marm_obj_compare.v: Fig/wrong_obj.v Fig/right_obj.v
	vp_Overlay $^ > $@

Fig/marm_inv_mod_right.v:
	Cat3d Inv/marm_3-10_inv_mod.H Inv/wrong_marm_3-10_inv_mod.H Vel/MarmBg.H axis=4 | \
	Window3d n3=1 n4=1 min1=600 max1=6100 | \
	Cabs | Pow pow=-1 | Transp | Scale dscale=1e-3 | \
	Grey color=j d1=0.012 d2=0.012 o1=0 o2=0.7 g2num=0.5 d2num=0.5 d1num=0.5 g1num=0.5 \
	newclip=1 gainpanel=a wantscalebar=1 grid=1 gridfat=5 bclip=1.5 eclip=6.0 \
	barlabel="Velocity [km/s]" minval=1.5 maxval=6.0 title=" " \
	label1="Depth [km]" label2="Distance [km]" labelfat=5 axisfat=5 \
	screenratio=0.4 out=tmp.v > /dev/null
	vppen < tmp.v scale=0.75 vpstyle=n > $@
	rm tmp.v

Fig/marm_inv_mod_wrong.v:
	Cat3d Inv/marm_3-10_inv_mod.H Inv/wrong_marm_3-10_inv_mod.H Vel/MarmBg.H axis=4 | \
	Window3d n3=1 n4=1 f4=1 min1=600 max1=6100 | \
	Cabs | Pow pow=-1 | Transp | Scale dscale=1e-3 | \
	Grey color=j d1=0.012 d2=0.012 o1=0 o2=0.7 g2num=0.5 d2num=0.5 d1num=0.5 g1num=0.5 \
	newclip=1 gainpanel=a wantscalebar=1 grid=1 gridfat=5 bclip=1.5 eclip=6.0 \
	barlabel="Velocity [km/s]" minval=1.5 maxval=6.0 title=" " \
	label1="Depth [km]" label2="Distance [km]" labelfat=5 axisfat=5 \
	screenratio=0.4 out=tmp.v > /dev/null
	vppen < tmp.v scale=0.75 vpstyle=n > $@
	rm tmp.v

Fig/marm_bg.v:
	Cat3d Inv/marm_3-10_inv_mod.H Inv/wrong_marm_3-10_inv_mod.H Vel/MarmBg.H axis=4 | \
	Window3d n3=1 n4=1 f4=2 min1=600 max1=6100 | \
	Cabs | Pow pow=-1 | Transp | Scale dscale=1e-3 | \
	Grey color=j d1=0.012 d2=0.012 o1=0 o2=0.7 g2num=0.5 d2num=0.5 d1num=0.5 g1num=0.5 \
	newclip=1 gainpanel=a wantscalebar=1 grid=1 gridfat=5 bclip=1.5 eclip=6.0 \
	barlabel="Velocity [km/s]" minval=1.5 maxval=6.0 title=" " \
	label1="Depth [km]" label2="Distance [km]" labelfat=5 axisfat=5 \
	screenratio=0.4 out=tmp.v > /dev/null
	vppen < tmp.v scale=0.75 vpstyle=n > $@
	rm tmp.v


%.pdf: %.v
	pstexpen $< $*.ps color=y fat=1.5 fatmult=1.5 invras=y
	ps2pdf -dEPSCrop -dAutoFilterColorImages=false  -dColorImageFilter=/FlateEncode \
	-dAutoFilterGrayImages=false  -dGrayImageFilter=/FlateEncode  -dAutoFilterMonoImages=false \
	-dMonoImageFilter=/CCITTFaxEncode $*.ps $@
	rm $*.ps

SRCS = $(wildcard Fig/*.v)
OBJS = $(patsubst %.v,%.pdf,$(SRCS))
PDF: $(OBJS)
	make $^
