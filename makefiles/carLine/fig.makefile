ds = d1=0.125 d2=0.125 o1=0 o2=0 g2num=10 d2num=10
clip = bclip=1.5 eclip=1.7
barvals = minval=1.5 maxval=1.7
ratio = 0.5
scale = 0.8

Fig/%.v: Vel/%.H
	Window3d n3=1 min1=0 max1=128000 < $< | \
	Cabs | Pow pow=-1 | Transp | Scale dscale=1e-3 | \
	Grey color=j $(ds) \
	newclip=1 gainpanel=a wantscalebar=1 grid=0 gridfat=5 $(clip) \
	barlabel="Velocity [km/s]" title=" " $(barvals) \
	label1="Depth [mm]" label2="Distance [mm]" labelfat=5 axisfat=5 \
	screenratio=$(ratio) out=tmp.v > /dev/null
	vppen < tmp.v scale=$(scale) vpstyle=n > $@
	rm tmp.v

Fig/%.v: Inv/%.H
	Window3d n3=1 min1=0 max1=128000 < $< | \
	Cabs | Pow pow=-1 | Transp | Scale dscale=1e-3 | \
	Grey color=j $(ds) \
	newclip=1 gainpanel=a wantscalebar=1 grid=0 gridfat=5 $(clip) \
	barlabel="Velocity [km/s]" title=" " $(barvals) \
	label1="Depth [mm]" label2="Distance [mm]" labelfat=5 axisfat=5 \
	screenratio=$(ratio) out=tmp.v > /dev/null
	vppen < tmp.v scale=$(scale) vpstyle=n > $@
	rm tmp.v

SRCS = $(wildcard Fig/*.v)
OBJS = $(patsubst %.v,%.pdf,$(SRCS))
PDF: $(OBJS)
	make $^

%.pdf: %.v
	pstexpen $< $*.ps color=y fat=1.5 fatmult=1.5 invras=y background=black
	ps2pdf -dEPSCrop -dAutoFilterColorImages=false  -dColorImageFilter=/FlateEncode \
	-dAutoFilterGrayImages=false  -dGrayImageFilter=/FlateEncode  -dAutoFilterMonoImages=false \
	-dMonoImageFilter=/CCITTFaxEncode $*.ps $@
	rm $*.ps
