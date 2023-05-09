ds = d1=0.01 d2=0.025 o1=0 o2=0 g2num=0.5 d2num=0.5
clip = bclip=2 eclip=2.5
barvals = minval=2 maxval=2.5
ratio = 0.2
scale = 0.53



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
