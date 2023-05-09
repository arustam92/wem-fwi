
grey_stack = label2="Distance [m]" label1="Depth[m]" labelfat=5 axisfat=5
grey_res = label2="Distance [m]" label1="Time[s]" labelfat=5 axisfat=5
grey_freq_slice = o1=3 d1=0.166 label1="Frequency [Hz]" label2="Velocity[m/s]" labelfat=5 axisfat=5
grey_par = newclip=1 gainpanel=a wantscalebar=1 barlabel="Velocity [km/s]" title=" "
clip_101 = bclip=1.7 eclip=2.2

Fig/%_fslice.v : Inv/%.H
	Window3d n1=1 min1=2500 n2=1 min2=800 < $< | \
	Real | Pow pow=-1 | Scale dscale=1e-3 | \
	Graph color=j $(grey_freq_slice) grid=1 movie=1 \
	min2=1.7 max2=2.2 title="$(title)" plotfat=10 \
	out=$@ > /dev/null

Fig/constSlowBg120-15Hz_final.v : Inv/constOnePass120-15Hz_eFWI_cmplx_eps_1_inv_mod.H Inv/constOnePass120-15Hz-1_eFWI_cmplx_eps_3_inv_mod.H Inv/constOnePass120-15Hz-1-3_eFWI_cmplx_eps_5_inv_mod.H Inv/constOnePass120-15Hz-1-3-5_eFWI_cmplx_eps_10_inv_mod.H Vel/constSlowBg120_ext_15Hz.H
	Cat3d $^ axis=4 | \
	Window3d n1=1 min1=2500 n2=1 min2=800 | \
	Real | Pow pow=-1 | Scale dscale=1e-3 | \
	Graph o1=1 d1=0.166 label1="Frequency [Hz]" label2="Velocity[km/s]" labelfat=5 axisfat=5 grid=1 \
	min2=1.6 max2=2.25 title=" " plotfat=10 \
	legend=1 legendtitle=" " legendfat=5 curvelabel="eps = 1":"eps = 3":"eps = 5":"eps = 10":"starting model" \
	out=$@ > /dev/null

Fig/%_panel.v : Inv/%.H
	Window3d n1=1 min1=2500 n2=1 min2=800 < $< | Real | Pow pow=-1 | \
	Grey color=j o1=3 d1=0.166 label1="Frequency [Hz]" label2="Iterations" labelfat=5 axisfat=5 \
	newclip=1 gainpanel=a wantscalebar=1 title="Real" plotfat=10 \
	out=t1.v grid=1 > /dev/null
	Window3d n1=1 min1=2500 n2=1 min2=800 < $< | Imag |\
	Grey color=j o1=3 d1=0.166 label1="Frequency [Hz]" label2="Iterations" labelfat=5 axisfat=5 \
	newclip=1 gainpanel=a wantscalebar=1 title="Imaginary" plotfat=10 \
	out=t2.v grid=1 > /dev/null
	vp_Movie t1.v t2.v > $@
	rm t1.v t2.v

Fig/%_fslice_dual.v : Inv/%.H
	Window3d n1=1 min1=2500 < $< | Real |\
	Grey color=j $(grey_freq_slice) $(grey_par) \
	title="$(title)" \
	out=$@ > /dev/null

Fig/%_stack.v : Inv/%.H
	Transp plane=23 < $< | \
	Real | Stack3d normalize=1 | Pow pow=-1 | Scale dscale=1e-3 | Transp | \
	Grey color=j $(grey_stack) $(grey_par) \
	bclip=$(bclip) eclip=$(eclip) title="$(title)" minval=$(bclip) maxval=$(eclip) \
	out=$@ > /dev/null

Fig/%_fwi.v : Inv/%.H
	Window3d n3=7 < $< | \
	Real | Pow pow=-1 | Scale dscale=1e-3 | Transp | \
	Grey color=j $(grey_stack) $(grey_par) \
	bclip=$(bclip) eclip=$(eclip) title="$(title)" minval=$(bclip) maxval=$(eclip) \
	out=$@ > /dev/null

##################################################### Objective function #####################################################
Fig/%_comp1.v : Inv/%_comp1.H
	Scale < $< | \
	Graph grid=1 min2=0 max2=1 title="Scaled obj. function" label1=" " plotfat=10 plotcol=6 \
	out=$@ > /dev/null
Fig/%_comp2.v : Inv/%_comp2.H
	Scale < $< | \
	Graph grid=1 min2=0 max2=1 title=" " label1=" " plotfat=10 plotcol=5 \
	out=$@ > /dev/null
Fig/%_cat.v : Fig/%_comp1.v Fig/%_comp2.v
	vp_Overlay $^ > $@

Fig/%_log.v : Inv/%.H
	Scale < $< > t.H
	Math file1=t.H exp="@LOG(file1)" > t2.H
	Graph < t2.H grid=1 out=$@ > /dev/null
	Rm t.H t2.H


Fig/%_log_comp1.v : Inv/%_comp1.H
	Scale < $< > t1.H
	Math file1=t1.H exp="@LOG(file1)/@LOG(10)" > t2.H
	Graph < t2.H grid=1 min2=-1 max2=0 title="Log obj. function" label1=" " plotfat=10 plotcol=6 \
	out=$@ > /dev/null
	Rm t1.H t2.H
Fig/%_log_comp2.v : Inv/%_comp2.H
	Scale < $< > t2.H
	Math file1=t2.H exp="@LOG(file1)/@LOG(10)" > t1.H
	Graph < t1.H grid=1 min2=-1 max2=0 title=" " label1=" " plotfat=10 plotcol=5 \
	out=$@ > /dev/null
	Rm t2.H t1.H
Fig/%_log_cat.v : Fig/%_log_comp1.v Fig/%_log_comp2.v
	vp_Overlay $^ > $@
##############################################################################################################################

Fig/t1.v : Inv/constOnePass120-15Hz_eFWI_cmplx_eps_1_obj_comp1.H Inv/constOnePass120-15Hz-1_eFWI_cmplx_eps_3_obj_comp1.H Inv/constOnePass120-15Hz-1-3_eFWI_cmplx_eps_5_obj_comp1.H Inv/constOnePass120-15Hz-1-3-5_eFWI_cmplx_eps_10_obj_comp1.H
	Cat3d axis=1 $^ | Scale | \
	Graph grid=1 min2=0 max2=1 title=" " label1=" " plotfat=10 plotcol=6 \
	out=$@ > /dev/null

Fig/t2.v : Inv/constOnePass120-15Hz_eFWI_cmplx_eps_1_obj_comp2.H Inv/constOnePass120-15Hz-1_eFWI_cmplx_eps_3_obj_comp2.H Inv/constOnePass120-15Hz-1-3_eFWI_cmplx_eps_5_obj_comp2.H Inv/constOnePass120-15Hz-1-3-5_eFWI_cmplx_eps_10_obj_comp2.H
	Cat3d axis=1 $^ | Scale | \
	Graph grid=1 min2=0 max2=1 title=" " label1=" " plotfat=10 plotcol=5 \
	out=$@ > /dev/null

Fig/constSlowBg120-15Hz_obj_final.v : Fig/t1.v Fig/t2.v
	vp_Overlay $^ > $@
	rm Fig/t1.v Fig/t2.v

Fig/%_res.v : Inv/%.H
	Window3d n3=1 f3=25 $(win) < $< | \
	Grey $(grey_res) pclip=100 gainpanel=a title="$(title)" \
	grid=1 min1=1 max1=3 out=$@ > /dev/null

Fig/%_first_res.v : Inv/%.H
	Window3d n3=1 f3=25 n4=1 < $< | \
	Grey $(grey_res) pclip=100 gainpanel=a title=" " \
	grid=1 min1=1 max1=3 out=$@ > /dev/null

Fig/%_pred.v : Inv/%.H
	Cat3d Dat/constData_onepass_cmplx.H $< axis=4 | \
	Window3d n3=1 f3=25 | \
	Grey $(grey_res) pclip=100 gainpanel=a title="$(title)" \
	out=$@ > /dev/null

################################################### REPORT ###################################################################
mod_dim = d2=0.01 d1=0.01
########## FWI ###########
Fig/wave_t.v: Wav/constWav_3-15.H
	Graph < $< grid=1 title=" " plotfat=10 labelfat=5 axisfat=5 \
	max1=1.6 g1num=.4 label1="Time[s]" label2="Amplitude" out=$@ > /dev/null

Fig/wave_f.v: Wav/constWav_3-15.H
	Spectra < $< | \
	Graph grid=1 title=" " plotfat=10 labelfat=5 axisfat=5 \
	min1=2 max1=20 g1num=5 label1="Frequency[Hz]" label2="Amplitude" out=$@ > /dev/null

Fig/fwi_inv_mod.v: Inv/fwi_constOnePass120-3-8Hz_inv_mod.H
	Window3d max2=1700 < $< | Real | Pow pow=-1 | Scale dscale=1e-3 | Transp | \
	Grey color=j $(mod_dim) grid=1 newclip=1 title=" " wantscalebar=1 eclip=2.2 bias=1.8 barlabel="Velocity[km/s]" \
	label2="Distance[km]" label1="Depth[km]" labelfat=5 axisfat=5 \
	out=$@ > /dev/null

Fig/fwi_obj.v: Inv/fwi_constOnePass120-3-8Hz_obj.H
	Scale < $< | \
	Graph min2=0 max2=1 grid=1 title=" " plotfat=10 labelfat=5 axisfat=5 \
	label2="Scaled obj. function" label1="Iterations" \
	out=$@ > /dev/null

########## RESIDUAL MATCHING ###########
Fig/first_res.v: Inv/auglag2_constOnePass120-3-15Hz_eFWI_cmplx_eps_1_residual_comp1.H
	Window3d n4=1 n3=1 f3=25 min1=1 max1=3 < $< | \
	Grey pclip=100 grid=1 labelfat=5 axisfat=5 d2=.01 title=" " \
	label1="Time[s]" label2="Distance[km]" out=$@ > /dev/null
Fig/first_grad.v: Inv/auglag2_constOnePass120-3-15Hz_eFWI_cmplx_eps_1_gradient.H
	Window3d n1=1 min1=2500 max2=1600 n4=1 < $< | Real | \
	Grey o2=3 pclip=99 grid=1 labelfat=5 axisfat=5 d2=.166 d1=0.01 title=" " \
	label1="Depth[km]" label2="Frequency[Hz]" out=$@ > /dev/null
Fig/first_grad_ext.v: Inv/auglag2_constOnePass120-3-15Hz_eFWI_cmplx_eps_1_gradient.H
	Window3d n1=1 min1=2500 n4=1 < $< | Real | \
	Grey o2=3 pclip=99 grid=1 labelfat=5 axisfat=5 d2=.166 d1=0.01 title=" " \
	label1="Depth[km]" label2="Frequency[Hz]" out=$@ > /dev/null
Fig/first_grad_tau.v: Inv/auglag2_constOnePass120-3-15Hz_eFWI_cmplx_eps_1_gradient_tau.H
	Window3d n2=1 min2=2500 max3=1600 min1=2 max1=4 < $< | Transp | \
	Grey o2=-1 pclip=100 grid=1 labelfat=5 axisfat=5 d1=.01 title=" " \
	label1="Depth[km]" label2="Time[s]" out=$@ > /dev/null
Fig/first_grad_tau_ext.v: Inv/auglag2_constOnePass120-3-15Hz_eFWI_cmplx_eps_1_gradient_tau.H
	Window3d n2=1 min2=2500 min1=2 max1=4 < $< | Transp | \
	Grey o2=-1 pclip=100 grid=1 labelfat=5 axisfat=5 d1=.01 title=" " \
	label1="Depth[km]" label2="Time[s]" out=$@ > /dev/null
########## RESIDUAL MATCHING ###########

# SRCS = $(wildcard Fig/*.v)
SRCS = Fig/wave_t.v Fig/wave_f.v Fig/fwi_inv_mod.v Fig/fwi_obj.v \
Fig/first_res.v Fig/first_grad.v Fig/first_grad_ext.v Fig/first_grad_tau.v Fig/first_grad_tau_ext.v

OBJS = $(patsubst %.v,%.pdf,$(SRCS))
PDF: $(OBJS)
	make $^

%.pdf: %.v
	pstexpen $< $*.ps color=y fat=1.5 fatmult=1.5 invras=y background=black
	ps2pdf -dEPSCrop -dAutoFilterColorImages=false  -dColorImageFilter=/FlateEncode \
	-dAutoFilterGrayImages=false  -dGrayImageFilter=/FlateEncode  -dAutoFilterMonoImages=false \
	-dMonoImageFilter=/CCITTFaxEncode $*.ps $@
	rm $*.ps
