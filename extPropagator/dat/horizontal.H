Input from a tty

Vel:   user@e50d93c8f46e   Tue Apr  9 05:01:16 2019
		sepstr_ready_out() : sets next in="stdin"
		Read  from param: #n1=200
		Read  from param: #n2=200
		Read  from param: #n3=1
		Read  from param: #big=0
		Read  from param: #dir=0
		Read  from param: #transp=0
		Read  from param: #slow=0
		Read  from param: #d1=0.5
		Read  from param: #d2=0.5
		Read  from param: #d3=1
		Read  from param: #o1=0
		Read  from param: #o2=0
		Read  from param: #o3=0
		n1=200
		n2=200
		n3=1
		d1=0.5
		d2=0.5
		d3=1
		o1=0
		o2=0
		o3=0
		Read  from param: #xa1=0
		Read  from param: #xb1=99.5
		Read  from param: #xa2=0
		Read  from param: #xb2=99.5
		Read  from param: #xa3=0
		Read  from param: #xb3=99.5
		Read  from param: #xa4=0
		Read  from param: #xb4=99.5
		Read  from param: #xa5=0
		Read  from param: #xb5=99.5
		Read  from param: #xa6=0
		Read  from param: #xb6=99.5
		Read  from param: #xa7=0
		Read  from param: #xb7=99.5
		Read  from param: #xa8=0
		Read  from param: #xb8=99.5
		Read  from param: #z1=50
		Read  from param: #z2=100
		Read  from param: #z3=100
		Read  from param: #z4=100
		Read  from param: #z5=100
		Read  from param: #z6=100
		Read  from param: #z7=100
		Read  from param: #z8=100
		Read  from param: #xasmo=0
		Read  from param: #xbsmo=99.5
		Read  from param: #z1smo=100
		Read  from param: #z2smo=100
		esize=4
		Read  from param: #const1=1
		Read  from param: #const2=0
		Read  from param: #const3=0
		Read  from param: #const4=0
		Read  from param: #const5=0
		Read  from param: #const6=0
		Read  from param: #const7=0
		Read  from param: #const8=0
		Read  from param: #alfa1=0
		Read  from param: #alfa2=0
		Read  from param: #alfa3=0
		Read  from param: #alfa4=0
		Read  from param: #alfa5=0
		Read  from param: #alfa6=0
		Read  from param: #alfa7=0
		Read  from param: #alfa8=0
		Read  from param: #vc=500
		Read  from param: #vr1=1000
		Read  from param: #vr2=500
		Read  from param: #vr3=500
		Read  from param: #vr4=500
		Read  from param: #vr4=500
		Read  from param: #vr5=500
		Read  from param: #vr6=500
		Read  from param: #hgrad=0
		Read  from param: #vr7=500
		Read  from param: #vr8=500
		Read  from param: #xan=0
		Read  from param: #zan=0
		Read  from param: #ranx=-1
		Read  from param: #ranz=-1
		Read  from param: #exan=1
		Read  from param: #dvan=0
		Read  from param: #widthan=1
		Read  from param: #rankx=-1
		Read  from param: #rankz=-1
		Read  from param: #dvank=0
		Read  from param: #exank=1
		Read  from param: #smooth1=1000
		Read  from param: #smooth2=500
		Read  from param: #smooth3=500
		Read  from param: #smooth4=500
		Read  from param: #smooth5=500
		Read  from param: #smooth6=500
		Read  from param: #smooth7=500
		Read  from param: #smooth8=500
		Read  from param: #v1smo=500
		Read  from param: #v2smo=500
		Read  from param: #exsmo=1
		Read  from param: #alfasmo=0
		Read  from param: #dzepsi=0
		Read  from param: #sin_period=5
		Read  from param: #sin_amp=5
		Read  from param: #sin_z=-0
		Read  from param: #sin_x=-0
		Read  from param: #sin_z_x0=0
		Read  from param: #sin_vel_up=0
		Read  from param: #sin_vel_down=0
		data_format="xdr_float"


Pow:   user@e50d93c8f46e   Tue Apr  9 05:01:16 2019
		sepstr_ready_out() : sets next in="stdin"
		esize=4
	#Using @SGN(file1)*file1^-1.000000 to produce real data

		#file1="in"
		data_format="xdr_float"


Pad:   user@e50d93c8f46e   Tue Apr  9 05:01:16 2019
		sets next: in="/home/user/project/scratch/horizontal.H@"
	 	n1=300 o1=-25 d1=0.5 label1="" 

		#beg=50
		#end=50
	 	n2=200 o2=0 d2=0.5 label2="" 

		data_format="xdr_float"

