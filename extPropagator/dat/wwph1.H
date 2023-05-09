Input from a tty

Vel:   user@e50d93c8f46e   Wed Apr 17 00:00:16 2019
		sepstr_ready_out() : sets next in="stdin"
		Read  from param: #n1=30
		Read  from param: #n2=200
		Read  from param: #n3=1
		Read  from param: #big=0
		Read  from param: #dir=0
		Read  from param: #transp=0
		Read  from param: #slow=0
		Read  from param: #d1=10
		Read  from param: #d2=10
		Read  from param: #d3=1
		Read  from param: #o1=0
		Read  from param: #o2=0
		Read  from param: #o3=0
		n1=30
		n2=200
		n3=1
		d1=10
		d2=10
		d3=1
		o1=0
		o2=0
		o3=0
		Read  from param: #xa1=0
		Read  from param: #xb1=1990
		Read  from param: #xa2=0
		Read  from param: #xb2=1990
		Read  from param: #xa3=0
		Read  from param: #xb3=1990
		Read  from param: #xa4=0
		Read  from param: #xb4=1990
		Read  from param: #xa5=0
		Read  from param: #xb5=1990
		Read  from param: #xa6=0
		Read  from param: #xb6=1990
		Read  from param: #xa7=0
		Read  from param: #xb7=1990
		Read  from param: #xa8=0
		Read  from param: #xb8=1990
		Read  from param: #z1=250
		Read  from param: #z2=300
		Read  from param: #z3=300
		Read  from param: #z4=300
		Read  from param: #z5=300
		Read  from param: #z6=300
		Read  from param: #z7=300
		Read  from param: #z8=300
		Read  from param: #xasmo=0
		Read  from param: #xbsmo=1990
		Read  from param: #z1smo=300
		Read  from param: #z2smo=300
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
		Read  from param: #vc=1000
		Read  from param: #vr1=2000
		Read  from param: #vr2=1000
		Read  from param: #vr3=1000
		Read  from param: #vr4=1000
		Read  from param: #vr4=1000
		Read  from param: #vr5=1000
		Read  from param: #vr6=1000
		Read  from param: #hgrad=0
		Read  from param: #vr7=1000
		Read  from param: #vr8=1000
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
		Read  from param: #smooth1=2000
		Read  from param: #smooth2=1000
		Read  from param: #smooth3=1000
		Read  from param: #smooth4=1000
		Read  from param: #smooth5=1000
		Read  from param: #smooth6=1000
		Read  from param: #smooth7=1000
		Read  from param: #smooth8=1000
		Read  from param: #v1smo=1000
		Read  from param: #v2smo=1000
		Read  from param: #exsmo=1
		Read  from param: #alfasmo=0
		Read  from param: #dzepsi=0
		Read  from param: #sin_period=100
		Read  from param: #sin_amp=100
		Read  from param: #sin_z=-0
		Read  from param: #sin_x=-0
		Read  from param: #sin_z_x0=0
		Read  from param: #sin_vel_up=0
		Read  from param: #sin_vel_down=0
		data_format="xdr_float"


Pow:   user@e50d93c8f46e   Wed Apr 17 00:00:16 2019
		sepstr_ready_out() : sets next in="stdin"
		esize=4
	#Using @SGN(file1)*file1^-1.000000 to produce real data

		#file1="in"
		data_format="xdr_float"


Pad:   user@e50d93c8f46e   Wed Apr 17 00:00:16 2019
		sepstr_ready_out() : sets next in="stdin"
	 	n1=30 o1=0 d1=10 label1="" 

	 	n2=300 o2=-500 d2=10 label2="" 

		#beg=50
		#end=50
		data_format="xdr_float"


Transp:   user@e50d93c8f46e   Wed Apr 17 00:00:16 2019
		sets next: in="/home/user/project/scratch/slow.H@"
Transposing plane 1 and 2 
Using 1000 megabytes 
		junkME=0
	n2=30  o2=0.000000  d2=10.000000   unit2="Undefined"
		gff="-1"
		hff="-1"
hff=-1 gff=-1  n9=1 
	n1=300  o1=-500.000000  d1=10.000000   unit1="Undefined"
		hff="-1"
		gff="-1"
		esize=4
		data_format="xdr_float"


../bin/checkWEM.x:   user@e50d93c8f46e   Wed Apr 17 00:00:24 2019
		sets next: in="/home/user/project/scratch/ph1.H@"
	 	n1=1500 o1=0 d1=0.002 label1="Undefined" 

	 	n2=200 o2=0 d2=10 label2="Undefined" 

	 	n3=2 o3=0 d3=10 label3="Undefined" 

	 	n4=1 o4=0 d4=1 label4="none" 

	 	n5=1 o5=0 d5=1 label5="none" 

	 	n6=1 o6=0 d6=1 label6="none" 

	 	n7=1 o7=0 d7=1 label7="none" 

	 	n8=1 o8=0 d8=1 label8="none" 

		data_format="xdr_float"



Window3d:   user@e50d93c8f46e   Wed Apr 17 00:03:27 2019
		sets next: in="/home/user/project/scratch/wph1.H@"
#Window #f1=0  #j1=1 #n1=750
#Window #f2=0  #j2=1 #n2=200
#Window #f3=0  #j3=1 #n3=2
		junkME=0
	n2=200  o2=0.000000  d2=10.000000   label2="Undefined"   unit2="Undefined"
	n3=2  o3=0.000000  d3=10.000000   label3="Undefined"   unit3="Undefined"
		gff="-1"
		hff="-1"
hff=-1 gff=-1  n9=1 
	n1=750  o1=0.000000  d1=0.002000   label1="Undefined"   unit1="Undefined"
		hff="-1"
		gff="-1"
		esize=4
		data_format="xdr_float"



Pad:   user@e50d93c8f46e   Wed Apr 17 00:12:15 2019
		sets next: in="/home/user/project/scratch/wwph1.H@"
	 	n1=1500 o1=0 d1=0.002 label1="" 

		#end=750
	 	n2=200 o2=0 d2=10 label2="" 

	 	n3=2 o3=0 d3=10 label3="" 

		data_format="xdr_float"

