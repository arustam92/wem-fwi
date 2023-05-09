
valgrind = valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --track-origins=yes --log-file=val.txt

ph1.H: slow.H wave.H ${B}/checkWEM.x
	${B}/checkWEM.x < slow.H wave=wave.H par=geometry.par > $@

adjWEM1.H: slow.H
	$(B)/check_adjWEM.x < slow.H data=ph1.H par=geometry.par > $@

dotWEM: slow.H wave.H
	$(B)/dotWEM.x < slow.H wave=wave.H par=geometry.par > /dev/null


ph2.H: slow2.H wave.H
	$(B)/checkWEM.x < slow2.H wave=wave.H par=geometry.par wfld=wfld.H > $@


born2.H: slow.H wave.H dslow2.H
	$(B)/checkBorn.x < slow.H dvel=dslow2.H wave=wave.H par=geometry.par > $@
