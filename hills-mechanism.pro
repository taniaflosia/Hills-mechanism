
pro hills_mechanism,mass,abin,vtan

;============================================================================================================;
; Initial Setting
;============================================================================================================;

g	= 6.673*1.d-11	  		;konstanta gravitasi [Nm2/kg2] atau [m3/s2kg]
au	= 1.49597871*1.d11		;satuan astronomi [m]
sm	= 1.9891*1.d30			;massa matahari [kg]
sr  	= 6.9550826*1.d8		;radius matahari [m]
pc	= 206265*AU			;jarak parsec [m]
mnt 	= double(60.)			;satu menit [s]
hr	= double(3600.)			;satu jam [s]
day 	= double(24.*3600.)		;satu hari [s]
yr	= double(365.*24.*3600.)	;satu tahun [s]
myr	= 1.d6*365.25*24*3600		;satu juta tahun [s]

n = 3					;jumlah benda
m = mass*sm				;massa bintang ganda dan black hole [kg]	
r = dblarr(n,3)			;posisi n benda dalam 3 dimensi
v = dblarr(n,3)			;kecepatan
a = dblarr(n,3)			;percepatan
jk= dblarr(n,3)			;jerk

old_r = dblarr(n,3)
old_v = dblarr(n,3)
old_a = dblarr(n,3)
old_j = dblarr(n,3)

;============================================================================================================;
; Time Input
;============================================================================================================;

t  = 0.
dt = mnt
tf = 100*yr

;============================================================================================================;
; Initial Position, Eccentricity, Semimajor Axis
;============================================================================================================;

ebin = 0.6		;eksentrisitas sistem bintang ganda
abin = abin*au		;setengah sumbu panjang bintang ganda [m]
dbin = 0.01*pc		;posisi awal bintang ganda [m]

r[0,*] = [dbin,0,0]
r[1,*] = [dbin+abin*(1+ebin),0,0]
r[2,*] = [0,0,0]

file="Output 1"+string(vtan)
set_plot,'ps'
device,filename=file+".ps",xsize=6,ysize=6,/inch,/color,/encap
loadct,39

plot,title="Plot Interaksi Bintang Ganda dan Black Hole",xtitle="Sumbu X Galaksi [AU]",$
	ytitle="Sumbu Z Galaksi [AU]",r[*,0]/AU,r[*,2]/AU,psym=2,xrange=[-2000,6000],$
	yrange=[-4000,4000],/iso,background=255,color=0

;============================================================================================================;
; Binary Binding Energy, Initial Velocity of Star & System
;============================================================================================================;

miu   = m[0]*m[1]/(m[0]+m[1])
Ebind = -G*m[0]*m[1]/(2*abin)
print,"Energi ikat awal		=",Ebind

rbin = r[1,*]-r[0,*]					;jarak relatif, vektor
rrel = sqrt(total(rbin*rbin))				;jarak relatif, skalar
vrel = sqrt((2*Ebind/miu)+(2*G*m[0]*m[1]/(miu*rrel)))	;kecepatan relatif, skalar

V1c	   = 49.827757		;[km/s]
V2c	   = 177.70573 		;[km/s]

v[0,*] = [0, vrel,vtan]
v[1,*] = [0,-vrel,vtan]
v[2,*] = [0,0,0]

;============================================================================================================;
; Energy Calculation - Different formula with Ebin
;============================================================================================================;

;Energi Ikat Sistem Bintang Ganda

epot=0.d0
ekin=0.d0
for i=0,1 do begin
	for j=i+1,n-2 do begin
		rji=dblarr(3)
		for k=0,2 do rji[k]=r[j,k]-r[i,k]
		r2=0.d0
		for k=0,2 do r2=r2+(rji[k]*rji[k])
		epot=epot-(G*m[i]*m[j]/sqrt(r2))
	endfor
	for k=0,2 do begin
		v[i]=dblarr(2)
		v[i]=sqrt(total(v[i,k]*v[i,k]))
	endfor
	for j=i+1,n-2 do begin
		if i eq j then continue
		mtot=m[i]+m[j]
		vcm=(m[i]*v[i]+m[j]*v[j])/mtot
	endfor
	v[i]=v[i]-vcm
	ekin=ekin+(0.5*miu*v[i]*v[i])
endfor

Ebin_in = ekin+epot
Rstar_1 = sqrt(total(r[0,*]*r[0,*]))-sqrt(total(r[2,*]*r[2,*]))
Rstar_2 = sqrt(total(r[1,*]*r[1,*]))-sqrt(total(r[2,*]*r[2,*]))
Vstar_1 = sqrt(total(v[0,*]*v[0,*]))-sqrt(total(v[2,*]*v[2,*]))
Vstar_2 = sqrt(total(v[1,*]*v[1,*]))-sqrt(total(v[2,*]*v[2,*]))

openw,lun,file+".dat",/get_lun
printf,lun,"# Data Timestep [Day], Ebin [1.d40 Joule], Rstar 1 [AU], Rstar 2 [AU], Vstar 1 [km/s], Vstar 2 [km/s]"
printf,lun,"# Batas Minimum Kecepatan 		= ",V1c," km/s"
printf,lun,"# Batas Maksimum Kecepatan		= ",V2c," km/s"
printf,lun,"# Nilai Kecepatan Tangensial	= ",Vtan/1000.," km/s"
printf,lun,"# Massa Komponen Bintang Ganda	= ",m[0]/SM," M sun"
printf,lun,"# Setengah Sumbu Panjang		= ",abin/AU," AU"
format='(f12.2,3x,f9.2,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.2)'
printf,lun,format=format,t/yr,Ebin_in/1.d40,Rstar_1/AU,Rstar_2/AU,Vstar_1/1.d3,Vstar_2/1.d3

;============================================================================================================;
; Initial Acceleration, Initial Jerk
;============================================================================================================;

for i=0,n-1 do begin
	for k=0,2 do a[i,k]=0.d0
	for k=0,2 do jk[i,k]=0.d0
endfor

for i=0L,n-1 do begin
	for j=i+1,n-1 do begin
	rji=dblarr(3)
	vji=dblarr(3)
		for k=0,2 do rji[k] = r[j,k]-r[i,k]
		for k=0,2 do vji[k] = v[j,k]-v[i,k]
		r2=0.d0
		for k=0,2 do r2=r2+(rji[k]*rji[k])
		r3 = r2*sqrt(r2)
		rv = 0.d0
		for k=0,2 do rv = rv+(rji[k]*vji[k])
		rv=rv/r2
		for k=0,2 do begin
		a[i,k] = a[i,k]+(G*m[j]*rji[k]/r3)			;persamaan newton 6.6
		a[j,k] = a[j,k]-(G*m[i]*rji[k]/r3)
		jk[i,k] = jk[i,k]+(G*m[j]*(vji[k]-3*rv*rji[k])/r3)	;turunan persamaan newton 6.4
		jk[j,k] = jk[j,k]-(G*m[i]*(vji[k]-3*rv*rji[k])/r3)
		endfor
	endfor
endfor

;============================================================================================================;
; Integrasi Hermite
;============================================================================================================;

iter=0


REPEAT BEGIN

	for i=0,n-1 do begin
		for k=0,2 do begin
		old_r[i,k]=r[i,k]
		old_v[i,k]=v[i,k]
		old_a[i,k]=a[i,k]
		old_j[i,k]=jk[i,k]
		r[i,k]=r[i,k]+(v[i,k]*dt)+(a[i,k]*dt*dt/2)+(jk[i,k]*dt*dt*dt/6)		;persamaan 6.13, deret taylor
		v[i,k]=v[i,k]+(a[i,k]*dt)+(jk[i,k]*dt*dt/2)				;persamaan 6.14, deret taylor
		endfor
	endfor

	for i=0,n-1 do begin
		for k=0,2 do a[i,k]=0.d0
		for k=0,2 do jk[i,k]=0.d0
	endfor

	;Energi Ikat Sistem Bintang Ganda

	epot=0.d0
	ekin=0.d0
	for i=0,1 do begin
		for j=i+1,n-2 do begin
			rji=dblarr(3)
			for k=0,2 do rji[k]=r[j,k]-r[i,k]
			r2=0.d0
			for k=0,2 do r2=r2+(rji[k]*rji[k])
			epot=epot-(G*m[i]*m[j]/sqrt(r2))
		endfor
		for k=0,2 do begin
		v[i]=dblarr(2)
		v[i]=sqrt(total(v[i,k]*v[i,k]))
		endfor
		for j=i+1,n-2 do begin
			if i eq j then continue
			mtot=m[i]+m[j]
			vcm=(m[i]*v[i]+m[j]*v[j])/mtot
		endfor
		v[i]=v[i]-vcm
		ekin=ekin+(0.5*miu*v[i]*v[i])
	endfor

	Ebin_out=ekin+epot

	for i=0,n-1 do begin
		for j=i+1,n-1 do begin
			rji=dblarr(3)
			vji=dblarr(3)
			for k=0,2 do begin
			rji[k]=r[j,k]-r[i,k]
			vji[k]=v[j,k]-v[i,k]
			endfor
			r2=0.d0
			for k=0,2 do r2=r2+(rji[k]*rji[k])
			r3=double(r2*sqrt(r2))
			rv=0.d0
			for k=0,2 do rv=rv+rji[k]*vji[k]
			rv=rv/r2
			for k=0,2 do begin
				a[i,k]=a[i,k]+(G*m[j]*rji[k]/r3)			;persamaan newton 6.6
				a[j,k]=a[j,k]-(G*m[i]*rji[k]/r3)
				jk[i,k]=jk[i,k]+(G*m[j]*(vji[k]-3*rv*rji[k])/r3)	;turunan persamaan newton 6.4
				jk[j,k]=jk[j,k]-(G*m[i]*(vji[k]-3*rv*rji[k])/r3)
			endfor
		endfor
	endfor

	for i=0,n-1 do begin
		for k=0,2 do begin
		v[i,k]=old_v[i,k]+(old_a[i,k]+a[i,k])*dt/2+(old_j[i,k]-jk[i,k])*dt*dt/12	;persamaan hermite
		r[i,k]=old_r[i,k]+(old_v[i,k]+v[i,k])*dt/2+(old_a[i,k]-a[i,k])*dt*dt/12		;persamaan hermite
		endfor
	endfor

	t = t+dt

	;Posisi dan Kecepatan Relatif Terhadap Black Hole

	Rstar_1 = sqrt(total(r[0,*]*r[0,*]))-sqrt(total(r[2,*]*r[2,*]))
	Rstar_2 = sqrt(total(r[1,*]*r[1,*]))-sqrt(total(r[2,*]*r[2,*]))
	Vstar_1 = sqrt(total(v[0,*]*v[0,*]))-sqrt(total(v[2,*]*v[2,*]))
	Vstar_2 = sqrt(total(v[1,*]*v[1,*]))-sqrt(total(v[2,*]*v[2,*]))

	;Jarak Relatif Terhadap Black Hole

	rplot=dblarr(n,3)
	rplot[0,*]=r[0,*]-r[2,*]
	rplot[1,*]=r[1,*]-r[2,*]
	rplot[2,*]=r[2,*]-r[2,*]

	iter=iter+1
	time=double(t/yr)
	if iter ge 1440 then begin
		printf,lun,format=format,time,Ebin_out/1.d40,Rstar_1/AU,Rstar_2/AU,Vstar_1/1.d3,Vstar_2/1.d3
		oplot,rplot[*,0]/AU,rplot[*,2]/AU,psym=3,color=0
		iter=0
	endif

ENDREP UNTIL t GE tf

free_lun, lun

;write_gif,file+".gif",tvrd()
device,/close
set_plot,'x'
print,"End"

end

;============================================================================================================;
; End of Program
;============================================================================================================;
