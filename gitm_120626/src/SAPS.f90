 subroutine  SAPS_potential(MLAT,MLT,ASYH,potential)
 implicit none
 real ASYH
 real ph_i,peak,a1,b1,a2,b2,ag,bg,cg,F,G
 real,parameter:: a=-6.47,b=9.48,phi0=-3.4801, xa1=-0.7353, ya1=0.0027, &
 xb1=1.0262, yb1=-0.0023, xa2=0.091, ya2=-0.0023, xb2=0.1036, &
 yb2=0.0016, xga=0.04777, yga=0.00053, zga=0.00623, xgb=0.85795, ygb=0.006883, &
 zgb=0.061588, xgc=0.07129, ygc=-0.00795, zgc=-0.077798
 integer MLT,MLAT,deltamlt00
 real potential
 real::pi=3.1415926
 ! character (len=*), intent(in) :: dir
!integer :: i,j


 if (MLT >12) then
 deltamlt00=MLT-24
 else
 deltamlt00=MLT
 end if

 ph_i=(pi/12)*MLT+phi0
 peak=A+B*log10(ASYH)
 A1=Xa1+ya1*AsyH
 b1=Xb1+yb1*AsyH
 A2=xa2+ya2*AsyH
 b2=xb2+yb2*AsyH
 Ag=xga+yga*AsyH+zga*deltamlt00
 bg=xgb+ygb*AsyH+zgb*deltamlt00
 cg=xgc+ygc*AsyH+zgc*deltamlt00

  F=(A1*cos(ph_i)+b1*sin(ph_i))+(A2*cos(2*ph_i)+b2*sin(ph_i))
  G=Ag+Bg*MLAT+cg*MLAT*MLAT
 potential=-1*peak*F*G
 !           open(unit=115, file=dir//"/raaa.DAT")
 !write(115,*) potential(MLT,MLAT)


 end subroutine SAPS_potential