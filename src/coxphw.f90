SUBROUTINE WEIGHTEDCOX(cards, parms, IOARRAY, DFBETA)
!DEC$ ATTRIBUTES DLLEXPORT :: weightedcox

IMPLICIT DOUBLE PRECISION (A-H,O-Z)  


!This version allows for offset (1st col of X matrix), control by 16th entry of parms
!no pre-sorting necessary if score weights are available in each of several tied event rows

real*8, dimension (16) :: parms
integer N, IP,JCODE,Irobust,ISEP,ITER,IMAXIT,IMAXHS,Ioffset
real*8, dimension (int(parms(1))) :: BX, T1, t2, TMSF, Doffset, doffset2
real*8, dimension (int(parms(1)),int(parms(2))) :: X, XMSF
!real*8, dimension (int(parms(1)),int(parms(2))) :: X, XMSF, bresx
real*8, dimension (int(parms(2)+parms(14))) :: B, B0, FD, TDERR,BMSF,zw1, xx, xfd, yy, dinit 
real*8, dimension (int(parms(2)+parms(14)),int(parms(2)+parms(14))) :: SD, VM, vmlw, vmls, WK, fish, ainv, vminv
!integer, dimension (int(parms(1))) :: ibresc, IC, ICMSF, patid
integer, dimension (int(parms(1))) :: IC, ICMSF, patid
integer, dimension (int(parms(2)+parms(14))) :: IFLAG
real*8, dimension (int(parms(15)), int(parms(2)+parms(14))) :: dfbeta
real*8, dimension (int(parms(1)),int((2*parms(2)+4+2*(parms(14))))) :: cards
real*8, dimension (int((3+3*(parms(2)+parms(14)))),int((parms(2)+parms(14)))) :: IOARRAY
real*8, dimension (14) :: DER, EREST
logical, dimension (int(parms(2)+parms(14)),int(parms(2)+parms(14))) :: mask
real*8, dimension (int(parms(1)),int(parms(14)+1)) :: ft
integer ngv, ntde
integer, dimension (int(parms(14)+1)) :: ftmap
real*8, dimension (int(parms(1)), int(parms(14)+parms(2))) :: score_weights

INTRINSIC DABS, DSQRT               

! cards contains in its last column the ID numbers of the patients (starting from 1)
! parms(15) contains the number of patients (the max of patid)

!open(unit=6, file="fgcssan.txt")

! ntde = parms(14)
ifail=0              
N=parms(1)
IP=parms(2)
Irobust=parms(3)
imaxit=Parms(4)
imaxhs=parms(5)
step=parms(6)
xconv=parms(7)
gconv=parms(8)

ngv=parms(13)
ntde=parms(14)
numbpatients=int(parms(15))
ioffset=int(parms(16))
patid=int(cards(:,int(ioffset+(2*parms(2)+4+2*(parms(14))))))

parms(10)=-11

ilastlike=0
!ilastlike=1   ! would always compute correct variance (lin-sasieni), but needs more iterations

dinit=ioarray(2,:)
iflag=ioarray(1,:)

t1=cards(:,ip+1+ioffset)
t2=cards(:,ip+2+ioffset)
ic=cards(:,ip+3+ioffset)
score_weights=cards(:,(ip+4+ioffset):(2*ip+3+ioffset+ntde))

x=cards(:,(ioffset+1):(ioffset+ip))
Doffset = 0.
if (ioffset .eq. 1) then
 Doffset = cards(:,1)
end if
doffset2 = doffset
 
if (ntde .gt. 0) then 
! do j=1,ntde
!  ftmap(j)=ioarray(4,ip+j)
!  write(6,*) ftmap(j)
!  do i=1,n
!   ft(i,j)=cards(i,(2*ip+3+ntde+j))
!   write(6,*) ft(i,j)
!  end do
  ft(:,1:ntde)=cards(:,(2*ip+3+ntde+1+ioffset):(2*ip+3+ntde*2+ioffset))
  ftmap(1:ntde)=ioarray(4,(ip+1):(ip+ntde))
!  end do
else
 ft=0
 ftmap=0
end if

!bresx=x
!ibresc=ic

do i=1,n-1,1
 if (ic(i+1)-1 .gt. -0.0001) then !if next one is event
  if (dabs(t2(i)-t2(i+1)) .lt. 0.0001) then !and if times are equal
   score_weights(i+1,:)=score_weights(i,:) ! then copy the weight down, not necessary provided the score weights are available for each of several simultaneous events
  end if
 end if
end do

! now we don't need bresx and ibresc any more

!do i=n-1,1,-1
! if (ic(i+1)-1 .gt. -0.0001) then
!  if (dabs(t2(i)-t2(i+1)) .lt. 0.0001) then
!   score_weights(i,:)=score_weights(i+1,:) ! falsch! muss runterkopiert nicht raufkopiert werden!
!   bresx(i,:)=bresx(i,:)+bresx(i+1,:)
!   ibresc(i)=ibresc(i)+ibresc(i+1)
!   ibresc(i+1)=0
!  end if
! end if
!end do

mask=.FALSE.
do j=1,Ip+ntde
 mask(j,j)=.TRUE.
end do

isep=0
XL=0.
xl0=xl-.000002

b0=0.

where(iflag .eq. 0)
 b0=dinit
elsewhere(iflag .eq. 2)
 b0=dinit
 iflag=1
endwhere

isflag=sum(iflag)
!write(6,*) "ISFLAG", isflag
b(:)=b0(:)

ITER=0
iconv=0
JCODE=0



!do while((isflag .gt. 0) .and.  ((iter .gt. 1) .and. (dabs(xl0-xl) .gt. crili)) .and. (iter .lt. imaxit))
do while((iconv .eq. 0) .and. (iter .lt. imaxit))
! write(6,*) iter, b
 iter=iter+1
 b0(:)=b(:)
 XL0=XL
 FD_sumabs0=sum(abs(fd))
 parms(10)=-10
! write(6,*) "Vor 1. LIKE", b
 if (iter .eq. 1) then
  CALL LIKE(N,IP,X,T1,t2,IC,XL,FD,vm,B,JCODE,ngv,score_weights,ntde,ft,ftmap,ilastlike,doffset,ainv, vminv)
 end if
! write(6,*) "Nach 1. LIKE"

 parms(10)=-9
 parms(8)= real (JCODE)
 IF (JCODE .GE. 1) RETURN
! IFAIL=0
! wk=-sd
! EPS=.000000000001D0
! CALL INVERT(WK,IP+ntde,IP+ntde,VM,IP+ntde,EPS,IFAIL)                             
 IF(ITER.EQ.1) then
  parms(12)=xl
  zw=0.
  zw1=matmul(vm, fd)
  zw=dot_product(fd,zw1)
  parms(7)=zw
 end if
 parms(11)=xl
 parms(10)=iter
 parms(9)=isep
 If (ISFLAG.ne.0) then
  IF(IFAIL.ne.0) then
!   "Save" Variance matrix if INVERT failed
 !  WRITE(6,*) 'Inversion failed', ITER,IFIRTH
   ICONV=0
   parms(8)=3
   return
  else  
   DO I=1,(IP+ntde)                                                      
    IF (IFLAG(I).EQ.1) then 
     TT=dot_product(vm(I,:),fd(:)*iflag(:))
     IF (DABS(TT).GT.STEP) TT=TT/DABS(TT)*STEP
     B(I)=B(I)+TT
    end if
   end do
!   half step if new log-likelihood is less than old one
   ICONV=0
   IHS=0

   CALL LIKE(N,IP,X,T1,t2,IC,XL,FD,vm,B,JCODE, ngv, score_weights,ntde,ft,ftmap,ilastlike,doffset,ainv, vminv)
   FD_sumabs=sum(abs(FD))
   do while(((FD_sumabs .gt. FD_sumabs0) .AND. (ITER.ne.1)) .AND. (ihs .le. imaxhs) .and. ((ngv .EQ. IP+ntde) .OR. (ngv .EQ. 0))) 
    IHS=IHS+1
    where (iflag .eq. 1)
     b=(b+b0)/2
    end where
    CALL LIKE(N,IP,X,T1,t2,IC,XL,FD,vm,B,JCODE, ngv, score_weights,ntde,ft,ftmap,ilastlike,doffset,ainv, vminv)
   end do
  end if
 end if
 ICONV=1
 if (isflag .gt. 0) then
  XX=dabs(B-B0)     
  XFD=dabs(fd)                                              
  IF(any(XX.GT.xconv)) ICONV=0
  if(any(xfd .gt. gconv)) iconv=0
 end if
end do




!wk=-sd
!EPS=.000000000001D0


!CALL INVERT(WK,IP+ntde,IP+ntde,VM,IP+ntde,EPS,IFAIL)       
ilastlike=1
CALL LIKE(N,IP,X,T1,t2,IC,XL,FD,vm,B,JCODE, ngv, score_weights,ntde,ft,ftmap,ilastlike,doffset,ainv, vminv)

vmls=vm
wk=vmls
do j=1,ip+ntde
  ioarray(3,j)=b(j)
  do j2=1,ip+ntde
   ioarray((3+j),j2)=vmls(j,j2)
  end do
 end do

if (irobust .eq. 1 .or. irobust .eq. 3) then   ! Lin-Wei-Varianz
 CALL dfbetaresid_lw(N,IP,X,T1,t2,IC,B,JCODE,ngv,score_weights,ntde,ft,ftmap,numbpatients,patid,ainv,dfbeta,doffset,sumabsoff)
! CALL dfbetaresid_lw(N,IP,X,T1,t2,IC,B,JCODE,ngv,score_weights,ntde,ft,ftmap,numbpatients,patid,vm,dfbeta,sumabsoff)
 vmlw=matmul(transpose(dfbeta),dfbeta)
 do j=1,ip+ntde
  do j2=1,ip+ntde
   ioarray((3+ip+ntde+j),j2)=vmlw(j,j2)
  end do
 end do
 wk=vmlw
end if
if (irobust .ge. 2) then   ! Jackknife-Varianz
 CALL dfbetaresid(N,IP,X,T1,t2,IC,B,JCODE,ngv,score_weights,ntde,ft,ftmap,numbpatients,patid,vm,dfbeta,doffset,imaxit,xconv)
 vm=real(numbpatients-1)/real(numbpatients)*matmul(transpose(dfbeta),dfbeta)
 do j=1,ip+ntde
  do j2=1,ip+ntde
   ioarray((3+2*(ip+ntde)+j),j2)=vm(j,j2)
  end do
 end do
 wk=vm
end if


CALL INVERT(WK,IP+ntde,IP+ntde,fish,IP+ntde,EPS,IFAIL)       

 

!do j=1,(ip+ntde)
! ioarray(3,j)=b(j)
! do k=1,(ip+ntde)
!  ioarray(3+j,k)=vm(j,k)
! end do
!end do
!return

!stderr=dsqrt(pack(VM,mask))
!parms(10)=-8



yy=pack(vm,mask)
yy=dabs(yy)
if (any(yy .gt. 10000)) isep=1



zw=0.
zw1=matmul(fish, b)
zw=dot_product(b,zw1)
parms(9)=zw
parms(7)=sum(abs(FD))
parms(8)=jcode
parms(11)=xl
parms(10)=iter
parms(16)=sumabsoff
cards(:,1)=doffset
DFBETA = dfbeta
!cards((1:numbpatients), (1:(ip+ntde))) = dfbeta(:)

!close(unit=6)

RETURN              

end


SUBROUTINE plusone(a)

real*8 a

a=a+1.

return
end



SUBROUTINE INVRT(A,IA)                          


!                                                                       
!...original matrix=a inverse matrix =a (on exit)                                 
!...note that a is changed on exit                                      
!                                                                       
 INTEGER IA,n
 real*8 eps                                             
 real*8, dimension (IA,ia) :: A, B, WK
 INTRINSIC DABS                                                    
                                                                       
 wk=a

 IFAIL=0
 b=a
 N=ia
 
 CALL vert(b, IA, N, WK)
 a=b
    
 RETURN
END  

SUBROUTINE INVERT(A,IA,N,B,IB,EPS,IFAIL)                          
!DEC$ ATTRIBUTES DLLEXPORT :: invert


!                                                                       
!...original matrix=a inverse matrix =b                                 
!...note that a is changed on exit                                      
!...eps is a small quantity used to see if matrix singular              
!...ifail on exit ifail=0 ok, ifail=1 matrix nearly singular            
!                                                                       
 INTEGER IA,N,ib,ifail
 real*8 eps                                             
 real*8, dimension (IA,N) :: A, B, WK
 INTRINSIC DABS                                                    
                                                                       
 wk=a

 IFAIL=0
 b=a

 CALL vert(b, IA, N, WK)

    
 RETURN
END  

function deter(ain, IA, n)

 INTEGER e, ia, n
 real*8, dimension (IA, N) :: AIN
 real*8, dimension(3+n*(n+1)) :: a
 call fact(ain, a, IA, N)
 deter=det(e,a,n)
 deter=deter*10**e

 return
end function deter

SUBROUTINE LIKE(N,IP,X,T1,t2,IC,XL,FD,VM,B,JCODE, ngv, score_weights, ntde,ft, ftmap,ilastlike,offset,ainv,vminv) 
!DEC$ ATTRIBUTES DLLEXPORT :: like

 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 real*8, dimension (IP+ntde,IP+ntde) :: DINFO, DINFOI, SDa, sdb, SDI, WK, help, vm, ainv, vminv
 real*8, dimension (IP+ntde,IP+ntde,IP+ntde) :: dabl
 real*8 SEBX, zeitp
 real*8, dimension (IP+ntde) :: XEBX
 !real*8, dimension (IP+ntde) :: XEBX, bresxges
 real*8, dimension (IP+ntde,IP+ntde) :: XXEBX
! real*8, dimension (N+1, IP, Ip, IP) :: XXXEBX
 real*8, dimension (IP+ntde) :: FD, B, h1, h2, h3
 real*8, dimension (N) :: EBX, BX, T1, t2, WKS, hh0, hh1, hh2, offset
! integer, dimension (N) :: IC,ibresc
 integer, dimension (N) :: IC
 integer, dimension (IP+NTDE) :: iflag
 !real*8, dimension (N,IP) :: X, bresx
 real*8, dimension (N,IP) :: X
 real*8, dimension (N, ip+ntde) :: xges, score_weights
 logical ifastmode
 integer ngv, ntde
 logical, dimension (N) :: maske
 real*8, dimension (N,ntde+1) :: ft
 integer, dimension (ntde+1) :: ftmap

 intrinsic dexp, dsqrt 

 dlowest=0.000000001
! write(6,*) "in LIKE"
 XL=0.

 
 ipges=ip+ntde
                                        
 ! bx=matmul(x,b)
 ! ebx=dexp(bx)


 xl=0.
 fd(:)=0.
 sda(:,:)=0.
 sdb(:,:)=0.
 sdi(:,:)=0.
 dabl(:,:,:)=0.


! Likelihood (XL) is only correct if all or none of the variables is weighted


! do i=1,n
!  do j=1,ip
!   xges(i,j)=x(i,j)
!  end do
! end do
 xges(:,1:ip)=x

if ((maxval(t1) .lt. minval(t2)) .and. (ntde .eq. 0)) then
 ifastmode=.true.
else
 ifastmode=.false.
end if


if (ifastmode .eqv. .false.) then
 do i=1,N
  if (ic(i) .ne. 0) then   
!   write(6,*) i
!   do i2=1,N
!    zeitp(i2)=t2(i)-0.00001
!   end do
   zeitp=t2(i)-0.00001
   where ((t1 .lt. zeitp) .and. (t2.ge. zeitp))
    maske=.true. 
   elsewhere
    maske=.false.
   end where
!   bresxges(1:ip)=bresx(i,1:ip)
   if (ntde .gt. 0) then
    do j=(ip+1),(ip+ntde)
!    do i2=1,n
!     xges(i2,j)=x(i2,ftmap(j-ip))*ft(i,j-ip)
!    end do
     xges(:,j)=x(:,ftmap(j-ip))*spread(ft(i,j-ip),1,n)
!     bresxges(j)=ft(i,j-ip)*bresx(i,ftmap(j-ip))
    end do
!   bresxges(ip+1:ip+ntde)=ft(i,1:ntde)*bresx(i,ftmap)
   end if
!   do i2=1,n
!    write(6,*) xges(i2,:)
!   end do
!   write(6,*) bresxges

   bx=matmul(xges,b)+offset
   ebx=dexp(bx)
   sebx=sum(ebx,1,maske)
!   write(6,*) "B",b
!   write(6,*) "BX",bx
!   write(6,*) "EBX",ebx
!   write(6,*) "SEBX", sebx
   do j=1,ipges
    hh0=xges(:,j)*ebx
    xebx(j)=sum(hh0,1,maske)
!   write(6,*) "XEBX(",j,")",xebx(j)
    do k=1,ipges
     hh1=hh0*xges(:,k)
     xxebx(j,k)=sum(hh1,1,maske)
!    write(6,*) "XXEBX(",j,",",k,")", xxebx(j,k)
    end do
   end do

   if (sebx .gt. dlowest) then
    dlogsebx=dlog(sebx)
   else
    dlogsebx=dlog(dlowest)
   endif
    

   if (ngv .eq. ipges) then 
    XL=XL+(dot_product(xges(i,:),b)-ic(i)*DLOGSEBX)*score_weights(i,1)
   else
    XL=XL+(dot_product(xges(i,:),b)-ic(i)*DLOGSEBX)
   endif
   do j=1,ipges
     FD(J)=FD(J)+(xges(i,J)-ic(i)*XEBX(J)/SEBX)*score_weights(i,j)                           
     do k=1,ipges
       SDa(J,K)=SDA(J,K)-ic(i)*((xxebx(j,k)-XEBX(J)/SEBX*XEBX(K))/SEBX)*dsqrt(score_weights(i,j))*dsqrt(score_weights(i,k))
       if (ilastlike .eq.1) then
        SDb(J,K)=SDB(J,K)-ic(i)*((xxebx(j,k)-XEBX(J)/SEBX*XEBX(K))/SEBX)*(score_weights(i,j)*score_weights(i,k))
       endif
     end do
   end do
  endif
 end do
 end if
 
 
 if (ifastmode .eqv. .true.) then
  xl=0.
  fd=0.
  sda=0.
  sdb=0.
  xebx=0.
  xxebx=0.
  sebx=0.
  ic_sum=0
  bx=matmul(xges,b)+offset
  ebx=dexp(bx)
 
  do i=N,1,-1
   if (i .gt. 1) then   !look ahead because of Breslow tie correction
    if (t2(i) .eq. t2(i-1)) then
     ic_current=0
     ic_sum = ic_sum+ic(i)
    else
     ic_current=ic_sum+ic(i)
     ic_sum=0
    end if
   else
    ic_current=ic_sum+ic(i)
    ic_sum=0
   end if
   sebx=sebx+ebx(i)
   do j=1,ipges
    hhh0=xges(i,j)*ebx(i)
    xebx(j)=xebx(j)+hhh0
    do k=1,ipges
     hhh1=hhh0*xges(i,k)
     xxebx(j,k)=xxebx(j,k)+hhh1
    end do
   end do

   if (sebx .gt. dlowest) then
    dlogsebx=dlog(sebx)
   else
    dlogsebx=dlog(dlowest)
   endif
    

   if (ngv .eq. ipges) then 
    XL=XL+(bx(i)*ic(i)-ic_current*DLOGSEBX)*score_weights(i,1)
   else
    XL=XL+(bx(i)*ic(i)-ic_current*DLOGSEBX)
   endif
   do j=1,ipges
     FD(J)=FD(J)+(xges(i,J)*ic(i)-ic_current*XEBX(J)/SEBX)*score_weights(i,j)                           
     do k=1,ipges
       SDa(J,K)=SDA(J,K)-ic_current*((xxebx(j,k)-XEBX(J)/SEBX*XEBX(K))/SEBX)*dsqrt(score_weights(i,j))*dsqrt(score_weights(i,k))
       if (ilastlike .eq.1) then
        SDb(J,K)=SDB(J,K)-ic_current*((xxebx(j,k)-XEBX(J)/SEBX*XEBX(K))/SEBX)*(score_weights(i,j)*score_weights(i,k))
       endif
     end do
   end do
  end do
 
 
 
 end if
 
 wk=-sda
 vminv=wk
 EPS=.000000000001D0
 ifail=0
 CALL INVERT(WK,ipges,Ipges,ainv,Ipges,EPS,IFAIL)
 vm=ainv
 if (ilastlike .eq. 1) then 
        jcode=ifail 
        sdb=-sdb
        sdi=vm
        wk=matmul(sdi,sdb)
        vm=matmul(wk,sdi)
 endif
! if (irobust .ne. 0) then
!  vm=sdi
! end if 
 RETURN
END 

subroutine dfbetaresid_lw(N,IP,X,T1,t2,IC,B,JCODE,ngv,score_weights,ntde,ft,ftmap,numbpatients,patid,vm,dfbeta,doffset,sumabsoff)

 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 real*8, dimension (IP+ntde,IP+ntde) :: DINFO, DINFOI, SD, SDI, WK, help, vm
 real*8, dimension (IP+ntde,IP+ntde,IP+ntde) :: dabl
 real*8 zeitp
! real*8, dimension (IP+ntde) :: XEBX, bresxges
 real*8, dimension (n) :: SEBX
 real*8, dimension (n,IP+ntde) :: XEBX, U_work
 real*8, dimension (IP+ntde,IP+ntde) :: XXEBX
! real*8, dimension (N+1, IP, Ip, IP) :: XXXEBX
 real*8, dimension (IP+ntde) :: FD, B, h1, h2, h3
 real*8, dimension (N) :: EBX, BX, T1, t2, WKS, hh0, hh1, hh2, doffset
 !integer, dimension (N) :: IC,ibresc
 !real*8, dimension (N,IP) :: X, bresx
 integer, dimension (N) :: IC
 real*8, dimension (N,IP) :: X
 real*8, dimension (N, ip+ntde) :: xges, score_weights
 integer ngv, ntde,numbpatients
 integer, dimension (N) :: patid
 logical, dimension (N) :: maske
 real*8, dimension (N,ntde+1) :: ft
 integer, dimension (ntde+1) :: ftmap
 real*8, dimension (numbpatients,ip+ntde) :: dfbeta

 intrinsic dexp, dabs 

 dlowest=0.000000001
 XL=0.

 ipges=ip+ntde

 dfbeta(:,:)=0.

 xges(:,1:ip)=x
 
 
 ! brauchen sebx(i), xebx(i,k)
sebx(:)=0.
xebx(:,:)=0.
!doffset=0.
sumabsoff=sum(dabs(doffset))
if (sumabsoff .lt. 0.0001) then
 doffset=0.
end if
!doffset=0.

do i=1, N
 if (ic(i) .ne. 0) then  
  zeitp=t2(i)-0.00001
!  where ((t1 .lt. zeitp) .and. (t2.ge. zeitp) .and. (patid .ne. ipatient))   ! corr GH 100702, ipatient gibts hier noch nicht
   where ((t1 .lt. zeitp) .and. (t2.ge. zeitp))
   maske=.true. 
  elsewhere
   maske=.false.
  end where
   if (ntde .gt. 0) then
    do j=(ip+1),(ip+ntde)
     xges(:,j)=x(:,ftmap(j-ip))*spread(ft(i,j-ip),1,n)
!     bresxges(j)=ft(i,j-ip)*bresx(i,ftmap(j-ip))
    end do
   end if

   bx=matmul(xges,b)
   bx(:)=bx(:)+doffset(:)
   ebx(:)=dexp(bx(:))

   where (maske .eqv. .false.)
    bx=0.
    ebx=0.
   end where
   sebx(i)=sum(ebx,1,maske)
   do j=1,ipges
    hh0=xges(:,j)*ebx
    xebx(i,j)=sum(hh0,1,maske)
   end do
 end if 
end do

xges(:,1:ip)=x(:,1:ip)

do i=1,n
 if (ntde .gt. 0) then
  xges(i,(ip+1):(ip+ntde))=x(i,ftmap(1:ntde))*ft(i,1:ntde)
 end if

 if (ic(i) .ne. 0) then
  u_work(i,:)=score_weights(i,:)*(xges(i,:)-xebx(i,:)/sebx(i))
 else
  u_work(i,:)=0.
 end if
 do ih=1,n
  zeitp=t2(ih)-0.00001
  if ((ic(ih) .ne. 0) .and. (t1(i) .lt. zeitp) .and. (t2(i) .ge. zeitp)) then
   if (ntde .gt. 0) then
    xges(i,(ip+1):(ip+ntde))=x(i,ftmap(1:ntde))*ft(ih,1:ntde)
   end if
   bxi=dot_product(xges(i,:),b)
   bxi=bxi+doffset(i)
   u_work(i,:)=u_work(i,:)-score_weights(ih,:)*dexp(bxi)/sebx(ih)*(xges(i,:)-xebx(ih,:)/sebx(ih)) !corr. 091008 offset(i)
!   u_work(i,:)=u_work(i,:)-score_weights(ih,:)*exp(dot_product(xges(i,:),b))/sebx(ih)*(xges(i,:)-xebx(ih,:)/sebx(ih)) !corr. 091008 offset(i)
  end if
 end do

end do

do ipatient=1,numbpatients
 where (patid .eq. ipatient)
  maske =.true.
 elsewhere
  maske =.false.
 end where
 do j=1,IPges
  dfbeta(ipatient,j)=sum(u_work(:,j),1,maske)
 end do
end do
dfbeta=matmul(dfbeta,vm)
sumabsoff=sum(dabs(doffset))
!sumabsoff=numbpatients   !checked, numbpatients is ok

RETURN
end                                                            

subroutine dfbetaresid(N,IP,X,T1,t2,IC,B,JCODE,ngv,score_weights,ntde,ft,ftmap,numbpatients,patid,ainv,dfbeta,offset, imaxit, xconv)

 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 real*8, dimension (IP+ntde,IP+ntde) :: DINFO, DINFOI, SD, SDI, WK, help, vm
 real*8, dimension (IP+ntde,IP+ntde,IP+ntde) :: dabl
 real*8 SEBX, zeitp, xconv
! real*8, dimension (IP+ntde) :: XEBX, bresxges
 real*8, dimension (IP+ntde) :: XEBX
 real*8, dimension (IP+ntde,IP+ntde) :: XXEBX, sda, ainv
! real*8, dimension (N+1, IP, Ip, IP) :: XXXEBX
 real*8, dimension (IP+ntde) :: FD, B, h1, h2, h3, bsave, step
 real*8, dimension (N) :: EBX, BX, T1, t2, WKS, hh0, hh1, hh2, offset
 !integer, dimension (N) :: IC,ibresc
 !real*8, dimension (N,IP) :: X, bresx
 integer, dimension (N) :: IC
 real*8, dimension (N,IP) :: X
 real*8, dimension (N, ip+ntde) :: xges, score_weights
 integer ngv, ntde,numbpatients
 integer, dimension (N) :: patid
 logical, dimension (N) :: maske
 real*8, dimension (N,ntde+1) :: ft
 integer, dimension (ntde+1) :: ftmap
 real*8, dimension (numbpatients,ip+ntde) :: dfbeta

 intrinsic dexp 

 dlowest=0.000000001
 XL=0.

 ipges=ip+ntde

 dfbeta(:,:)=0.

 xges(:,1:ip)=x
 bsave=b
do ipatient=1,numbpatients
 fd(:)=0.
 sda(:,:)=0.
 b=bsave 
 iconv = 0
 iter=0
 do while((iconv .eq. 0) .and. (iter .lt. imaxit))
        ! write(6,*) iter, b
      iter=iter+1
     fd(:)=0.
     SDA(:,:)=0.
     do i=1,N
      if (ic(i) .ne. 0 .and. (patid(i) .ne. ipatient)) then  
    !  if (ibresc(i) .ne. 0) then   
       zeitp=t2(i)-0.00001
       where ((t1 .lt. zeitp) .and. (t2.ge. zeitp) .and. (patid .ne. ipatient))
    !   where ((t1 .lt. zeitp) .and. (t2.ge. zeitp))
        maske=.true. 
       elsewhere
        maske=.false.
       end where
    !   bresxges(1:ip)=bresx(i,1:ip)
       if (ntde .gt. 0) then
        do j=(ip+1),(ip+ntde)
         xges(:,j)=x(:,ftmap(j-ip))*spread(ft(i,j-ip),1,n)
    !     bresxges(j)=ft(i,j-ip)*bresx(i,ftmap(j-ip))
        end do
       end if

       bx=matmul(xges,b)+offset
       ebx=dexp(bx)

       where (maske .eqv. .false.)
        bx=0.
        ebx=0.
       end where
       sebx=sum(ebx,1,maske)
       do j=1,ipges
        hh0=xges(:,j)*ebx
        xebx(j)=sum(hh0,1,maske)
        do k=1,ipges
         hh1=hh0*xges(:,k)
         xxebx(j,k)=sum(hh1,1,maske)
        end do
       end do

       if (sebx .gt. dlowest) then
        dlogsebx=dlog(sebx)
       else
        dlogsebx=dlog(dlowest)
       endif
        
       do j=1,ipges
        FD(J)=FD(J)+(Xges(i,J)-ic(i)*XEBX(J)/SEBX)*score_weights(i,j)    
         do k=1,ipges
           SDa(J,K)=SDA(J,K)-ic(i)*((xxebx(j,k)-XEBX(J)/SEBX*XEBX(K))/SEBX)*dsqrt(score_weights(i,j))*dsqrt(score_weights(i,k))
         end do
       end do

      end if
     end do
     wk=-sda
     EPS=.000000000001D0
     ifail=0
     CALL INVERT(WK,ipges,Ipges,ainv,Ipges,EPS,IFAIL)
     step=-matmul(ainv,fd)
     IF(any(abs(step) .GT. XCONV)) then
       ICONV=0
     else 
      iconv=1
     end if
     b=b-step
 end do
 do j=1,ipges
 !TT=dot_product(vm(I,:),fd(:)*iflag(:))
  dfbeta(ipatient,:)=bsave-b
!  dfbeta(ipatient,j)=fd(j)*sum(vm(j,:))
 end do
end do !ipatient
b=bsave


 RETURN


end






!      ________________________________________________________

! Code converted using TO_F90 by Alan Miller
! Date: 2006-04-13  Time: 10:09:09

!     |                                                        |
!     |  COMPUTE THE DETERMINANT OF A GENERAL FACTORED MATRIX  |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         A     -FACT'S OUTPUT                           |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         DET,E --DETERMINANT IS DET*10.**E (E INTEGER)  |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,ALOG10,DLOG10                |
!     |________________________________________________________|

FUNCTION det(e,a,nin)

INTEGER, INTENT(OUT)            :: e
integer, intent(in)             :: nin
REAL*8, INTENT(IN)              :: a(3+nin*(nin+1))
REAL*8 :: d,f,g
DOUBLE PRECISION :: c
INTEGER :: h,i,j,k,l,m,n

intrinsic DABS, dlog10

d = a(1)
IF ( DABS(d) == 1230 ) GO TO 10
!WRITE(6,*) 'ERROR: MUST FACTOR BEFORE COMPUTING DETERMINANT'
!STOP
10    e = 0
IF ( d < 0. ) GO TO 70
n = a(2)
IF ( n == 1 ) GO TO 80
d = 1.
f = 2.**64
g = 1./f
h = 64
m = n + 1
j = 0
k = 4
l = 3 - m + m*n
DO  i = k,l,m
  j = j + 1
  IF ( a(i) > j ) d = -d
  d = d*a(i+j)
  20         IF ( DABS(d) < f ) GO TO 30
  e = e + h
  d = d*g
  GO TO 20
  30         IF ( DABS(d) > g ) CYCLE
  e = e - h
  d = d*f
  GO TO 30
END DO
d = d*a(l+m)
IF ( e /= 0 ) GO TO 50
det = d
RETURN
50    IF ( d == 0. ) GO TO 90
c = DLOG10(DABS(d)) + e*DLOG10(2.d0)
e = c
c = c - e
IF ( c <= 0.d0 ) GO TO 60
c = c - 1
e = e + 1
60    f = 10.**c
IF ( d < 0. ) f = -f
det = f
RETURN
70    det = 0.
RETURN
80    det = a(5)
RETURN
90    e = 0
GO TO 70
END FUNCTION det

!

! Code converted using TO_F90 by Alan Miller
! Date: 2006-04-13  Time: 10:10:03

!      ________________________________________________________
!     |                                                        |
!     |     FACTOR A GENERAL MATRIX WITH PARTIAL PIVOTING      |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         A     --ARRAY CONTAINING MATRIX                |
!     |                 (LENGTH AT LEAST 3 + N(N+1))           |
!     |                                                        |
!     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
!     |                                                        |
!     |         N     --DIMENSION OF MATRIX STORED IN A        |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         A     --FACTORED MATRIX                        |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS                              |
!     |    PACKAGE SUBROUTINES: PACK                           |
!     |________________________________________________________|

SUBROUTINE fact(ain,a,la,n)

INTEGER, INTENT(IN)                      :: la
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(in OUT)                     :: a(3+n*(n+1))
real*8, intent(in)              ::  ain(n,n)
REAL*8 :: r,s,t
INTEGER :: e,f,g,h,i,j,k,l, m, o,p

intrinsic dabs

!IF ( la > n ) CALL pack(a,la,n)
!nicht notwendig da la immer = N

do i=1,n
 do j=1,n
  a((i-1)*n+j)=ain(i,j)
 end do
end do


r = 0.
o = n + 1
p = o + 1
l = 5 + n*p
i = -n - 3
!     ---------------------------------------------
!     |*** INSERT PIVOT ROW AND COMPUTE 1-NORM ***|
!     ---------------------------------------------
10    l = l - o
IF ( l == 4 ) GO TO 30
s = 0.
DO  k = 1,n
  j = l - k
  t = a(i+j)
  a(j) = t
  s = s + DABS(t)
END DO
IF ( r < s ) r = s
i = i + 1
GO TO 10
30    a(1) = 1230
a(2) = n
a(3) = r
i = 5 - p
k = 1
40    i = i + p
IF ( k == n ) GO TO 110
e = n - k
m = i + 1
h = i
l = i + e
!     ---------------------------------------
!     |*** FIND PIVOT AND START ROW SWAP ***|
!     ---------------------------------------
DO  j = m,l
  IF ( DABS(a(j)) > DABS(a(h)) ) h = j
END DO
g = h - i
j = i - k
a(j) = g + k
t = a(h)
a(h) = a(i)
a(i) = t
k = k + 1
IF ( t == 0. ) GO TO 100
!     -----------------------------
!     |*** COMPUTE MULTIPLIERS ***|
!     -----------------------------
DO  j = m,l
  a(j) = a(j)/t
END DO
f = i + e*o
70    j = k + l
h = j + g
t = a(h)
a(h) = a(j)
a(j) = t
l = e + j
IF ( t == 0. ) GO TO 90
h = i - j
!     ------------------------------
!     |*** ELIMINATE BY COLUMNS ***|
!     ------------------------------
m = j + 1
DO  j = m,l
  a(j) = a(j) - t*a(j+h)
END DO
90    IF ( l < f ) GO TO 70
GO TO 40
100   a(1) = -1230
GO TO 40
110   IF ( a(i) == 0. ) a(1) = -1230
RETURN
END SUBROUTINE fact


!

! Code converted using TO_F90 by Alan Miller
! Date: 2006-04-13  Time: 10:10:05

!      ________________________________________________________
!     |                                                        |
!     |   REARRANGE THE ELEMENTS OF A REAL ARRAY SO THAT THE   |
!     |  ELEMENTS OF A SQUARE MATRIX ARE STORED SEQUENTIALLY   |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         A     --REAL ARRAY CONTAINING SQUARE MATRIX    |
!     |                                                        |
!     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
!     |                                                        |
!     |         N     --DIMENSION OF MATRIX STORED IN A        |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         A     --MATRIX PACKED AT START OF ARRAY        |
!     |________________________________________________________|

SUBROUTINE packna(a,la,n)

REAL*8, INTENT(OUT)                        :: a(1)
INTEGER, INTENT(IN)                      :: la
INTEGER, INTENT(IN)                      :: n

INTEGER :: h,i,j,k,l, o

h = la - n
IF ( h == 0 ) RETURN
IF ( h > 0 ) GO TO 10
! WRITE(6,*) 'ERROR: LA ARGUMENT IN PACK MUST BE .GE. N ARGUMENT'
STOP
10    i = 0
k = 1
l = n
o = n*n
20    IF ( l == o ) RETURN
i = i + h
k = k + n
l = l + n
DO  j = k,l
  a(j) = a(i+j)
END DO
GO TO 20
END SUBROUTINE packna


!

! Code converted using TO_F90 by Alan Miller
! Date: 2006-04-13  Time: 10:09:46

!      ________________________________________________________
!     |                                                        |
!     |                INVERT A GENERAL MATRIX                 |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         V     --ARRAY CONTAINING MATRIX                |
!     |                                                        |
!     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
!     |                                                        |
!     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
!     |                                                        |
!     |         W     --INTEGER WORK ARRAY WITH AT LEAST N-1   |
!     |                      ELEMENTS                          |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         V     --INVERSE                                |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS                              |
!     |________________________________________________________|

SUBROUTINE vert(v,lv,n,w)

INTEGER, INTENT(IN OUT)                  :: lv
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN OUT)                     :: v(lv,N)
REAL*8, INTENT(OUT)                     :: w(N)
REAL*8 :: s,t
INTEGER :: i,j,k,l,m, p

!Anm GH bei den Dimensionen die Indizes verändert, waren v(lv,1) und w(1) vorher

intrinsic dabs

IF ( n == 1 ) GO TO 110
l = 0
m = 1
10    IF ( l == n ) GO TO 90
k = l
l = m
m = m + 1
!     ---------------------------------------
!     |*** FIND PIVOT AND START ROW SWAP ***|
!     ---------------------------------------
p = l
IF ( m > n ) GO TO 30
s = DABS(v(l,l))
DO  i = m,n
  t = DABS(v(i,l))
  IF ( t <= s ) CYCLE
  p = i
  s = t
END DO
w(l) = p
30    s = v(p,l)
v(p,l) = v(l,l)
IF ( s == 0. ) GO TO 120
!     -----------------------------
!     |*** COMPUTE MULTIPLIERS ***|
!     -----------------------------
v(l,l) = -1.
s = 1./s
DO  i = 1,n
  v(i,l) = -s*v(i,l)
END DO
j = l
50    j = j + 1
IF ( j > n ) j = 1
IF ( j == l ) GO TO 10
t = v(p,j)
v(p,j) = v(l,j)
v(l,j) = t
IF ( t == 0. ) GO TO 50
!     ------------------------------
!     |*** ELIMINATE BY COLUMNS ***|
!     ------------------------------
IF ( k == 0 ) GO TO 70
DO  i = 1,k
  v(i,j) = v(i,j) + t*v(i,l)
END DO
70    v(l,j) = s*t
IF ( m > n ) GO TO 50
DO  i = m,n
  v(i,j) = v(i,j) + t*v(i,l)
END DO
GO TO 50
!     -----------------------
!     |*** PIVOT COLUMNS ***|
!     -----------------------
90    l = w(k)
DO  i = 1,n
  t = v(i,l)
  v(i,l) = v(i,k)
  v(i,k) = t
END DO
k = k - 1
IF ( k > 0 ) GO TO 90
RETURN
110   IF ( v(1,1) == 0. ) GO TO 120
v(1,1) = 1./v(1,1)
RETURN
120   continue
!WRITE(6,*) 'ERROR: MATRIX HAS NO INVERSE'
STOP
END SUBROUTINE vert
