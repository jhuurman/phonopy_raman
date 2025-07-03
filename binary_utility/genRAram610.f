C==============================
C  Ver. 6.10
C==============================
C For VASP
C Generates RAMFILE from OUTCARs in folders /opt_pos_Aa' etc.,
C and from files RAMDISCAR and RAMPOSCAR.
C----------------------------------------
      program GenRam
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 APRIM(3,3),eps(3),x(3)
      real*8 Rold(:,:),u(:,:),uu(3),vv(3)
      integer LI(1000),LIL(1000),NAT(1000),m4,IDIR
      character*12 CONTCAR, DISPCAR
      character*60 text1,text2,text3,text,text0,FirstLine,hftext
      character*74 info(4)
      character*4 ch,name(:),chh,ch1(:),ch4
      character*4 names(:),na(:),NameAtom(:)
      character*5 outname(:)
      character out*24,alf*26,vasp*12,run*12
      integer i,k,JJ,KSTEP(:)
      integer iatom(:),na1(:),na2(:)
      real*8 RR(:,:),DISP(:,:),EPSIL(:,:,:),EPSIL0(3,3)
      allocatable :: Rold,u,name,iatom,names,na,na1,na2,ch1
      allocatable :: RR,DISP,EPSIL,outname,KSTEP,NameAtom
      data alf/'abcdefghijklmnoprqstuvwxyz'/
C----------------------------------------
C Display on sreen
      write(*,*)' Program <genram610.f>'
      write(*,*)' generates Raman files for PHONON from '
      write(*,*)' rather vasprun.xml then OUTCAR files in'
      write(*,*)' folders /ra_pos_aa, /ra_pos_ab, etc.'
      write(*,*)' Current folder must contain:'
      write(*,*)' - file RAMDISCAR. Displacements along three'
      write(*,*)'        X,Y,Z directions'
      write(*,*)' - file RAMPOSCAR. Positions of optimized supercell.'
      write(*,*)' '
      write(*,*)' Using script "kopia" files vasprun.xml from folders'
      write(*,*)'        /ra_pos_aa, /ra_pos_Aa, etc,'
      write(*,*)'        should be copied to AXML folder'
      write(*,*)'        as files: aa.xml, ab.xml, etc.'
      write(*,*)' '
      write(*,*)' Output files are:'
      write(*,*)' - RAMAN.out, dielectric constant tensors,'
      write(*,*)' - RAMFILE, containing Raman tensors.'
      write(*,*)'        It is used in PHONON for Raman scattering'
      write(*,*)' '
      write(*,*)' Remark: Resulting Raman tensors are:'
      write(*,*)' - multiplied by volume of primitive unit cell,'
      write(*,*)' - devided by 4pi,'
      write(*,*)' '
      pause ' Press any key, or type go, to continue'
C------------------------------------------------------------------
C Four different methods of extracting Raman tensors.
C Below only method 4 is in use.
C      write(*,*)' Define import data for Raman tensor'//
C     *                             ' (recommended ** 4 **):'
C      info(1)='  1. Microscopic dielectric constants from OUTCARs'// 
C     *                                        ' (less digits)(DIEL)'
C      info(2)='  2. Microscopic dielectric constants from'//
C     *                                        ' vasprun.xml (AXML)'
C      info(3)='  3. Macroscopic dielectric constants from'//
C     *                              ' OUTCARs (less digits) (DIEL)'
C      info(4)='  4. Macroscopic dielectric constants from'//
C     *                                        ' vasprun.xml (AXML)'
C      write(*,999)info(1)
C      write(*,999)info(2)
C      write(*,999)info(3)
C      write(*,999)info(4)
C      write(*,*)' '
C  40  write(*,*)' Type MIC-DIEL(1), MIC-AXML(2), MAC-DIEL(3),'//
C     *                            ' MAC-AXML(4), Cancel (0):'
C      read(*,*)IDIR
C      write(*,2040)IDIR
C      if(IDIR.lt.0.or.IDIR.gt.4) goto 40
C      if(IDIR.eq.0) goto 10
C 2040 format(1X,'You typed: ',1I2)
C-----------------------------------------------------------------
C read RAMPOSCAR
      IDIR=4  ! fixed selection index
      open(10,file='RAMPOSCAR',status='OLD',ERR=30)
      inxx=0
      KLI=0
      goto 20
  22  rewind(10)
  20  read(10,'(A)')text1
      read(10,'(A)')text2
      do J=1,3
!old        read(10,*)(APRIM(J,I),I=1,3)
        read(10,*)(APRIM(I,J),I=1,3)
      enddo
      KLI=KLI+1
      read(10,*,END=21,ERR=21)(LI(I),I=1,KLI)
      if(inxx.eq.1) goto 23
      goto 22
   21 KLI=KLI-2
      inxx=1
      goto 22
   23 LICZBA=KLI
       read(10,'(A)')text3
      NMAX=0
      do I=1,LICZBA
        NMAX=NMAX+LI(I)
      enddo
      allocate(Rold(NMAX,3))
      do I=1,NMAX
        read(10,*)(Rold(I,J),J=1,3)
      enddo
      close(10)
C-------------------------------------------------------
C read RAMDISCAR from PHONON ver.5.11
      open(11,file='RAMDISCAR',status='OLD',ERR=31)
      read(11,'(A)')text
      read(11,*)KI
      allocate(name(KI),iatom(KI),u(KI,3))
      do J=1,KI
        read(11,*,END=31,ERR=31)name(J),iatom(J),(u(J,I),I=1,3)
      enddo
C-------------------------------------------------------
      allocate(outname(KI),RR(3,KI),DISP(3,KI),EPSIL(3,3,KI))
!--------------------
      read(11,*,END=31,ERR=31)natosc,kki
      allocate(KSTEP(kki),NameAtom(kki))
      i=0
      k=1
      iniatom=iatom(1)
      do J=1,KI
        if(iniatom.eq.iatom(J))then
           i=i+1
        else
           k=k+1
           i=1
        endif
        KSTEP(k)=i
        iniatom=iatom(J)
      enddo
      if(k.ne.kki) write(*,*)' Scatter of degree of freedom is incorect'
      if(NMAX.ne.natosc)then
         write(*,*)' Wrong number of atoms in unit cell'
      endif
      allocate(ch1(kki),na1(kki),na2(kki))
      do i=1,kki
        read(11,*,END=31,ERR=31)ch1(i),na1(i)
      enddo
        do i=1,kki-1
           na2(i)=na1(i+1)-1
        enddo
      na2(kki)=natosc
      allocate(names(natosc))
      do k=1,kki
        NameAtom(k)=ch1(k)
        do i=na1(k),na2(k)
           names(i)=ch1(k)
        enddo
      enddo
      deallocate(ch1,na1,na2)
      close(11)
C-------------------------------------------------------
C     RAMAN.out
C-------------------------------------------------------
C Name of folders where are OUTCARs
      open(10,file='RAMAN.out',ERR=9)
C      write(10,999)info(IDIR)
      write(10,'(A)')' '//text
      write(10,*)' '
      write(10,'(A)')
     *' Dielectric constant tensors from the following files:'
C      write(10,'(A)')' - zero atomic displacments "zero.xml"'
      write(10,'(A)')'  - one displaced atom (from position X.Y,Z)'//
     *                         ' along x, y, or z direction,'
      write(10,'(A)')'     with amplitude in angstroms.'
      write(10,'(A)')
     *' Atom position in fractional coordinates with respect to metrix:'
! Write metrix
      do J=1,3
        write(10,1001)(APRIM(J,I),I=1,3)
      enddo
C-------------------------------------------------------
! Avoid Zero dielectric constant
! It is not needed for the linear term of dielectric
!      constant expansion with atomic the atomic displacements.
!      It is needed for the second order term
!       goto 89  
!       iuni=14
!       open(iuni,status='scratch') 
!       select case(IDIR)
!         case(1,3)
!           write(10,'(A)')' zero.opt'
!           out='DIEL/zero.opt'                 ! for "genram610.f", DIEL
!         case(2,4)
!           write(10,'(A)')' zero.xml'
!           out='AXML/zero.xml'                 ! for "genram610.f", AXML
!       end select
!       call readDiel(iuni,out,nHF,inx,IDIR)
!       rewind(iuni)
!       do I=1,3
!         read(iuni,*)(eps(JJ),JJ=1,3)
!         write(10,1002)(eps(JJ),JJ=1,3)
!         do JJ=1,3
!           EPSIL0(I,JJ)=eps(JJ)
!         enddo
!       enddo
!       close(iuni)
!   89  continue
C-------------------------------------------------------
      I1=0
      NN=0
      chh=name(1)
      do J=1,KI
       if(name(J).eq.chh)then
         NN=NN+1
       else
         NN=1
       endif
       I1=I1+1
       chh=name(I1)
       ch4=name(J)
       m4=LEN_TRIM(ch4)
       select case(IDIR)
         case(1,3)
           out='DIEL/'//ch4(1:m4)//alf(NN:NN)//'.opt'     ! for "genram610.f", DIEL
         case(2,4)
           out='AXML/'//ch4(1:m4)//alf(NN:NN)//'.xml'     ! for "genran610.f", AXML
       end select
       outname(J)=ch4(1:m4)//alf(NN:NN)
C------------------------------------------
C Write displacement
       do JJ=1,3
         uu(JJ)=u(J,JJ)
       enddo
       call MULMV(3,APRIM,uu,vv)
       write(10,1003)outname(J),
     *               (Rold(iatom(J),JJ),JJ=1,3),(vv(JJ),JJ=1,3)
        do ji=1,3
           RR(ji,J)=Rold(iatom(J),ji)
        enddo
! displacements in angstroms
        do ji=1,3
           DISP(ji,J)=vv(ji)
        enddo
C------------------------------------------
C Read HF forces from out file
       iuni=14
       open(iuni,status='scratch') 
       call readDiel(iuni,out,nHF,inx,IDIR)
       if(inx.eq.1) goto 9
C-----
       rewind(iuni)
       do I=1,3
         read(iuni,*)(eps(JJ),JJ=1,3)
         write(10,1002)(eps(JJ),JJ=1,3)
         do JJ=1,3
           EPSIL(I,JJ,J)=eps(JJ)
         enddo
       enddo
       close(iuni)
      enddo
      deallocate(names)
      deallocate(name,iatom,u)
      deallocate(Rold)
      close(10)
C-----------------------------------------------------------------------
C Writes Raman matrices 3x3, and composses Raman Tensor 3x9
      call Raman_Tensor(APRIM,KI,outname,KSTEP,NameAtom,kki,RR,
     *                           DISP,EPSIL,EPSIL0,IDIR,info,text)
C------------------------------------------------------------------------
      deallocate(outname,RR,DISP,EPSIL)
      deallocate(KSTEP,NameAtom)
C-------------------------------------------------------
  999 format(A)
 1000 format(1X,3F20.16)
 1001 format(1X,3F22.16)
 1002 format(1X,3F20.7)
 1003 format(1X,1A5,2X,3F10.6,2X,3F10.6)
      goto 10
   30 write(*,*)'***** RAMPOSCAR not found *****'
      goto 9
   31 write(*,*)'***** RAMDISCAR not found *****'
    9 pause ' Press any key to continue'
   10 END  
C===========================================================
C===========================================================
C                   SUBROUTINES
C===========================================================
C===========================================================
      SUBROUTINE readDiel(iuni,out,nHF,inx,IDIR)
C Subroutine cuts out HF forces from out=/ph_pos_Aa/OUTCAR
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 eps(3)
      integer i,j,nHF,IDIR
      character*24  out
      character*132 str
      character*6   strb
      nHF=3
      inx=0
C---------------------------------------------------------
      select case(IDIR)
C-------------------
! File OUTCAR, MICOROSCOPIC DIELECTRIC CONSTANT
       case(1)
        open(15,file=out,status='old',err=30)
!       Find microscopic dielelectric constant
  31    read(15,1001,END=11)str
        if(str(2:37).NE.'MICROSCOPIC STATIC DIELECTRIC TENSOR')
     *               GOTO 31
        read(15,1001,END=11)str
        do j=1,3
         read(15,*,END=11,ERR=11)(eps(i),i=1,3)
         write(iuni,1000)eps
        enddo
   11   close(15)
C-------------------
! File "vasprun.xml"
       case(2)
        open(15,file=out,status='old',err=30)
!       Find microscopic dielectric constant.
        read(15,1001,END=21)str
        do 
          if(str(3:30).EQ.'<varray name="epsilon_rpa" >') exit
          if(str(3:30).EQ.'<varray name="epsilon_scf" >') exit
        end do
        do j=1,3
         read(15,*,END=21,ERR=21)strb,(eps(i),i=1,3)
         write(iuni,*)eps
        enddo
   21   close(15)
C-------------------
! File OUTCAR, MACOROSCOPIC DIELECTRIC CONSTANT
       case(3)
        open(15,file=out,status='old',err=30)
!       Find macroscopic dielectric constant
  61    read(15,1001,END=51)str
        if(str(2:37).NE.'MACROSCOPIC STATIC DIELECTRIC TENSOR')
     *               GOTO 61
!       Find second dielectric constant
  62    read(15,1001,END=51)str
        if(str(2:37).NE.'MACROSCOPIC STATIC DIELECTRIC TENSOR')
     *               GOTO 62
        read(15,1001,END=51)str
        do j=1,3
         read(15,*,END=51,ERR=51)(eps(i),i=1,3)
         write(iuni,1000)eps
        enddo
   51   close(15)
C-------------------
! File "vasprun.xml", MACOROSCOPIC DIELECTRIC CONSTANT
       case(4)
        open(15,file=out,status='old',err=30)
!       Find first dielectric constant
        do 
          read(15,1001,END=71) str
          if(str(3:26).EQ.'<varray name="epsilon" >') exit
          if(str(3:30).EQ.'<varray name="epsilon_scf" >') exit
        end do
        do j=1,3
         read(15,*,END=71,ERR=71)strb,(eps(i),i=1,3)
         write(iuni,*)eps
        enddo
   71   close(15)
C--------------------
      end select
C---------------------------------------------------------
 1000 format(6X,3F22.7)
 1001 format(A)
 1019 format(1A6,3F20.7)
 1020 format(6X,3F20.7)
      goto 9
   30 write(*,*)' File: '//out//' not found'
      close(15)
      inx=1
   9  END
C===========================================================
C===========================================================
      SUBROUTINE Raman_Tensor(APRIM,KI,outname,KSTEP,NameAtom,kki,RR,
     *                     DISP,EPSIL,EPSIL0,IDIR,info,text)
C-----------------------------------------------------------
      real*8 APRIM(3,3),pi,pi4,VPRIM
      character*5 outname(KI)
      character*4 NameAtom(kki)
      character*60 text
      character*74 info(4)
      integer i,II,j,k,l,KK,n,m,i1,i2,kkk,KSTEP(kki),IDIR
      integer low(:),high(:)
      real*8 RR(3,KI),DISP(3,KI),EPSIL(3,3,KI),EPSIL0(3,3)
!     real*8 EPSI0(3,3)
      real*8 EPSI(:,:,:),RRR(:,:),RAM(:,:,:)
      real*8 A(3,3),AI(3,3),B(3),BB(3,3),C(3)  ! for matrix equation A.B=C, B - unknown
      allocatable :: EPSI,RRR,RAM,low,high
C-----------------------------------------------------------
      CALL DETER3(APRIM,VPRIM)
      allocate(EPSI(3,3,KI))
C Deviding dielectric constant by 4pi:
      pi=2.0D00*dasin(1.0D00)
      pi4=4.0D00*pi
      do k=1,KI
      do i=1,3
      do j=1,3
        EPSI(i,j,k)=VPRIM*EPSIL(i,j,k)/pi4
      enddo
      enddo
      enddo
!      do i=1,3
!      do j=1,3
!        EPSI0(i,j)=EPSIL0(i,j)/pi4
!      enddo
!      enddo
C-----------------------
! Open RAMFILE file
      open(10,file='RAMFILE')
!      write(10,1000)info(IDIR)
      write(10,'(A)')'! '//text
      write(10,207)
      write(10,199)
      do J=1,3
        write(10,1001)(APRIM(J,I),I=1,3)
      enddo
 1000 format('!',A)
 1001 format('!',3F22.16)
C------------------------
C Legend
      write(10,200)
      write(10,201)
      write(10,202)
      write(10,203)
      write(10,204)
      write(10,205)
      write(10,206)
      write(10,207)
      write(10,208)
      write(10,209)
      write(10,210)
      write(10,211)
      write(10,212)

  199 format('! Unit cell metrix:')
  200 format('! Assumption:')
  201 format('! Electric polarizibility tensor H = (eps - 1)/4pi')
  202 format('! is linearly dependent on atomic displacement.')
  203 format('! The H is represented by H = B.U')
  204 format('!     where U is atom displacement along X,Y, or Z,')
  205 format('!     B=[H(+U)-H(-U)]/[2*U*4*pi] defines')
  206 format('!     1-st order Raman scattering (in [1/A]).')
  207 format('!')
  208 format('! H is polarizibility tensor:')
  209 format('!                            H11  H12  H13')
  210 format('!                            H21  H22  H23')
  211 format('!                            H31  H32  H33')
  212 format('!')
C------------------------
C For a given atom deal with all displacements from "low" till "high"
      allocate(low(kki),high(kki))
      low(1)=1
      high(1)=KSTEP(1)
      do i=2,kki
        low(i) =low(i-1) +KSTEP(i-1)
        high(i)=high(i-1)+KSTEP(i)
      enddo
C---------------------------------------------------------------------
C For Raman tensor RAM=0.0
      allocate(RRR(3,kki),RAM(3,3*3,kki))
      do KK=1,kki
      do m=1,3*3
      do n=1,3
         RAM(n,m,KK)=0.0D00
      enddo
      enddo
      enddo
C---------------------------------------------------------------------
C Writes 6 displacements for one independent atom in a loop KK
      do KK=1,kki
        write(10,220)
        write(10,221)
        write(10,222)(outname(j),j=low(KK),high(KK))
        write(10,223)(RR(j,low(KK)),j=1,3)
        write(10,207)
       if(low(KK)  .le.high(KK))write(10,224)(DISP(j,low(KK)  ),j=1,3)
       if(low(KK)+1.le.high(KK))write(10,227)(DISP(j,low(KK)+1),j=1,3)
       if(low(KK)+2.le.high(KK))write(10,225)(DISP(j,low(KK)+2),j=1,3)
       if(low(KK)+3.le.high(KK))write(10,227)(DISP(j,low(KK)+3),j=1,3)
       if(low(KK)+4.le.high(KK))write(10,226)(DISP(j,low(KK)+4),j=1,3)
       if(low(KK)+5.le.high(KK))write(10,227)(DISP(j,low(KK)+5),j=1,3)

  220 format('!-------------------------------------------------------')
  221 format('! Partial Raman Tensors:')
  222 format('! Atoms: ',48A5)
  223 format('! Atomic Position (X,Y,Z):',3F10.6)
  224 format('! Displacement X (in A):   ',3F10.6)
  225 format('! Displacement Y (in A):   ',3F10.6)
  226 format('! Displacement Z (in A):   ',3F10.6)
  227 format('!                          ',3F10.6)
C-----------------------------------
C Copy displacements
        do j=1,3
          RRR(j,KK)=RR(j,low(KK))
        enddo
C-----------------------------------
! Raman tensor
       do II=low(KK),high(KK)-1,2     !  loop over II
          kr=(II-low(KK))/2+1
C---------------------
          do n=1,3
          do m=1,3
            BB(n,m)=(EPSI(n,m,II)-EPSI(n,m,II+1))/
     *                  (DISP(kr,II)-DISP(kr,II+1))
          enddo
          enddo
C---------------------
      select case(kr)
        case(1) 
          write(10,228)(BB(1,m),m=1,3)
        case(2) 
          write(10,229)(BB(1,m),m=1,3)
        case(3) 
          write(10,230)(BB(1,m),m=1,3)
      end select
          write(10,231)(BB(2,m),m=1,3)
          write(10,231)(BB(3,m),m=1,3)
          write(10,207)
  228 format('! Raman Tensor X: ',3F20.7,3X,3F20.7,3X,3F20.7)
  229 format('! Raman Tensor Y: ',3F20.7,3X,3F20.7,3X,3F20.7)
  230 format('! Raman Tensor Z: ',3F20.7,3X,3F20.7,3X,3F20.7)
  231 format('!                 ',3F20.7)
          do m=1,3
            RAM(1,3*(kr-1)+m,KK)=BB(1,m)
            RAM(2,3*(kr-1)+m,KK)=BB(2,m)
            RAM(3,3*(kr-1)+m,KK)=BB(3,m)
          enddo
       enddo    ! loop II
      enddo     ! loop KK
C------------------------------------
  237 format('!======================================================')
      write(10,237) 
      write(10,235)
      write(10,236)
      write(10,220)
      write(10,232)
  235 format('! RAMAN TENSORS')
  236 format('! The Raman tensors are multiplied by volume of ',
     *       'primitive unit cell Vprim, and devided by 4pi.')
  232 format('! Name of Atom with Displacements, Atomic ',
     *           'Position (X,Y,Z), and Raman tensor:')
  233 format(6X,1A4,1X,3F10.6)
  234 format(1X,3F16.4,1X,3F16.4,1X,3F16.4)
      do KK=1,kki
        write(10,233)NameAtom(KK),(RRR(j,KK),j=1,3)
        khl=3
        do m=1,3
           write(10,234)(RAM(m,j,KK),j=1,3*khl)
        enddo
      enddo      
      deallocate(low,high)
      deallocate(RRR,RAM)
      close(10)
      deallocate(EPSI)
      END
C===========================================================
C===========================================================
C     ALGEBRAIC SUBROUTINES
C===========================================================
C===========================================================
      SUBROUTINE INVER3(A,B)
C-----------------------------------------------------------
C  Calculates from A(3x3) matrix its invers B(3x3) and AB=1
C-----------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(3,3),B(3,3),EPS
      EPS=1.0D-12
      CALL DETER3(A,V)
      IF((V+EPS).EQ.0)WRITE(*,*)'Determinant equal zero'
      B(1,1)= A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(2,1)=-A(2,1)*A(3,3)+A(3,1)*A(2,3)
      B(3,1)= A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(1,2)=-A(1,2)*A(3,3)+A(3,2)*A(1,3)
      B(2,2)= A(1,1)*A(3,3)-A(3,1)*A(1,3)
      B(3,2)=-A(1,1)*A(3,2)+A(3,1)*A(1,2)
      B(1,3)= A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(2,3)=-A(1,1)*A(2,3)+A(2,1)*A(1,3)
      B(3,3)= A(1,1)*A(2,2)-A(2,1)*A(1,2)
      DO 100 J=1,3
      DO 100 I=1,3
  100 B(I,J)=(B(I,J)+EPS)/(V+EPS)
      END
C===========================================================
C===========================================================
      SUBROUTINE MULMV(N,A,V,U)
C-----------------------------------------------------------
C  Multiplies matrices A.v=u. A is (NxN) matrix.
C  u,v are N column vectors.
C-----------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(N,N),U(N),V(N)
      DO 100 I=1,N
      U(I)=0.0D00
      DO 100 J=1,N
  100 U(I)=U(I)+A(I,J)*V(J)
      END
C===========================================================
C===========================================================
      SUBROUTINE MULMV3(A,V,U)
C-----------------------------------------------------------
C  Multiplies matrices A.v=u. A is (3x3) matrix.
C  u,v are 3 column vectors.
C-----------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(3,3),U(3),V(3)
      DO 100 I=1,3
      U(I)=0.0D00
      DO 100 J=1,3
  100 U(I)=U(I)+A(I,J)*V(J)
      END
C===========================================================
C===========================================================
      SUBROUTINE DETER3(A,V)
C-----------------------------------------------------------
C  Calculates determinant of matrix A of (3x3)
C-----------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(3,3)
      V=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+
     *       A(1,3)*A(2,1)*A(3,2)-A(3,1)*A(2,2)*A(1,3)-
     *       A(3,2)*A(2,3)*A(1,1)-A(3,3)*A(2,1)*A(1,2)
      END
C===========================================================
c********************** END OF PROGRAM**********************
C===========================================================
C===========================================================
