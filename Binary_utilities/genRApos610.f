C=========================================
C  Ver. 6.10
C=========================================
C For VASP
C Generates a series of 'pos_Aa' file to be used by VASP to 
C generate HF forces
C To adapt the script 'run' to your machine write a proper name of
C setting card at the first line of 'run'
C----------------------------------------
      program GenPos
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(lmax=1000)
      real*8 APRIM(3,3),X(3)
      real*8 Rnew(:,:),Rold(:,:),u(:,:)
      integer LI(lmax),LIL(lmax),NAT(lmax),m1,m4,idirect,i
      character*60 text1,text2,text3,text0,text,hftext
      character*60 Line1,Line2
      character*4 txt
      character*4 ch,name(:),chh,chr(lmax),ch4,chh4
      character pos*9,out*9,alf*26,run*12,ttx*40
      character config*12,hfposcar*8
      logical existence
      integer iatom(:)
      allocatable :: Rnew,Rold,u,name,iatom
      data alf/'abcdefghijklmnoprqstuvwxyz'/
C----------------------------------------
C Displayed on screen
      write(*,*)' Program <genRApos610>'
      write(*,*)' generates all POSCARs of VASP'
      write(*,*)' needed to find dielectric tensors.'
      write(*,*)' Necessary displacements must be given in RAMDISCAR'
      write(*,*)' which is <project>.e44 of Phonon'
      write(*,*)' Names of new POSCAR files with atoms displaced'
      write(*,*)' are generated automatically.'
      write(*,*)' This program reads CONTCAR (not POSCAR).'
      write(*,*)' You need'
      write(*,*)'     CONTCAR'
      write(*,*)'     RAMDISCAR'
C      write(*,*)'     <configure>'
      write(*,*)' files in current directory.'
      write(*,*)' Program <genRApos610> copies CONTCAR to RAMPOSCAR.'
      write(*,*)' It also generates <runRA> script, to calculate'
      write(*,*)' dielectric tensors.' 
      write(*,*)' Use: <chmod a+x runRA> before running this script.'
      pause ' Press any key, or type go, to continue'
C----------------------------------------
C read CONTCAR, and eliminate line with atomic symbols, if needed
      open(10,file='CONTCAR',status='OLD',ERR=30)
      idirect=0
      do i=1,7
        read(10,'(A)')txt
      enddo
      if(txt(1:1).eq.'D'.OR.txt(1:1).eq.'d')then
        idirect=7
      else
        read(10,'(A)')txt
        if(txt(1:1).eq.'D'.OR.txt(1:1).eq.'d')then
          idirect=8
        else
          goto 32
        endif
      endif
      rewind(10)
C-----
      inxx=0
      KLI=0
      goto 20
  22  rewind(10)
  20  read(10,'(A)')text1
      read(10,'(A)')text2
      do J=1,3
        read(10,*)(APRIM(J,I),I=1,3)
      enddo
      if(idirect.eq.8) read(10,'(A)')txt
C-----
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
      allocate(Rnew(NMAX,3),Rold(NMAX,3))
      do I=1,NMAX
        read(10,*)(Rold(I,J),J=1,3)
      enddo
      close(10)
C-------------------------------------------------------
C Copy CONTCAR to RAMPOSCAR
      open(10,file='CONTCAR',ERR=30)
      open(11,file='RAMPOSCAR',ERR=30)
       read(10,'(A)')hftext
       write(11,'(A)')hftext
       read(10,'(A)')hftext
       write(11,'(A)')hftext
       do J=1,3
         read(10,*)(X(I),I=1,3)
         write(11,1000)(X(I),I=1,3)
       enddo
       if(idirect.eq.8) read(10,'(A)')txt
       read(10,*)(LIL(II),II=1,LICZBA)
       write(11,1001)(LIL(II),II=1,LICZBA)
       read(10,'(A)')hftext
       write(11,'(A)')hftext
       do I=1,NMAX
         read(10,*)(X(J),J=1,3)
         write(11,1002)(X(J),J=1,3)
       enddo
       close(10)
       close(11)
C-------------------------------------------------------
C Open <configure> file
C      config='configure'
C      open(13,file=config,status='OLD',ERR=33)
C     do i=1,5
C        read(13,'(A)',END=33,ERR=33)txt
C      enddo
C      read(13,'(A)')Line1
       Line1='#!/bin/bash '
C      read(13,'(A)')Line2
C      close(13)
C-------------------------------------------------------
C read RAMDISCAR from PHONON ver.6.10
C write to script "kopia"
      open(11,file='RAMDISCAR',status='OLD',ERR=31)
      open(14,file='kopia')
      write(14,'(A)')Line1
      write(14,'(A)')'#'
      write(14,'(A)')'mkdir AXML'
      write(14,'(A)')'#'
      read(11,'(A)')text
      read(11,*)KI
      allocate(name(KI),iatom(KI),u(KI,3))
      do J=1,KI
        read(11,*,END=31,ERR=31)name(J),iatom(J),(u(J,I),I=1,3)
      enddo
      close(11)
C-------------------------------------------------------
C Write series of POSCAR: pos_aa, pos_ab, ... with diplacements
      I=0
      NN=0
      chh=name(1)
      do J=1,KI
       if(name(J).eq.chh)then
         NN=NN+1
       else
         NN=1
       endif
       I=I+1
       chh=name(I)
       ch4=name(J)
       m4=LEN_TRIM(ch4)
       pos='pos_'//ch4(1:m4)//alf(NN:NN)
C--------
       chh4=ch4(1:m4)//alf(NN:NN)
       m5=LEN_TRIM(chh4)
       ttx='cp ra_pos_'//
     *      chh4(1:m5)//'/vasprun.xml AXML/'//chh4(1:m5)//'.xml'
       write(14,'(A)')ttx
C------
       open(10,file=pos)
        write(10,'(A)')text1
        write(10,'(A)')text2
        do M=1,3
          write(10,1000)(APRIM(M,II),II=1,3)
        enddo
        write(10,1001)(LI(II),II=1,LICZBA)
        write(10,'(A)')text3
        do K=1,NMAX
          do LLL=1,3
            Rnew(K,LLL)=Rold(K,LLL)
          enddo
          if(K.eq.iatom(J))then
            do L=1,3
              Rnew(K,L)=Rold(K,L)+u(J,L)
            enddo
          endif
          write(10,1002)(Rnew(K,L),L=1,3)
        enddo
       close(10)
      enddo
      close(14)
      deallocate(name,iatom,u)
      deallocate(Rnew,Rold)
 1000 format(1X,3F22.16)
 1001 format(999I4)
 1002 format(3F20.16)
C-------------------------------------------------------
C-------------------------------------------------------
C Write script <runRA> 
      run='runRA'
      open(12,file=run,ERR=30)
      write(12,'(A)')Line1
      write(12,'(A)')'#'
!       write(12,'(A)')' VASP='//Line2
!       write(12,'(A)')'#'
      write(12,'(A)')'for d in pos_*; do'
      write(12,'(A)')'('
      write(12,'(A)')'   mkdir ra_$d'
      write(12,'(A)')'   cd ra_$d'
      write(12,'(A)')'   cp ../$d POSCAR'
      write(12,'(A)')'   ln -s ../KPOINTS .'
      write(12,'(A)')'   ln -s ../POTCAR .'
      write(12,'(A)')'   ln -s ../INCAR .'
      inquire(file="vdw_kernel.bindat",exist=existence)
      if(existence) then
        write(12,'(A)')'   ln -s ../vdw_kernel.bindat .'
      end if
!       write(12,'(A)')'   $VASP'
      write(12,'(A)')')'
      write(12,'(A)')'done'
      close(12)
C-------------------------------------------------------
C-------------------------------------------------------
      goto 9
   30 write(*,*)'***** CONTCAR not found *****'
      goto 9
   31 write(*,*)'***** RAMDISCAR not found *****'
      goto 9
   32 write(*,*)'***** Wrong format in CONTCAR *****'
      goto 9
   33 write(*,*)'***** Wrong format in configure *****'
    9 end  
C=========================================
