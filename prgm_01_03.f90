      Program prgm_01_01
!
!     This program reads a 3x3 matrix from a user-provided input file. After the
!     file is opened and read, it is closed and then printed.
!
!     H. P. Hratchian, 2019.
!
      implicit none
      integer,parameter::inFileUnitA=10,inFileUnitB=20
      integer::errorFlag,i
      real,dimension(3,3)::matrixInA,matrixInB,matrixOutAB
      character(len=128)::fileNameA,fileNameB
!
!
!     Start by asking the user for the name of the data file.
!
      write(*,*)' What is the name of the input data file?'
      read(*,*) fileNameA
!
!     Open the data file and read matrixInA from that file.
!
      open(unit=inFileUnitA,file=TRIM(fileNameA),status='old',  &
        iostat=errorFlag)
      if(errorFlag.ne.0) then
        write(*,*)' There was a problem opening the input file.'
!        goto 999
        stop
      endIf
      do i = 1,3
        read(inFileUnitA,*) matrixInA(1,i),matrixInA(2,i),matrixInA(3,i)
      endDo
      close(inFileUnitA)
!
!     Asking for the user for the name of the second data file.
!
      write(*,*)' What is the name of the second input data file?'
      read(*,*) fileNameB
!
!     Open the data file and read matrixInA from that file.
!
      open(unit=inFileUnitB,file=TRIM(fileNameB),status='old',  &
        iostat=errorFlag)
      if(errorFlag.ne.0) then
        write(*,*)' There was a problem opening the input file.'
!        goto 999
        stop
      endIf
      do i = 1,3
        read(inFileUnitB,*) matrixInB(1,i),matrixInB(2,i),matrixInB(3,i)
      endDo
      close(inFileUnitB)
!
!     Call the subroutine PrintMatrix to print matrixInA.
!
      call PrintMatrix3x3(matrixInA)
      call PrintMatrix3x3(matrixInB)
!
      call MatrixProduct3x3(matrixInA,matrixInB,matrixOutAB)
!
!      write(*,*)' The product of the two matrices is... '
!
      call PrintMatrix3x3(matrixOutAB)
!
!  999 continue
      End Program prgm_01_01


      Subroutine PrintMatrix3x3(matrix)
!
!     This subroutine prints a 3x3 real matrix. The output is written to StdOut.
!
      implicit none
      real,dimension(3,3),intent(in)::matrix
      integer::i
!
!     Format statements.
!
 1000 format(3(2x,f5.1))
!
!     Do the printing job.
!
      write(*,*)' Printing Matrix'
!
! ADDED CODE HERE
!
!     Karnamohit, 01/25/2019.
!
      do i = 1,3
        write(*,1000)matrix(i,:)
      end do
!
      return
      End Subroutine PrintMatrix3x3


      Subroutine MatrixProduct3x3(matrix1,matrix2,matrix12)
!
!     This subroutine multiplies two 3x3 matrices using the MatMul intrinsic
!     function.
!
      implicit none
      integer::i,j,k,l
      real,dimension(3,3),intent(in)::matrix1,matrix2
      real,dimension(3,3),intent(out)::matrix12
      real,parameter,dimension(2)::dim1=shape(matrix1)
      real,parameter,dimension(2)::dim2=shape(matrix2)
!
!      i = dim1(1)
!      j = dim1(2)
!      k = dim2(1)
!      l = dim2(2)
!
!      if (i.ne.k .OR. j.ne.l) then
!        write(*,*)' The dimensions of the matrices do not conform.'
!        stop
!      else
!
!
       matrix12 = MatMul(matrix1,matrix2)
!
!
!      end if
!
      return
      End Subroutine MatrixProduct3x3
