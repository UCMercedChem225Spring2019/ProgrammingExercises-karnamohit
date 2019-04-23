      program pgrm_03_01
!
!     This program computes the inverse square-root of a matrix.
!
!     At execution time, the program expects 2 command line arguments: (1) nDim;
!     and (2) the matrix being raised to the (1/2) power.
!
!
      use MatOps
!
      implicit none
      integer,parameter::unitIn=10
      integer::i,iError,nDim,lenSym
      real,dimension(:),allocatable::inputSymMatrix
      real,dimension(:,:),allocatable::inputSqMatrix,invSqrtInputMatrix
      character(len=256)::cmdlineArg
!
!
!     Begin by reading the leading dimension of the matrix and the input file
!     name from the command line. Then, open the file and read the input matrix,
!     inputSymMatrix
!
      call Get_Command_Argument(1,cmdlineArg)
      read(cmdlineArg,'(I)') nDim
      lenSym = (nDim*(nDim+1))/2
      allocate(inputSymMatrix(lenSym),inputSqMatrix(nDim,nDim),  &
        invSqrtInputMatrix(nDim,nDim))
      call Get_Command_Argument(2,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) inputSymMatrix(i)
      endDo
      close(Unit=unitIn)
!
!     Form the square-root of inputSymMatrix. The result is loaded into square
!     (full storage) matrix invSqrtInputMatrix.
!
      write(*,*)' The matrix loaded (column) upper-triangle packed:'
      call SymmetricPacked2Matrix_UpperPac(nDim,inputSymMatrix,  &
        inputSqMatrix)
      write(*,*)' Input Matrix:'
      call Print_Matrix_Full_Real(inputSqMatrix,nDim,nDim)
      call InvSQRT_SymMatrix(nDim,inputSymMatrix,invSqrtInputMatrix)
      write(*,*)' Inverse SQRT Matrix:'
      call Print_Matrix_Full_Real(invSqrtInputMatrix,nDim,nDim)
      write(*,*)' Matrix product that should be the identity.'
      call Print_Matrix_Full_Real(MatMul(MatMul(invSqrtInputMatrix,  &
        invSqrtInputMatrix),inputSqMatrix),nDim,nDim)
!
      end program pgrm_03_01
      Module MatOps
      Implicit None
!
      Contains
!
      Subroutine Print_Matrix_Full_Real(AMat,M,N)
!
!     This subroutine prints a real matrix that is fully dimension -
!     i.e.,
!     not stored in packed form. AMat is the matrix, which is
!     dimensioned
!     (M,N).
!
!     The output of this routine is sent to unit number 6 (set by the
!     local
!     parameter integer IOut).
!
!
!     Variable Declarations
!
      implicit none
      integer,intent(in)::M,N
      real,dimension(M,N),intent(in)::AMat
!
!     Local variables
      integer,parameter::IOut=6,NColumns=5
      integer::i,j,IFirst,ILast
!
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
!
      Do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        Do i = 1,M
          write(IOut,2010) i,(AMat(i,j),j=IFirst,ILast)
        endDo
      endDo
!
      Return
      End Subroutine Print_Matrix_Full_Real


      Subroutine SymmetricPacked2Matrix_LowerPac(N,ArrayIn,AMatOut)
!
!     This subroutine accepts an array, ArrayIn, that is (N*(N+1))/2
!     long.
!     It then converts that form to the N-by-N matrix AMatOut taking
!     ArrayIn to be in lower-packed storage form. Note: The storage mode
!     also assumes the lower-packed storage is packed by columns.
!
      Implicit None
      Integer,Intent(In)::N
      Real,Dimension((N*(N+1))/2),Intent(In)::ArrayIn
      Real,Dimension(N,N),Intent(Out)::AMatOut
!
      Integer::i,j,k
!
!     Loop through the elements of AMatOut and fill them appropriately
!     from
!     Array_Input.
!
      k = 0
      do i=1,N
        do j=i,N
          AMatOut(j,i) = ArrayIn(k + j)
          AMatOut(i,j) = AMatOut(j,i)
        enddo
        k = k+N-i
      enddo
!
! *************************************************************************
! WRITE CODE HERE TO READ THE ARRAY ELEMENTS FROM THE INPUT FILE.
! *************************************************************************
!
!
      Return
      End Subroutine SymmetricPacked2Matrix_LowerPac


      Subroutine SymmetricPacked2Matrix_UpperPac(N,ArrayIn,AMatOut)
!
!     This subroutine accepts an array, ArrayIn, that is (N*(N+1))/2
!     long.
!     It then converts that form to the N-by-N matrix AMatOut taking
!     ArrayIn to be in upper-packed storage form. Note: The storage mode
!     also assumes the upper-packed storage is packed by columns.
!
      Implicit None
      Integer,Intent(In)::N
      Real,Dimension((N*(N+1))/2),Intent(In)::ArrayIn
      Real,Dimension(N,N),Intent(Out)::AMatOut
!
      Integer::i,j,k
!
!     Loop through the elements of AMatOut and fill them appropriately
!     from
!     Array_Input.
!
      k = 0
      do i=1,N
        do j=1,i
          AMatOut(i,j) = ArrayIn(k + j)
          AMatOut(j,i) = AMatOut(i,j)
        enddo
        k = k+i
      enddo
!
! *************************************************************************
! WRITE CODE HERE TO READ THE ARRAY ELEMENTS FROM THE INPUT FILE.
! *************************************************************************
!
!
      Return
      End Subroutine SymmetricPacked2Matrix_UpperPac


      Subroutine InvSQRT_SymMatrix(nDim,ArrayIn,invSqrtofArrayIn)
!
!     This subroutine accepts an array, ArrayIn, which is an
!     upper-packed symmetric matrix, (nDim x nDim) in size. The matrix
!     is in packed form, so is actually a 1-D array, nDim*(nDim+1)/2 in length.
!     It then diagonalizes this matrix to give back a matrix, raised to
!     (-1/2), (nDim x nDim) in length.
!
      implicit none
      Integer,Intent(In)::nDim
      Real,Dimension(nDim*(nDim+1)/2),Intent(In)::ArrayIn
      Real,Dimension(nDim,nDim),Intent(Out)::invSqrtofArrayIn
!
      Integer::i,IError
      Real,Dimension(nDim,nDim)::EVecs,TempMat
      Real,Dimension(nDim)::EVals
      Real,Dimension(3*nDim)::Temp_Vector
!
      Call SSPEV('V','U',nDim,ArrayIn,EVals,EVecs,nDim,  &
        Temp_Vector,IError)
      do i=1,nDim
        EVals(i) = 1.0/sqrt(EVals(i))
        invSqrtofArrayIn(i,i) = EVals(i)
      enddo
      invSqrtofArrayIn = matmul(EVecs,matmul(invSqrtofArrayIn,  &
                                             transpose(EVecs)))
!
      Return
      End Subroutine InvSQRT_SymMatrix


      Subroutine Sq2SymMatrix_UpCol(nDim,ArrayIn,symmtrzdArrayIn)
!
!     This subroutine accepts a symmetric square matrix ArrayIn of dimensions
!     (nDim x nDim) and symmetrizes it in upper-triangular column-packed
!     form.
!
      Integer,Intent(In)::nDim
      Real,Dimension(nDim,nDim),Intent(In)::ArrayIn
      Real,Dimension(nDim*(nDim+1)/2),Intent(Out)::symmtrzdArrayIn
!
      Integer::i,j,k
!
      k = 0
      do i=1,nDim
        do j=1,i
          symmtrzdArrayIn(k + j) = ArrayIn(j,i)
        enddo
        k = k+i
      enddo
!
!      k = 0
!      do i=1,N
!        do j=i,N
!          AMatOut(j,i) = ArrayIn(k + j)
!          AMatOut(i,j) = AMatOut(j,i)
!        enddo
!        k = k+N-i
!      enddo
!
      Return
      End Subroutine Sq2SymMatrix_UpCol
!
      End Module MatOps
