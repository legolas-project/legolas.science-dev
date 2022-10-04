! =============================================================================
!> Submodule containing the implementation of the QR-cholesky algorithm.
!! Using LAPACKS's <tt>zpbtrf</tt> and BLAS's <tt>zgbtrs</tt>, the original
!! problem is written as a standard eigenvalue problem as
!! $$ U^{-H}\mathcal{A}U^{-1}\textbf{Y} = \omega\textbf{Y}\ $$
!! where \(\mathcal{B} = U^HU\) is positive definite and
!! \(\textbf{Y}=U\textbf{X}\).
!! Eventually a call to LAPACK's <tt>zgeev</tt> routine is done to obtain
!! all eigenvalues and eigenvectors.
submodule (mod_solvers) smod_qr_cholesky
  use mod_global_variables, only: dim_subblock
  use mod_banded_matrix_hermitian, only: hermitian_banded_matrix_t
  use mod_transform_matrix, only: matrix_to_hermitian_banded
  implicit none

contains

  !> Solves the eigenvalue problem by rewriting it to a standard form
  !! through splitting of the B-matrix.
  !! @warning Throws an error if <tt>matrix_A</tt> or <tt>matrix_B</tt>
  !!          is not a square matrix. @endwarning
  module procedure qr_cholesky
    !> banded \(U^HU\)
    type(hermitian_banded_matrix_t) :: UU
    !> matrix \(U^{-H}A^U{-1}\)
    complex(dp) :: UAU(matrix_A%matrix_dim, matrix_A%matrix_dim)
    !> order of matrix \(U^{-H}A^U{-1}\)
    integer     :: N
    !> leading dimension of matrix \(U^{-H}A^U{-1}\)
    integer     :: ldUAU
    !> calculate left eigenvectors if "V", omit if "N"
    character   :: jobvl
    !> leading dimension of vl
    integer     :: ldvl
    !> calculate right eigenvectors if "V", omit if "N"
    character   :: jobvr
    !> leading dimension of vr
    integer     :: ldvr
    !> dimension of work array
    integer     :: lwork
    !> work array
    complex(dp), allocatable  :: work(:)
    !> info parameter, 0 on successful exit
    integer     :: info
    !> second work array
    real(dp), allocatable     :: rwork(:)
    !> dummy for left eigenvectors, jobvl = "N" so this is never referenced
    complex(dp) :: vl(2, 2)
    !> number of subdiagonals
    integer :: kl
    !> number of superdiagonals
    integer :: ku

    integer :: i

    ! check input sanity
    if (.not. (matrix_A%matrix_dim == matrix_B%matrix_dim)) then
      call log_message("A or B not square, or not compatible", level="error")
    end if
    ! set array dimensions
    N = size(UAU, dim=1)
    ldUAU = N
    ldvl = N
    ldvr = N
    ! compute B = U^HU
    call log_message("computing B = U^HU", level="debug")
    ! first convert B to banded
    call matrix_B%get_nb_diagonals(ku=ku, kl=kl)
    if (ku /= kl) then
      call log_message( &
        "B has unequal sub/super diagonals: " // str(kl) // "/" // str(ku), &
        level="error" &
      )
      return
    end if
    ! NOTE: We need a banded hermitian matrix with uplo = "U" to get U^HU!
    call matrix_to_hermitian_banded(matrix=matrix_B, diags=ku, uplo="U", banded=UU)
    call zpbtrf("U", UU%n, UU%kd, UU%AB, UU%kd+1, info)
    if (info /= 0) then ! LCOV_EXCL_START
      call log_message("zpbtrf failed: B is not positive definite", level="error")
    end if ! LCOV_EXCL_STOP
    ! compute U^{-H}AU^{-1}
    call log_message("computing B = U^{-H}AU^{-1}", level="debug")
    call matrix_to_array(matrix=matrix_A, array=UAU)
    ! first U^{-H}A by solving U^HX = A
    do i = 1, N
      call ztbsv("U", "C", "N", N, UU%kd, UU%AB, UU%kd+1, UAU(1, i), 1)
    end do
    ! second U^{-H}AU^{-1} by solving XU = U^{-H}A
    do i = 1, N
      call ztbsv("U", "T", "N", N, UU%kd, UU%AB, UU%kd+1, UAU(i, 1), N)
    end do
    ! calculate eigenvectors, we don't use the left ones
    jobvl = "N"
    if (should_compute_eigenvectors()) then
      jobvr = "V"
    else
      jobvr = "N"
    end if
    ! allocate rwork array
    allocate(rwork(2 * N))
    ! get lwork
    allocate(work(1))
    call zgeev( &
      jobvl, jobvr, N, UAU, ldUAU, omega, vl, ldvl, vr, ldvr, work, -1, rwork, info &
    )
    lwork = int(work(1))
    deallocate(work)
    ! allocate work array
    allocate(work(lwork))


    ! solve eigenvalue problem
    call log_message("solving evp using QR algorithm zgeev (LAPACK)", level="debug")
    call zgeev( &
      jobvl, jobvr, N, UAU, ldUAU, omega, vl, ldvl, vr, ldvr, work, lwork, rwork, info &
    )
    if (info /= 0) then ! LCOV_EXCL_START
      call log_message("LAPACK routine zgeev failed!", level="warning")
      call log_message( &
        "value for the info parameter: " // str(info), &
        level="warning", &
        use_prefix=.false. &
      )
    end if ! LCOV_EXCL_STOP

    ! tear down work arrays
    deallocate(work)
    deallocate(rwork)

    ! convert from Y to X by solving UX = Y
    if (should_compute_eigenvectors()) then
      call log_message("computing evs as X = U^{-1}Y", level="debug")
      do i = 1, N
        call ztbsv("U", "N", "N", N, UU%kd, UU%AB, UU%kd+1, vr(1, i), 1)
      end do
    end if

    ! tear down U^HU
    call UU%destroy()
  end procedure qr_cholesky

end submodule smod_qr_cholesky
