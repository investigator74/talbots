module ic
    use, intrinsic :: iso_c_binding
    use fourier

    implicit none

    namelist /param/ nz, period, lz, lx, nx, nk, nth, delta, a0_peak, xcp, alfa, g_amp, g_x0, g_x1, r0, r1, sigma, gamma, xe, c, &
        x_out, intrvl, in_type, recount, thout, central_mirror, cont, amp_only, it_todo, norm, ops, lambda, s2k
        real(c_double) h, lz, lx, hx, hth, delta, a0_peak, imp_x0, imp_xsp, g_amp, g_x0, g_x1, r0, r1, x_out, sigma, xe, gamma, c, c3, xcp, alfa, norm, &
        kappa, lambda, kk
        integer(c_int) :: nx, nk, nth, nz, iimp_x0, iimp_xend, ix_out, intrvl, it_todo, it_doiter, in_type, it_made = 0, it_flag =0, ie, ig0, ig1, ixe1, ixe2
    real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)
    logical(c_bool) recount, thout, central_mirror, period, amp_only, cont, ops, s2k

    !namelist /perenorm_param/ nz, period, lz, lx, nx, nk, nth, delta, a0_peak, xcp, alfa, g_amp, g_x0, g_x1, r0, r1, sigma, gamma, xe, c, &
    ! x_out, intrvl, in_type, recount, thout, central_mirror, cont, amp_only, it_todo, norm

        complex(c_double_complex), allocatable :: a1(:), a0(:), ak1(:), ak0(:), atmp(:), jk1(:), jk0(:), k(:), ex( :), dlt(:), tmp(:), aktmp(:), akzl(:), &
                                              k2(:), akz0(:), a0z0(:), a0z0cut(:)
    real(c_double), allocatable :: th0(:, :), th1(:, :), dthdz(:, :), fk1(:), fk2(:), rhs0(:, :), z(:), x(:), &
                                   a_amp_z0(:), a_amp_zl(:), a_spec_amp_z0(:), a_spec_amp_zl(:), g(:), &
                                   sum_abs2_a_plus_by_z(:), sum_abs2_a_plus_by_z_k(:)!, theta(:,:,:)
    integer(c_int), allocatable :: it(:)

    complex(c_double_complex), parameter :: im1 = (0.0d0, 1.0d0)
    real(c_double) :: start_time, finish_time, smirr, soutm

    !interface
    ! subroutine fn_for_fortran_to_call(ptr) &
    ! bind(c, name='fn_for_fortran_to_call')
    ! use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    ! implicit none
    ! type(c_ptr), intent(in), value :: ptr
    ! end subroutine fn_for_fortran_to_call
    !end interface

contains
    subroutine calc_idx()

        implicit none

        kk = 4.0d0*pi/lambda
        if (s2k .eq. .true.) then
            !lx = lx / sqrt(kk)

            !what depends on the new lx
            xcp = 0.5d0*lx! *norm** (1.0/6.0)
            alfa = 0.5d0*lx/100.0d0! *norm** (1.0/6.0)
            g_x0 = 0.35*lx! *norm** (1.0/6.0)
            g_x1 = 0.65*lx! *norm** (1.0/6.0)
            !xe = 0.166666666666667 * lx! *norm** (1.0/6.0)
            xe = 0.22*lx! *norm** (1.0/6.0)
            x_out = 0.5d0*lx! *norm** (1.0/6.0)

            s2k = .false.
        end if

        if (norm .ne. 0.0d0) then
            lx = lx*norm**(1.0/6.0)

            !what depends on the new lx
            xcp = 0.5d0*lx! *norm** (1.0/6.0)
            alfa = 0.5d0*lx/100.0d0! *norm** (1.0/6.0)
            g_x0 = 0.35*lx! *norm** (1.0/6.0)
            g_x1 = 0.65*lx! *norm** (1.0/6.0)
            !xe = 0.166666666666667 * lx! *norm** (1.0/6.0)
            xe = 0.22*lx! *norm** (1.0/6.0)
            x_out = 0.5d0*lx! *norm** (1.0/6.0)

            delta = delta/norm**(1.0/3.0)
            a0_peak = a0_peak/norm**(2.0/3.0)
            sigma = sigma/norm**(1.0/3.0)
            c = c/norm
            if (period == .true.) then
                lz = lz*(lx*lx)/lambda
            end if
            period = .false.
            norm = 0
        else
            if (period == .true.) then
                lz = lz*(lx*lx)/lambda
            end if
            period = .false.
        end if

        open (unit=1, file='input_fortran_real.dat')
        write (unit=1, nml=param)
        close (unit=1)

        !c3 = c*c*c
        c3 = c
        h = lz/nz
        nz = nz + 1

        hth = 2.0d0*pi/nth
        hx = lx/nx

        ixe1 = xe/hx + 1
        ixe2 = nx - ixe1 + 2
        xe = (ixe1 - 1)*hx !xe clarification

        !print *, 'ixe1 = ', ixe1
        !print *, 'ixe2 = ', ixe2
        !print *, 'xe = ', xe
        !pause

        ix_out = int(x_out/hx)
        if (ix_out <= 0) ix_out = 1
        if (ix_out > nx) ix_out = nx

        iimp_x0 = max(1, int(imp_x0/hx) + 1)
        imp_x0 = (iimp_x0 - 1)*hx !to (iimp_x0 - 1) * hx exactly equal to imp_x0
        iimp_xend = min(nx + 1, int((imp_x0 + imp_xsp)/hx) + 1) ! calculate for nx'= nx + 1
        imp_xsp = (iimp_xend - iimp_x0)*hx !for point precision
        if (iimp_xend == nx + 1) iimp_xend = iimp_xend - 1 ! last interval point not involved

        !what is the iteration
        if (it_flag == 0) it_made = 0
        it_doiter = it_made + it_todo
        !it_todo = 0
    end subroutine

    subroutine calc_theta(th, dthdz)

        implicit none

        real(c_double), intent(inout) :: th(:, :), dthdz(:, :)
        integer, dimension(size(th)) :: i
        integer ix

        i = (/1:size(th, 1)/)

        do ix = 1, 2
            th(:, ix) = hth*(i - 1)
            dthdz(:, ix) = delta
        end do

        !open(777, file='test.dat')
        !doix=1,nth
        ! write(777,*) ix-1, th(ix,1), th(ix,2)
        !enddo
        !close(777)
        !stop
    end subroutine
end module ic

program elektron2dsi
    use ic; use fourier

    implicit none

    interface
        subroutine init() bind(c, name='init')
        end subroutine init
        subroutine finish() bind(c, name='finish')
        end subroutine finish
        subroutine makea(atmp, ak)
            use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
            complex(c_double_complex), dimension(:), intent(inout) :: atmp
            complex(c_double_complex), dimension(:), intent(in) :: ak
        end subroutine makea
        subroutine makeak(ak, atmp)
            use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
            complex(c_double_complex), dimension(:), intent(inout) :: ak
            complex(c_double_complex), dimension(:), intent(in) :: atmp
        end subroutine makeak
        subroutine make_a0z0_ak1_atmp(a0z0, ak1, atmp, ak)
            use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
            complex(c_double_complex), dimension(:), intent(inout) :: a0z0, ak1, atmp
            complex(c_double_complex), dimension(:), intent(in) :: ak
        end subroutine make_a0z0_ak1_atmp
    end interface

    !type(c_ptr), intent(in), value :: ptr
    real(c_double) eff_tmp, eff_tmp_k, eff_tmp_b, eff_tmp_k_b, eff(2), &
        sum_eff, loss_on_the_way_plus, loss_on_the_way_plus_k, loss_on_the_way_minus, &
     int_abs2_a_plus_at_z0, int_abs2_a_plus_at_z0_k, int_abs2_a_plus_at_zl, int_abs2_a_plus_at_zl_k, int_abs2_a_plus_at_zl_on_mir, &
        int_abs2_a_plus_at_zl_out_mir, int_abs2_a_plus_at_zl_out_mir_k, int_abs2_a_minus_at_z0_out_mir_k, &
            int_abs2_a_minus_at_z0_k, int_abs2_a_minus_at_z0_on_mir, int_abs2_a_minus_at_z0_out_mir, int_abs2_a_minus_at_zl_on_mir, loss_on_the_way_minus_k, &
        int_abs2_a_minus_at_zl_on_mir_k, int_obrezki_z0, int_obrezki_zl, int_abs2_a_minus_at_z0, int_abs2_a_minus_at_zl, obrezki
    integer(c_int) iz, percent, rec_length, nr, first_it, ith
    character*6 str
    integer i

    !chk_arr = check_alloc()
    !if (chk_arr /= .true.) stop 'allocation error!'

    !read out the parameters, calculate a0(x), f(x), k2, z, x, etc.
    call init()

    !what does not depend on j
    dlt = im1*gamma*k2 + sigma
    ex = cdexp(-dlt*h)

    !open(17, file='test.dat')
    !do i=1,nk
    ! write(17, '(i, 2f17.8)') i, dlt(i)
    !enddo
    !close(17)
    !print *, size(dlt), size(ex), nk, size(k2)
    !pause
    !stop

    rec_length = c_double_complex*2*nx/4
    !open files
    open (221, file='$$$z0.bin', access='direct', recl=rec_length, err=101)
    open (222, file='$$$zl.bin', access='direct', recl=rec_length, err=101)
    open (3, file='eff.dat', err=101)
    open (35, file='eff_b.dat', err=101)
    open (53, file='eff_new.dat', err=101)
    open (33, file='eff_k.dat', err=101)
    open (335, file='eff_k_b.dat', err=101)
    open (533, file='eff_new_k.dat', err=101)
    open (788, file='it.dat', err=101)
    open (777, file='for_graphics.dat', err=101)

    call cpu_time(start_time)

    first_it = it_made + 1
    do_iter: do it_made = first_it, it_doiter
        percent = 0
        write (*, '(a11,\,i5,\,a2,i5,a,\)') 'iteration #', it_made, ': ', percent, '%'

        !initial conditions for the phase and its distunction at z=0
        call calc_theta(th0, dthdz)

        !theta(:, :, 1) = th0 !initial position of electrons

        !calculation of efficiency at point z=0
        !doix=1,nx
        ! eff(ix) = 1.0d0 / nth * sum(dthdz(:,ix) - delta) !counting efficiency
        !end do
        eff(1) = 1.0d0/nth*sum(dthdz(:, 1) - delta) !(xcp1)
        eff(2) = 1.0d0/nth*sum(dthdz(:, 2) - delta) !(xcp2)

        !initial current (z=0)
        jk0 = fk1*jf(th0(:, 1)) + fk2*jf(th0(:, 2))

        !ak0
        !ak0 = fs(a0)
        call makeak(ak0, a0)

        !setting the file index
        inquire (222, nextrec=nr)
        !write a(z=0,x)
        write (221, rec=nr) a0
        !print *, 'nr = ', nr

        if (((intrvl > 0) .and. (mod(it_made, intrvl) == 0)) &
            .or. (it_made == first_it) &
            .or. (it_made == it_doiter)) then
            a_amp_z0 = cdabs(a0)
            a_spec_amp_z0 = cdabs(fs(a0))
            write (str, 105) it_made
105         format(i0.6)
            str = trim(str)
            open (1, file='a_'//str//'.bin', form='binary', err=101)
            open (2, file='ak_'//str//'.bin', form='binary', err=101)
            open (88, file='jk_'//str//'.bin', form='binary', err=101)
            write (1, err=103) a0z0
            write (2, err=103) ak0
            write (88, err=103) jk0
            if (amp_only == .false.) then
                open (8, file='eff_'//str//'.bin', form='binary', err=101)
                write (8, err=103) eff
            end if
            if (thout == .true.) then
                open (10, file='th_'//str//'.dat', err=101)
                write (10, '(e17.8,\)', err=103) 0.0
                do ith = 1, nth
                    write (10, '(e17.8,\)', err=103) th0(ith, 1)
                end do
                do ith = 1, nth
                    write (10, '(e17.8,\)', err=103) th0(ith, 2)
                end do
                write (10, '(/,\)')
            end if
        end if

        !if (((intrvl > 0) .and.(mod(it_made, intrvl) == 0)) &
        ! .or. (it_made == first_it) &
        ! .or. (it_made == it_doiter)) then
        ! atmp=ifs(atmp)
        ! open(777, file='test_' // str // '.dat')
        ! do ix=1,nx
        ! write(777,'(9e27.18)') (ix-1)*hx, cdabs(a0(ix)), g(ix), atmp(ix), atmp(ix) * g(ix), a0( ix)
        ! enddo
        ! close(777)
        !endif

        !to calculate the efficiency we calculate the sum of abs(a(z=0))**2 by x
        int_abs2_a_plus_at_z0 = sum(cdabs(a0)*cdabs(a0))*hx
        int_abs2_a_plus_at_z0_k = 0.5d0*sum(cdabs(ak0)*cdabs(ak0))

        !sum of field amplitudes at z=0 (also for effectiveness)
        sum_abs2_a_plus_by_z = cdabs(a0)*cdabs(a0)
        sum_abs2_a_plus_by_z_k = 0.5d0*cdabs(ak0)*cdabs(ak0)

        do_z: do iz = 1, nz - 1
            rhs0 = rhs(ak0, th0)
            th1 = th0 + dthdz*h + h/2.0d0*rhs0*h !theta predictor
            jk1 = fk1*jf(th1(:, 1)) + fk2*jf(th1(:, 2)) !current predictor

            !predictor a (interpolation)
            ak1(1) = ak0(1) + h/2.0d0*(jk0(1) + jk1(1))
            ak1(2:nk) = ak0(2:nk)*ex(2:nk) + &
                        c3*(jk0(2:nk) + jk1(2:nk)*(-1.0d0 + dlt(2:nk)*h) + &
                            ex(2:nk)*(jk1(2:nk) - jk0(2:nk)*(1.0d0 + dlt(2:nk)*h)))/dlt(2:nk)/dlt(2:nk)/h

            !predictor a (keystone)
            !atmp = (ak0 + c3 * h / 2.0d0 * jk0) * cdexp(-dlt * h) !part a
            !ak1 = atmp + c3 * h / 2.0d0 * jk1 !predictor a

            !a1 = ifs(ak1) !back to reality
            call makea(a1, ak1)

            !corrector theta - 1
            th1 = th0 + dthdz*h + h/6.0d0*rhs0*h &
                  + h/3.0d0*rhs(ak1, th1)*h

            !theta(:, :, iz + 1) = th1 !write down the electron trajectory

            jk1 = fk1*jf(th1(:, 1)) + fk2*jf(th1(:, 2)) !current corrector

            !what depends on the current
            !jkd = jk1 - jk0

            !offset a (interpolation)
            ak1(1) = ak0(1) + h/2.0d0*(jk0(1) + jk1(1))
            ak1(2:nk) = ak0(2:nk)*ex(2:nk) + &
                        c3*(jk0(2:nk) + jk1(2:nk)*(-1.0d0 + dlt(2:nk)*h) + &
                            ex(2:nk)*(jk1(2:nk) - jk0(2:nk)*(1.0d0 + dlt(2:nk)*h)))/dlt(2:nk)/dlt(2:nk)/h

            !corrector a (keystone)
            !atmp = (ak0 + c3 * h / 2.0d0 * jk0) * cdexp(-dlt * h) !part a
            !ak1 = atmp + c3 * h / 2.0d0 * jk1 !corrector a

            dthdz = dthdz + h/2.0d0*(rhs0 + rhs(ak1, th1))
            !corrector theta - 1 end

            !calculation of efficiency at point z = iz*h
            !doix=1,nx
            ! eff(ix) = 1.0d0 / nth * sum(dthdz(:,ix) - delta) !counting efficiency
            !end do
            eff(1) = 1.0d0/nth*sum(dthdz(:, 1) - delta) !counting efficiency (xcp1)
            eff(2) = 1.0d0/nth*sum(dthdz(:, 2) - delta) !counting efficiency (xcp2)

            !back to reality
            !a1 = ifs(ak1)
            call makea(a1, ak1)

            !sum of field amplitudes at z=iz*h
            sum_abs2_a_plus_by_z = sum_abs2_a_plus_by_z + cdabs(a1)*cdabs(a1)
            sum_abs2_a_plus_by_z_k = sum_abs2_a_plus_by_z_k + 0.5d0*cdabs(ak1)*cdabs(ak1)

            if (((intrvl > 0) .and. (mod(it_made, intrvl) == 0)) &
                .or. (it_made == first_it) &
                .or. (it_made == it_doiter)) then
                write (1, err=103) a1
                write (2, err=103) ak1
                write (88, err=103) jk1
                if (amp_only == .false.) then
                    write (8, err=103) eff
                end if
                if (thout == .true.) then
                    write (10, '(e17.8,\)', err=103) iz*h
                    do ith = 1, nth
                        write (10, '(e17.8,\)', err=103) th1(ith, 1)
                    end do
                    do ith = 1, nth
                        write (10, '(e17.8,\)', err=103) th1(ith, 2)
                    end do
                    write (10, '(/,\)')
                end if
            end if

            th0 = th1 !for the next z step
            jk0 = jk1 !for the next step in z
            ak0 = ak1 !for the next step in z
            !dthdz0 = dthdz1 !for the next z step

            percent = int(real(iz)/real(nz - 2)*100 + 0.5)
            write (*, '(\,a6,i5,a)') '\b\b\b\b\b\b'c, percent, '%'
        end do do_z

        !calculation of the sum of efficiency for x at the point z=lz
        sum_eff = (eff(1) + eff(2))/2.0d0

        !to calculate the efficiency we calculate the sum of abs(a(z=lz))**2 by x
        int_abs2_a_plus_at_zl = sum(cdabs(a1)*cdabs(a1))*hx
        int_abs2_a_plus_at_zl_k = 0.5d0*sum(cdabs(ak1)*cdabs(ak1))

        int_abs2_a_plus_at_zl_on_mir = sum(cdabs(a1(ig0:ig1 - 1))*cdabs(a1(ig0:ig1 - 1)))*hx
   int_abs2_a_plus_at_zl_out_mir = (sum(cdabs(a1(1:ig0-1))*cdabs(a1(1:ig0-1))) + sum(cdabs(a1(ig1+1:nx))*cdabs(a1(ig1+ 1:nx)))) * hx

        aktmp = a1*g
        call sint(aktmp)
!subtract the uncut modes of the cut reflection
        int_abs2_a_plus_at_zl_out_mir_k = int_abs2_a_plus_at_zl_k - 0.5d0*sum(cdabs(aktmp)*cdabs(aktmp))

        aktmp(1:nk) = dcmplx(0)
        int_obrezki_zl = 0.5d0*dmysum(cdabs(aktmp)*cdabs(aktmp)) ! energy in cut to zl

        !double sum of field amplitudes in z and in x
        loss_on_the_way_plus = 2.0d0*sigma*sum(sum_abs2_a_plus_by_z)*hx*h
        loss_on_the_way_plus_k = 2.0d0*sigma*sum(sum_abs2_a_plus_by_z_k)*h

        !write a(z=l,x)
        write (222, rec=nr) a1
        write (788, *) it_made

        if (it_made == it_doiter) then
            a_amp_zl = cdabs(a1)
            a_spec_amp_zl = cdabs(fs(a1))
        end if

        !initial conditions for the next iteration
        !if (it_made < it_doiter) then
        if (recount == .false.) then
            a0 = r1*a1*g*r0
        else
            atmp = r1*a1*g!! a minus on aperture z \u003d lz (real)

            !int_abs2_a_minus_at_zl_on_mir = sum(cdabs(atmp(ig0:ig1-1))*cdabs(atmp(ig0:ig1-1))) * hx
            !int_abs2_a_minus_at_zl_on_mir = sum(cdabs(a1(ig0:ig1-1))*cdabs(a1(ig0:ig1-1))) * hx

            !fourier
            call makeak(akzl, atmp)

            !aktmp = dcmplx(0)
            !aktmp(1:nk) = akzl
            !call isint(aktmp)
            call makea(atmp, akzl)

            int_abs2_a_minus_at_zl_on_mir_k = 0.5d0*dmysum(cdabs(akzl)*cdabs(akzl)) !for efficiency (what is already reflected from the mirror)
            int_abs2_a_minus_at_zl_on_mir = sum(cdabs(atmp)*cdabs(atmp))*hx

            !open(777, file='test.dat')
            !do i=1,nx
            ! write(777, '(3e17.8)') (i-1)*hx, atmp(i)
            !enddo
            !close(777)
            !
            !print *, int_abs2_a_minus_at_zl_on_mir_k
            !print *, int_abs2_a_minus_at_zl
            !pause
            !stop

            !------------------------------------------------- -------------------------------------------------- -------------------------------------------------- --

            akz0 = akzl*cdexp(-dlt*lz) !what is returned to z0

            int_abs2_a_minus_at_z0_k = 0.5d0*dmysum(cdabs(akz0)*cdabs(akz0)) !a minus on aperture z=0 (real)

            call make_a0z0_ak1_atmp(a0z0, a0z0cut, ak0, akz0)

            int_abs2_a_minus_at_z0 = sum(cdabs(a0z0)*cdabs(a0z0))*hx

            !akz0 fashions of what is going to the beginning
            !ak0 cut mode of reflection cut

            loss_on_the_way_minus_k = int_abs2_a_minus_at_zl_on_mir_k - int_abs2_a_minus_at_z0_k

            aktmp = a0z0*g
            call sint(aktmp) ! uncut reflection cut modes
            int_abs2_a_minus_at_z0_out_mir_k = int_abs2_a_minus_at_z0_k - 0.5d0*sum(cdabs(aktmp)*cdabs(aktmp))

            aktmp(1:nk) = dcmplx(0)
            int_obrezki_z0 = 0.5d0*dmysum(cdabs(aktmp)*cdabs(aktmp)) ! energy in cut to z0

            int_abs2_a_minus_at_z0_on_mir = sum(cdabs(a0z0(ig0:ig1 - 1))*cdabs(a0z0(ig0:ig1 - 1)))*hx
                    int_abs2_a_minus_at_z0_out_mir = (sum(cdabs(a0z0(1:ig0-1))*cdabs(a0z0(1:ig0-1))) + sum(cdabs(a0z0(ig1+1:nx))*cdabs(a0z0(ig1+ 1:nx)))) * hx
            loss_on_the_way_minus = int_abs2_a_minus_at_zl_on_mir - (int_abs2_a_minus_at_z0_on_mir + int_abs2_a_minus_at_z0_out_mir)

            !reflection and cut at the first mirror
            a0 = a0z0cut
            open (388, file='a0.bin', form='binary', err=101)
            write (388) a0
            close (388)
        end if
        !endif

        !recording efficiency at this iteration
        eff_tmp = (loss_on_the_way_plus + int_abs2_a_plus_at_zl - int_abs2_a_plus_at_z0)/lx - 4.0d0*c3*sum_eff
        eff_tmp_k = loss_on_the_way_plus_k + int_abs2_a_plus_at_zl_k - int_abs2_a_plus_at_z0_k - 4.0d0*c3*sum_eff

        obrezki = int_obrezki_z0 + int_obrezki_zl

        write (3, 104, err=103) & ! eff.dat
            it_made, &
            int_abs2_a_plus_at_zl/lx, &
            loss_on_the_way_plus/lx, &
            int_abs2_a_plus_at_z0/lx, &
            sum_eff, &
            int_obrezki_z0, &
            int_obrezki_zl, &
            eff_tmp, &
            eff_tmp/(4.0d0*c3*sum_eff), &
            cdabs(a1(ix_out))

        eff_tmp_b = (loss_on_the_way_minus + int_abs2_a_minus_at_z0 - int_abs2_a_minus_at_zl_on_mir)/lx

        write (35, 108, err=103) & ! eff_b.dat
            it_made, &
            loss_on_the_way_minus/lx, &
            int_abs2_a_minus_at_zl_on_mir/lx, &
            int_abs2_a_minus_at_z0/lx, &
            eff_tmp_b

108     format(i, 4e17.8)

        write (33, 104, err=103) & ! eff_k.dat
            it_made, &
            int_abs2_a_plus_at_zl_k, &
            loss_on_the_way_plus_k, &
            int_abs2_a_plus_at_z0_k, &
            sum_eff, &
            int_obrezki_z0, &
            int_obrezki_zl, &
            eff_tmp_k, &
            eff_tmp_k/(4.0d0*c3*sum_eff), &
            cdabs(a1(ix_out))

104     format(i, 9e17.8)

        eff_tmp_k_b = loss_on_the_way_minus_k + int_abs2_a_minus_at_z0_k - int_abs2_a_minus_at_zl_on_mir_k

        write (335, 108, err=103) & ! eff_k_b.dat
            it_made, &
            loss_on_the_way_minus_k, &
            int_abs2_a_minus_at_zl_on_mir_k, &
            int_abs2_a_minus_at_z0_k, &
            eff_tmp_k_b

        write (53, 107, err=103) it_made, & ! eff_new.dat
            int_abs2_a_minus_at_z0_out_mir/lx, &
            int_abs2_a_plus_at_zl_out_mir/lx, &
            (loss_on_the_way_minus_k + loss_on_the_way_plus/lx), &
            obrezki, &
            4.0d0*c3*sum_eff, &
            loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir + &
                                       loss_on_the_way_plus + &
                                       int_abs2_a_minus_at_z0_out_mir)/lx + obrezki - &
            4.0d0*c3*sum_eff, &
            (loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir + &
                                        loss_on_the_way_plus + &
                                        int_abs2_a_minus_at_z0_out_mir)/lx + obrezki - &
             4.0d0*c3*sum_eff)/4.0d0/c3/sum_eff

        write (533, 107, err=103) it_made, & ! eff_new_k.dat
            int_abs2_a_minus_at_z0_out_mir_k, &
            int_abs2_a_plus_at_zl_out_mir_k, &
            (loss_on_the_way_minus_k + loss_on_the_way_plus_k), &
            obrezki, &
            4.0d0*c3*sum_eff, &
            loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir_k + &
                                       loss_on_the_way_plus_k + &
                                       int_abs2_a_minus_at_z0_out_mir_k + obrezki) - &
            4.0d0*c3*sum_eff, &
            (loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir_k + &
                                        loss_on_the_way_plus_k + &
                                        int_abs2_a_minus_at_z0_out_mir_k + obrezki) - &
             4.0d0*c3*sum_eff)/4.0d0/c3/sum_eff

107     format(i, 7e17.8)

        write (777, '(8e17.8)', err=103) & ! for_graphics.dat
            lz, &
            int_abs2_a_minus_at_z0_out_mir_k, &
            int_abs2_a_plus_at_zl_out_mir_k, &
            (loss_on_the_way_minus_k + loss_on_the_way_plus_k), &
            obrezki, &
            4.0d0*c3*sum_eff, &
            (int_abs2_a_minus_at_z0_out_mir_k + int_abs2_a_plus_at_zl_out_mir_k)/(4.0d0*c3*sum_eff), &
            delta

        !ix_out = 65
        if (recount == .false.) then
            write (*, '(a,\)') char(13)
                write(*,'(a,i6,a,f7.3,a,e17.8,a,e17.8,a)') 'iteration #', it_made, ': a(', x_out, ') = ', cdabs(a1(ix_out)), ' eff = ', sum_eff, ' carryover'
        else
            write (*, '(a,\)') char(13)
                write(*,'(a,i6,a,f7.3,a,e17.8,a,e17.8,a)') 'iteration #', it_made, ': a(', x_out, ') = ', cdabs(a1(ix_out)), ' eff = ', sum_eff, ' recount'
        end if

        !call fn_for_fortran_to_call(ptr)
    end do do_iter

    write (*, '(/,/)')

    !closing write files a(z=0,x) and a(z=l,x)
    close (1)
    close (3)
    close (53)
    close (221)
    close (222)
    close (788)
    if (amp_only == .false.) then
        close (2)
        close (8)
    end if
    if (thout == .true.) close (10)

    call cpu_time(finish_time)
    print *, 'computation time = ', finish_time - start_time, 'seconds'

    call write_result()
    call finish()

    print *, 'calculation finished.'
    pause
    stop
101 stop 'error of file open.'
102 stop 'error of file reading.'
103 stop 'error of file writing.'

contains
    function jf(th)
        implicit none
        real(c_double), intent(in) :: th(:)
        complex(c_double_complex) jf

        jf = 2.0d0/dble(nth)*sum(cdexp(-im1*th))
    end function jf

    function rhs(ak, th)
        implicit none
        complex(c_double_complex), intent(in) :: ak(:)
        real(c_double), intent(in) :: th(:, :)
        real(c_double), dimension(size(th, 1), size(th, 2)) :: rhs

        !rhs(:,1) = dreal(sum(fk1 * ak) * cdexp(im1 * th(:,1)))
        !rhs(:,2) = dreal(sum(fk2 * ak) * cdexp(im1 * th(:,2)))

        rhs(:, 1) = dreal(mysum(fk1*ak)*cdexp(im1*th(:, 1)))
        rhs(:, 2) = dreal(mysum(fk2*ak)*cdexp(im1*th(:, 2)))

        !atmp=ifs(ak)
        !rhs(:,1) = dreal(atmp(ixe1) * cdexp(im1 * th(:,1)))
        !rhs(:,2) = dreal(atmp(ixe2) * cdexp(im1 * th(:,2)))
    end function rhs

    function mysum(a)
        implicit none
        complex(c_double_complex), intent(in) :: a(:)
        complex(c_double_complex) :: mysum
        integer(c_int) i, n

        n = size(a)
        mysum = dcmplx(0)

        do i = n, 1, -1
            mysum = mysum + a(i)
        end do
    end function mysum

    function dmysum(a)
        implicit none
        real(c_double), intent(in) :: a(:)
        real(c_double) :: dmysum
        integer(c_int) i, n

        n = size(a)
        !dmysum = dcmplx(0)
        dmysum = 0.0d0

        do i = n, 1, -1
            dmysum = dmysum + a(i)
        end do
    end function dmysum
    !end subroutine calculate_fortran
end program elektron2dsi

subroutine init() bind(c, name='init')
    use, intrinsic :: iso_c_binding
    use ic

    implicit none

    integer i

    interface
        subroutine read_param() bind(c, name='read_param')
        end subroutine read_param
        function a0_fn_stat() result(a0_res)
            use ic
            complex(c_double_complex), dimension(nx) :: a0_res
        end function a0_fn_stat
        function fk_fn(xe) result(fk_res)
            use, intrinsic :: iso_c_binding, only: c_double
            use ic, only: nk
            real(c_double), dimension(nk) :: fk_res
            real(c_double) xe
        end function fk_fn
        function g_fn() result(g_res)
            use ic
            real(c_double), dimension(nx) :: g_res
        end function g_fn
        function k_fn() result(k)
            use ic, only: nx
            use, intrinsic :: iso_c_binding
            complex(c_double_complex), dimension(2*nx) :: k
        end function k_fn
        function k2_fn() result(k2_res)
            use ic
            complex(c_double_complex), dimension(nk) :: k2_res
        end function k2_fn
        function dn_fn() result(dn_res)
            use ic
            complex(c_double_complex), dimension(nk) :: dn_res
        end function dn_fn
    end interface

    call read_param()
    call calc_idx()
    call allocate_arrays()
    call calc_zxit()
    call sincost_init(nx)
    call fft_init(2*nx)
    call dst_init(nx, lx)

    ! initial conditions for a (z = 0)
    if (cont .eq. .true.) then
        open (1, file='a0.bin', form='binary', err=101)
        read (1) a0
        close (1)
    else
        a0 = a0_fn_stat()
        a0z0 = a0
    end if

    ! open (17, file = 'test.dat')
    ! do i = 1, nx
    ! write (17, '(3f17.8)') (i-1)*hx, a0 (i)
    ! enddo
    ! close (17)
    ! stop

    ! smooth f
    fk1(:) = fk_fn(xe)
    fk2(:) = fk_fn(lx - xe)

    !
    g = g_fn()

    ! k ** 2
    if (ops .eq. .false.) then
        k2 = k2_fn()
    else
        k2 = dn_fn()
    end if

    open (1, file='initag.dat')
    do i = 1, nx
        write (1, '(3e17.8,i10)') (i - 1)*hx, dreal(a0(i)), g(i)
    end do
    close (1)

    open (1, file='initfk.dat')
    do i = 1, nk
        write (1, '(2e17.8,i10)') fk1(i), fk2(i), int(k2(i))
    end do
    close (1)

    write (*, *) 'nz = ', nz
    write (*, *) 'h = ', h

    print *, 'lz = ', lz
    print *, 'lx = ', lx
    print *, 'c3 = ', c3

    !print *, size(k2), size(fk1), size(fk2), size(dlt), size(ex)
    !pause
    !stop

    return
101 stop 'error of file open.'
end subroutine init

subroutine finish() bind(c, name='finish')
    use fourier, only: sincost_destroy, fft_destroy

    call sincost_destroy()
    call fft_destroy()
    call deallocate_arrays()
end subroutine finish

subroutine write_result()
    use ic
    use, intrinsic :: iso_c_binding

    implicit none

    integer i, j

    call cpu_time(start_time)

    it_made = it_made - 1

    open (1, file='aend.dat', err=101)
    do i = 1, nx
        write (1, '(1p3e17.8)', err=103) (i - 1)*hx, a_amp_z0(i), a_amp_zl(i)
    end do
    close (1)

    !open(877, file = 'theta1.dat', err = 101)
    ! do i=1,nz
    ! write(877, '(e17.8,\)') (i - 1) * h
    ! do j=1,nth
    ! write(877, '(e17.8,\)') theta(j, 1, i)
    ! enddo
    ! write(877, '(/,\)')
    ! enddo
    !close(877)
    !
    !open(877, file = 'theta2.dat', err = 101)
    ! do i=1,nz
    ! write(877, '(e17.8,\)') (i - 1) * h
    ! do j=1,nth
    ! write(877, '(e17.8,\)') theta(j, 2, i)
    ! enddo
    ! write(877, '(/,\)')
    ! enddo
    !close(877)

    call cpu_time(finish_time)
    print *, 'writing time = ', finish_time - start_time, ' seconds'

    return
101 stop 'error of file open.'
102 stop 'error of file reading.'
103 stop 'error of file writing.'
end subroutine write_result

subroutine calc_zxit()
    use ic

    implicit none

    integer i

    do i = 1, nz
        z(i) = (i - 1)*h
    end do
    do i = 1, nx
        x(i) = (i - 1)*hx
    end do
    do i = (it_made + 1), it_doiter
        it(i) = i
    end do

    open (1, file='z.dat', err=101)
    do i = 1, nz
        write (1, *, err=103) z(i)
    end do
    close (1)

    open (1, file='x.dat', err=101)
    do i = 1, nx
        write (1, *, err=103) x(i)
    end do
    close (1)

    open (1, file='k.dat', err=101)
    do i = 1, nk
        write (1, *, err=103) i
    end do
    close (1)

    !open(1, file = 'it.dat', err = 101)
    !do i=(it_made + 1),it_doiter
    ! write(1,*,err = 103) it(i)
    !enddo
    !close(1)

    return
101 stop 'error of file open.'
103 stop 'error of file writing.'
end subroutine calc_zxit

function a0_fn_stat() result(a0_res)
    use ic, only: nx, hx, iimp_x0, a0_peak, pi, imp_xsp, imp_x0, iimp_xend, in_type, lx, central_mirror, xcp, alfa!, coeff
    use fourier

    implicit none

    complex(c_double_complex), dimension(nx) :: a0_res, c
    real(c_double), dimension(nx) :: a0env
    integer i, icp, ix(nx)

    if (in_type == 1) then
        !initial conditions for a (one pulse in the middle)
        if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
        a0_res(nx/2 + 2:) = 0.0d0
        do i = iimp_x0, nx/2 + 1
            a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))
        end do
        a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
    elseif (in_type == 2) then
        !initial conditions for a (symmetric pulses at the edges)
        if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
        if (iimp_xend < nx) a0_res(iimp_xend + 1:nx) = 0.0d0
        do i = iimp_x0, iimp_xend
            a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))
        end do
        a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
    elseif (in_type == 3) then
        !initial conditions from harmonics
        c = 0
        !seed = (/2147483562, 2147483398/)
        !seed = (/3, 2/)
        !call random_seed(size = n)
        !if (n /= 2) stop 'error of random at a0_fn_stat'
        !call random_seed (put = seed)
        !call random_number(rc)
        !c(2:10:2) = cmplx(0.1 * rc(1:5), 0.0d0)

        c(2) = 0.1
        c(3) = 0.1
        c(4) = 0.05
        c(5) = 0.05

        !print *, size(c(2:10:2))
        !do i=1,size(rc)
        ! write(*,'(a,i2,a,f6.4,a,i2,a,f6.4,a,f6.4)') 'rc(', i, ') = ', rc(i) , ' c(', 2*i, ') = ', dreal(c(2*i)), ' ', dimag(c(2*i))
        !enddo

        a0_res = ifs(c)
        a0_res = cmplx(dreal(a0_res), 0.0d0)

        !open(1, file = 'test.dat')
        !do i=1,nx
        ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i)), abs(a0_res(i))
        !enddo
        !close(1)
        !stop
    elseif (in_type == 4) then
        !test initial conditions for a (one pulse in the middle)
        do i = 1, nx
            a0_res(i) = a0_peak*dsin(1*pi/lx*(i - 1)*hx)
        end do
    elseif (in_type == 5) then
        !specialist. initial conditions for a

        c = dcmplx(0)

        !c(2) = dcmplx(0.1)
        c(3) = dcmplx(0.1)
        !c(4) = dcmplx(0.1)
        !c(6) = dcmplx(0.1)
        !c(8) = dcmplx(0.1)
        !c(10) = dcmplx(0.1)
        !c(12) = dcmplx(0.1)
        !c(14) = dcmplx(0.1)

        a0_res = ifs(c)

        !open(1, file = 'test.dat')
        !do i=1,nx
        ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
        !enddo
        !close(1)
        !stop
    elseif (in_type == 6) then
        !specialist. initial conditions for a
        !initial conditions for a (symmetric pulses at the edges)

        if (central_mirror == .false.) then
            icp = xcp/hx + 1
            iimp_xend = 2*icp - 1
            ix = 0; ix = (/0:nx/2 - 1/)

            a0_res(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
            a0_res(nx/2 + 2:nx) = a0_res(nx/2:2:-1)
        else
            icp = xcp/hx + 1
            iimp_xend = 2*icp - 1
            ix = 0; ix = (/0:nx - 1/)

            a0_res = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
        end if

        !open(1, file = 'test.dat')
        !do i=1,nx
        ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
        ! !write(1, *) ix(i)
        !enddo
        !close(1)
        !stop
    elseif (in_type == 7) then
        !specialist. initial conditions for a
        !initial conditions for a (symmetric pulses at the edges)

        c = dcmplx(0)

        !c(2) = dcmplx(0.1)
        c(4) = dcmplx(0.1)
        c(6) = dcmplx(0.1)
        c(8) = dcmplx(0.1)
        c(10) = dcmplx(0.1)
!c(12) = dcmplx(0.1)
        !c(14) = dcmplx(0.1)

        if (central_mirror == .false.) then
            icp = xcp/hx + 1
            iimp_xend = 2*icp - 1
            ix = 0; ix = (/0:nx/2 - 1/)

            a0env(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
            a0env(nx/2 + 2:nx) = a0env(nx/2:2:-1)
        else
            ix = (/1:nx/) - 1

            a0env = dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
        end if

        a0_res = ifs(c)*a0env

        !open(1, file = 'test.dat')
        !do i=1,nx
        ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
        ! !write(1, *) ix(i)
        !enddo
        !close(1)
        !stop
    else
        print *, 'error: wrong in_type'
        pause
        stop
    end if

end function a0_fn_stat

function fk_fn(xe) result(fk_res)
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use ic, only: nk, pi, lx

    implicit none

    real(c_double) :: fk_res(nk), xe
    integer(c_int) n(nk)

    n = (/0:nk - 1/)

    fk_res = dsin(pi*n*xe/lx)

    !open(1, file = 'test.dat')
    !do i=1,nx
    ! write(1,'(i,e17.8)') i-1, fk_res(i)
    !enddo
    ! close (1)
    ! stop
end function fk_fn

! function g_fn () result (g_res)
! use ic
!
! implicit none
!
! real (c_double), dimension (nx) :: g_res
! integer i
!
! i = g_x0 / hx + 1
! ! i = g_x0 / hx
!
! if (central_mirror == .false.) then
! g_res (1: i) = 1.0d0
! g_res (nx-i+2: nx) = 1.0d0
! g_res (i+1: nx-i+1) = 0.0d0
! g_res = g_res * g_amp
! else
! g_res (1: i) = 0.0d0
! g_res (nx-i+2: nx) = 0.0d0
! g_res (i+1: nx-i+1) = 1.0d0
! g_res = g_res * g_amp
! endif
! end function g_fn

function g_fn() result(g_res)
    use ic

    implicit none

    real(c_double), dimension(nx) :: g_res
    integer(c_int) i

    ig0 = g_x0/hx + 1
    ig1 = g_x1/hx + 2

    if (central_mirror == .false.) then
        g_res = 1.0d0
        g_res(ig0:ig1) = 0.0d0
        soutm = -1.0
        smirr = 0.0
    else
        g_res = 0.0d0
        g_res(ig0:ig1) = 1.0d0
        smirr = -1.0
        soutm = 0.0
    end if

    do i = 1, nx
        if (g_res(i) > 0.0) then
            smirr = smirr + 1.0
        else
            soutm = soutm + 1.0
        end if
    end do

    g_res = g_res*g_amp

    ! print *, 'smirr =', smirr, 'soutm =', soutm, 'ig0 =', ig0, 'ig1 =', ig1
    ! stop
    ! print *, i0, nx-i0+2, i1
    ! stop
end function g_fn

function k_fn() result(k)
    use ic, only: nx
    use, intrinsic :: iso_c_binding

    implicit none

    complex(c_double_complex), dimension(2*nx) :: k
    complex(c_double_complex) :: im1 = (0.0d0, 1.0d0)
    integer nn

    nn = 2*nx

    k = im1*(/0:nn/2 - 1, -nn/2:-1/)
end function k_fn

function k2_fn() result(k2_res)
    use ic

    implicit none

    complex(c_double_complex), dimension(nk) :: k2_res
    integer i
    real(c_double) w

    !k**2
    do i = 1, nk
        w = pi*(i - 1)/lx
        !k2_res(i) = w * w - was
        !k2_res(i) = - w * w ! become
        k2_res(i) = -w*w/kk! become
    end do

    open (1, file='k2_n.dat')
    do i = 1, nk
        write (1, '(i,2e17.8)') i, k2_res(i)
    end do
    close (1)
end function k2_fn

function dn_fn() result(dn_res)
    use ic, only: nk, c_double_complex, c_double, lambda, lx, im1, pi

    implicit none

    complex(c_double_complex), dimension(nk) :: dn_res
    complex(c_double_complex) k
    real(c_double) tmp
    integer i

    k = 2.0d0*pi/lambda

    dn_res(1) = dcmplx(1)
    do i = 1, nk
        tmp = 1.0d0 - (i - 1)*(i - 1)/4.0d0*(lambda/lx)*(lambda/lx)
        dn_res(i) = dsqrt(tmp) - 1.0d0
    end do

    dn_res = k*dn_res

    !tmp = 1.0d0 - 1.0d0 / 4.0d0 * (lambda / lx) * (lambda / lx)
    !if (tmp .ge. 0.0d0) then
    !    dn1 = dsqrt(tmp)
    !else
    !    dn1 = im1 * dsqrt(dabs(tmp))
    !endif
    !
    !do i=1,nk
    !    tmp = 1.0d0 - (i * i) / 4.0d0 * (lambda / lx) * (lambda / lx)
    !    if (tmp .ge. 0.0d0) then
    !        dn_res(i) = dsqrt(tmp) - dn1
    !    else
    !        dn_res(i) = im1 * dsqrt(dabs(tmp)) - dn1
    !    endif
    !end do

    open (1, file='delta_n.dat')
    do i = 1, nk
        write (1, '(i,2e17.8)') i, dn_res(i)
    end do
    close (1)
end function dn_fn

subroutine allocate_arrays()
    use ic

    implicit none

    integer(c_int) err_alloc

    allocate (g(nx), a1(nx), a0(nx), ak1(nk), ak0(nk), jk1(nk), jk0(nk), atmp(nx), &
              th0(nth, 2), th1(nth, 2), dthdz(nth, 2), fk1(nk), fk2(nk), rhs0(nth, 2), z(nz), x(nx), k2(nk), &
              a_amp_z0(nx), a_amp_zl(nx), aktmp(nx), akzl(nk), akz0(nk), a0z0(nx), a0z0cut(nx), &
              a_spec_amp_z0(nx), a_spec_amp_zl(nx), it(it_todo), &
              ex(nk), k(2*nk), dlt(nk), sum_abs2_a_plus_by_z(nx), sum_abs2_a_plus_by_z_k(nk), tmp(nx), &
              !theta(nth, 2, nz), &
              stat=err_alloc)

    if (err_alloc /= 0) then
        pause "allocation error"
        stop
    end if
end subroutine allocate_arrays

subroutine deallocate_arrays()
    use ic

    implicit none

    integer(c_int) err_dealloc

    deallocate (g, a1, a0, ak1, ak0, jk1, jk0, atmp, &
                th0, th1, dthdz, fk1, fk2, rhs0, z, x, k2, &
                a_amp_z0, a_amp_zl, aktmp, akzl, akz0, a0z0cut, &
                a_spec_amp_z0, a_spec_amp_zl, it, &
                ex, k, dlt, sum_abs2_a_plus_by_z, sum_abs2_a_plus_by_z_k, tmp, &
                !theta, &
                stat=err_dealloc)

    if (err_dealloc /= 0) stop "deallocation error"
end subroutine deallocate_arrays

subroutine read_param() bind(c, name='read_param')
    use ic

    implicit none

    open (unit=1, file='input_fortran.dat', status='old', err=101)
    read (unit=1, nml=param, err=102)
    close (unit=1)

    write (*, nml=param)

    return
101 print *, 'error of file open'; pause; stop
102 print *, 'error of reading file "input_fortran.dat"'; pause; stop
end subroutine read_param

!subroutine d(dadx)
!    use ic, only : k, nx, lx, pi, hx
!    use fourier
!    use, intrinsic :: iso_c_binding
!
!    implicit none
!
!    complex(c_double_complex), dimension(nx) :: dadx
!    complex(c_double_complex), dimension(2*nx) :: v
!
!    v(1) = (0.0d0,0.0d0)
!    v(2:nx) = dadx(2:nx)
!    v(nx+1) = (0.0d0,0.0d0)
!    v(2*nx:nx+2:-1) = -dadx(2:nx)
!
!    call fft(v)
!    v = 2.0 * pi / (2.0d0 * lx) * k * v
!    call ifft(v)
!    dadx = v(1:nx)
!
!end subroutine d

!function check_alloc()
!use, intrinsic :: iso_c_binding, only : c_bool
!    use ic, only : a1, a0, ak1, ak0, atmp, jk1, jk0, dadx, k, ex, &
!        th0, th1, dthdz, fk1, fk2, rhs0, z, x, k2, &
!        a_amp_z0, a_amp_zl, a_spec_amp_z0, a_spec_amp_zl, g, sum_abs2_dadx_by_x
!
!    implicit none
!
!    logical(c_bool) check_alloc
!
!    !print *, allocated(a0), size(a0)
!
!    if (.not. allocated(a0)) then
!        check_alloc = .false.
!        print *, 'arrays were not allocated!'
!        pause
!        stop
!    endif
!
!    check_alloc = .true.
!end function check_alloc

subroutine makea(atmp, ak)
    use ic, only: nk, nx
    use fourier
    use, intrinsic :: iso_c_binding, only: c_double_complex, c_int

    implicit none

    complex(c_double_complex), dimension(:), intent(inout) :: atmp
    complex(c_double_complex), dimension(:), intent(in) :: ak
    integer(c_int) n1, n2

    n1 = size(ak)
    n2 = size(atmp)

    if (n1 .ne. nk .or. n2 .ne. nx) then
        print *, 'error in "makea"'
        pause
        stop
    end if

    atmp = dcmplx(0.0d0)
    atmp(1:nk) = ak

    call isint(atmp)
end subroutine makea

subroutine makeak(ak, atmp)
    use ic, only: nk, nx, aktmp
    use fourier
    use, intrinsic :: iso_c_binding, only: c_double_complex, c_int

    implicit none

    complex(c_double_complex), dimension(:), intent(inout) :: ak
    complex(c_double_complex), dimension(:), intent(in) :: atmp
    integer(c_int) n1, n2

    n1 = size(ak)
    n2 = size(atmp)

    if (n1 .ne. nk .or. n2 .ne. nx) then
        print *, 'error in "makeak"'
        pause
        stop
    end if

    aktmp = atmp

    call sint(aktmp)
    ak = aktmp(1:nk)
end subroutine makeak

subroutine make_a0z0_ak1_atmp(a0z0, az0cut, ak0, ak)
    use ic, only: nk, nx, r0, g
    use fourier
    use, intrinsic :: iso_c_binding, only: c_double_complex, c_int

    implicit none

    interface
        subroutine makea(atmp, ak)
            use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
            complex(c_double_complex), dimension(:), intent(inout) :: atmp
            complex(c_double_complex), dimension(:), intent(in) :: ak
        end subroutine makea
        subroutine makeak(ak, atmp)
            use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
            complex(c_double_complex), dimension(:), intent(inout) :: ak
            complex(c_double_complex), dimension(:), intent(in) :: atmp
        end subroutine makeak
    end interface

    complex(c_double_complex), dimension(:), intent(inout) :: a0z0, ak0, az0cut
    complex(c_double_complex), dimension(:), intent(in) :: ak
    integer(c_int) n1, n2, n3, n4

    n1 = size(ak)
    n2 = size(az0cut)
    n3 = size(a0z0)
    n4 = size(ak0)

    if (n1 .ne. nk .or. n2 .ne. nx .or. n3 .ne. nx .or. n4 .ne. nk) then
        print *, 'error in "makea"'
        pause
        stop
    end if

    call makea(a0z0, ak) !before cutting

    az0cut = a0z0*r0*g !mirror reflection in z=0 and cut

    call sint(az0cut)

    ak0 = az0cut(1:nk)

    call makea(az0cut, ak0)
end subroutine make_a0z0_ak1_atmp
