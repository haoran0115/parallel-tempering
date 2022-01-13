program parallel_tempering
    implicit none
    ! parameters
    integer, parameter :: total_replica = 4
    integer, parameter :: n_iter        = 1000000
    integer, parameter :: n_record      = 1
    integer, parameter :: n_exchage     = 50
    real(8), parameter :: pi            = 3.14159265
    real(8), parameter :: kB            = 1.0 
    real(8), parameter :: eps           = 1.0 
    real(8), parameter :: stepsize      = 1.0
    real(8), parameter :: border_l      = 0.0
    real(8), parameter :: border_r      = 80.0
    real(8), parameter, dimension(total_replica) :: temp = [0.001, 0.5, 1.0, 2.0]

    ! variables
    real(8), dimension(total_replica) :: x_coords_o, x_coords_n ! coordinates of replica
    real(8), dimension(total_replica) :: x_temp = temp          ! temperature
    real(8), dimension(total_replica) :: x_disp                 ! displacement of replica
    real(8), dimension(total_replica) :: x_ratio                ! acceptance of replica
    real(8), dimension(total_replica) :: x_energy               ! H(x)
    logical, dimension(total_replica) :: x_ifmove               ! true/false, whether move
    integer, dimension(2)             :: x_choice               ! choice for Hamitonian swapping
    integer                           :: total_intervals        ! number of intervals in hist
    ! tmp variables
    integer :: i, j, k, l
    real(8) :: r1, xk, xl, exchange_ratio

    ! open .dat files
    open(1, file="coords.dat")
    open(2, file="info.dat")
    open(3, file="temp.dat")
    ! write informations about border and total intervals
    total_intervals = (border_r - border_l) / stepsize
    write(2, "(2F10.2, I10)") border_l, border_r, total_intervals
    ! write informations about temperature
    write(3, *) temp

    call random_seed()
    do i = 1, n_iter
        x_disp     = random_move()
        x_coords_n = x_coords_o + x_disp
        x_ratio    = fr(x_coords_o, x_coords_n, x_temp)
        x_ifmove   = accept_move(x_ratio)
        do j = 1, total_replica
            if (x_ifmove(j) .and. x_coords_n(j) .ge. border_l .and. &
              & x_coords_n(j) .le. border_r) then
                x_coords_o(j) = x_coords_n(j)
            end if
        end do

        ! perform Hamitonian swapping every n_exchage
        if (mod(i, n_exchage) .eq. 0) then
            call two_choice(x_choice)
            x_energy = u(x_coords_o)
            k = x_choice(1)
            l = x_choice(2)
            exchange_ratio = exp((x_energy(k)-x_energy(l))* &
              & (1/x_temp(k)-1/x_temp(l))/kB)
            call random_number(r1)
            if (r1 .le. exchange_ratio) then
                xk = x_coords_o(k)
                xl = x_coords_o(l)
                x_coords_o(k) = xl
                x_coords_o(l) = xk
            end if
        end if

        ! do record every n_record
        if (mod(i, n_record) .eq. 0) then
            write(1, "(4F10.2)") x_coords_o
        end if
    end do

    ! close .dat files
    close(1)
    close(2)
    close(3)

    contains

    ! potential energy function
    function u(x) result(r)
        real(8), dimension(total_replica), intent(in) :: x
        real(8), dimension(total_replica) :: r
        ! r = x * eps
        r = cos(x*pi/10)*5 + 1
    end function 
    
    ! metropolis function 
    function accept_move(N_ratio) result(r)
        logical, dimension(total_replica) :: r
        real(8), dimension(total_replica), intent(in) :: N_ratio
        real(8), dimension(total_replica) :: acc ! random number
        call random_number(acc(:))
        r = acc .le. N_ratio
    end function
    
    ! boltzmann ratio to change two replicas
    function fr(x1, x2, t) result(r)
        real(8), dimension(total_replica), intent(in) :: x1, x2, t
        real(8), dimension(total_replica) :: r
        r = exp(-1/(kB*t) * (u(x2) - u(x1)))
    end function

    ! random move
    function random_move() result(disp)
        real(8), dimension(total_replica) :: disp
        real(8), dimension(total_replica) :: rand1, rand2
        integer :: idx
        call random_number(rand1)
        call random_number(rand2)
        rand1 = rand1 - 0.5
        rand2 = rand2 - 0.5
        do idx = 1, total_replica
            if (rand1(idx) .le. 0.0) then
                disp(idx) = sign(stepsize, rand2(idx))
            end if
        end do
    end function

    ! generate two index (1, total_replica)
    subroutine two_choice(choice)
        integer, dimension(2) :: choice
        real(8), dimension(2) :: choice_r
        choice = [1, 1]
        do while (choice(1) .eq. choice(2))
            call random_number(choice_r(:))
            choice = choice_r * total_replica + 1.0
        end do
    end subroutine


end program


