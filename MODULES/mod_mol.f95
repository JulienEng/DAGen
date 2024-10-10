module mod_atom
        implicit none
        public

        type,public::atom_class
                integer,public::indx
                integer,public::group
                double precision,public::mass
                double precision,dimension(3),public::coord
                double precision::dist_com
                integer::sym_equi
                character(len=2),public::label
                integer,public::z_number
         end type

        type,public::mol_class
                character(len=666)                           :: file_xyz
                
                integer,public                               :: at
                integer,public                               :: sym
                
                type(atom_class),dimension(:),allocatable    :: atoms
                type(atom_class),dimension(:),allocatable    :: atoms_re

                double precision                             :: mass_tot
                double precision                             :: mass_red
                double precision,dimension(3)                :: cntr_mass
                double precision,dimension(3)                :: inertia_moment
                double precision,dimension(3,3)              :: inertia_mat
                double precision,dimension(3,3)              :: inertia_rot
                
                integer                                      :: n_groups
                integer                                      :: size_groups
                integer,dimension(:,:),allocatable           :: sym_groups
                integer                                      :: max_size

                double precision,dimension(:,:),allocatable  :: hessian
                double precision,dimension(:),allocatable    :: gradient
                double precision,dimension(:),allocatable    :: kappa_g
                double precision,dimension(:),allocatable    :: kappa_d
                double precision,dimension(:,:),allocatable  :: dnc
                double precision,dimension(:,:),allocatable  :: nm
                double precision,dimension(:),allocatable    :: xyz
                double precision,dimension(:),allocatable    :: rmsa
                double precision                             :: gap
                double precision                             :: energy
         contains
                 procedure,public::print_out  =>  mol_print_out
                 procedure,public::read_xyz  =>  mol_read_xyz
                 procedure,public::init_xyz => mol_init_xyz
                 procedure,public::get_sym => mol_get_sym
         end type mol_class

contains
        subroutine mol_read_xyz(this)
                class(mol_class),intent(inout)::this
                integer::i,j
                integer::k,l
                double precision::normvec

                open(98,file=trim(adjustl(this%file_xyz)))
                read(98,*) this%at
                allocate(this%atoms(this%at))
                read(98,*)

                do i=1,this%at
                read(98,*) this%atoms(i)%label, (this%atoms(i)%coord(j),j=1,3)
                end do
        end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        subroutine mol_init_xyz(this)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!!  Here we prepare the molecule and atom for the subsequent analysis.          !!!
!!!  Mass is assigned to each atom and their molecule is centred and rotated     !!!
!!! to aligne its Inertia axis. Inertia moments are calculated as well.          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                
                class(mol_class),intent(inout)::this
                integer::i,j
                integer::k,l
                double precision::normvec                                        !!! see normvec function in math.f95
                double precision,dimension(118)::dtbase_mass                     !!! database array with all the masses
                character(len=2),dimension(118)::dtbase_labels                   !!! matching database array with all the elements' name.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                dtbase_labels(1)='H'       ;  dtbase_mass(1)=1.0079
                dtbase_labels(2)='He'      ;  dtbase_mass(2)=4.0026
                dtbase_labels(3)='Li'      ;  dtbase_mass(3)=6.941
                dtbase_labels(4)='Be'      ;  dtbase_mass(4)=9.0122
                dtbase_labels(5)='B'       ;  dtbase_mass(5)=10.811
                dtbase_labels(6)='C'       ;  dtbase_mass(6)=12.0107
                dtbase_labels(7)='N'       ;  dtbase_mass(7)=14.0067
                dtbase_labels(8)='O'       ;  dtbase_mass(8)=15.9994
                dtbase_labels(9)='F'       ;  dtbase_mass(9)=18.9984
                dtbase_labels(10)='Ne'     ;  dtbase_mass(10)=20.1797
                dtbase_labels(11)='Na'     ;  dtbase_mass(11)=22.9897
                dtbase_labels(12)='Mg'     ;  dtbase_mass(12)=24.305
                dtbase_labels(13)='Al'     ;  dtbase_mass(13)=26.9815
                dtbase_labels(14)='Si'     ;  dtbase_mass(14)=28.0855
                dtbase_labels(15)='P'      ;  dtbase_mass(15)=30.9738
                dtbase_labels(16)='S'      ;  dtbase_mass(16)=32.065
                dtbase_labels(17)='Cl'     ;  dtbase_mass(17)=35.453
                dtbase_labels(18)='Ar'     ;  dtbase_mass(18)=39.0983
                dtbase_labels(19)='K'      ;  dtbase_mass(19)=39.948
                dtbase_labels(20)='Ca'     ;  dtbase_mass(20)=40.078
                dtbase_labels(21)='Sc'     ;  dtbase_mass(21)=44.9559
                dtbase_labels(22)='Ti'     ;  dtbase_mass(22)=47.867
                dtbase_labels(23)='V'      ;  dtbase_mass(23)=50.9415
                dtbase_labels(24)='Cr'     ;  dtbase_mass(24)=51.9961
                dtbase_labels(25)='Mn'     ;  dtbase_mass(25)=54.938
                dtbase_labels(26)='Fe'     ;  dtbase_mass(26)=55.845
                dtbase_labels(27)='Co'     ;  dtbase_mass(27)=58.6934
                dtbase_labels(28)='Ni'     ;  dtbase_mass(28)=58.9332
                dtbase_labels(29)='Cu'     ;  dtbase_mass(29)=63.546
                dtbase_labels(30)='Zn'     ;  dtbase_mass(30)=65.39
                dtbase_labels(31)='Ga'     ;  dtbase_mass(31)=69.723
                dtbase_labels(32)='Ge'     ;  dtbase_mass(32)=72.64
                dtbase_labels(33)='As'     ;  dtbase_mass(33)=74.9216
                dtbase_labels(34)='Se'     ;  dtbase_mass(34)=78.96
                dtbase_labels(35)='Br'     ;  dtbase_mass(35)=79.904
                dtbase_labels(36)='Kr'     ;  dtbase_mass(36)=83.8
                dtbase_labels(37)='Rb'     ;  dtbase_mass(37)=85.4678
                dtbase_labels(38)='Sr'     ;  dtbase_mass(38)=87.62
                dtbase_labels(39)='Y'      ;  dtbase_mass(39)=88.9059
                dtbase_labels(40)='Zr'     ;  dtbase_mass(40)=91.224
                dtbase_labels(41)='Nb'     ;  dtbase_mass(41)=92.9064
                dtbase_labels(42)='Mo'     ;  dtbase_mass(42)=95.94
                dtbase_labels(43)='Tc'     ;  dtbase_mass(43)=98
                dtbase_labels(44)='Ru'     ;  dtbase_mass(44)=101.07
                dtbase_labels(45)='Rh'     ;  dtbase_mass(45)=102.9055
                dtbase_labels(46)='Pd'     ;  dtbase_mass(46)=106.42
                dtbase_labels(47)='Ag'     ;  dtbase_mass(47)=107.8682
                dtbase_labels(48)='Cd'     ;  dtbase_mass(48)=112.411
                dtbase_labels(49)='In'     ;  dtbase_mass(49)=114.818
                dtbase_labels(50)='Sn'     ;  dtbase_mass(50)=118.71
                dtbase_labels(51)='Sb'     ;  dtbase_mass(51)=121.76
                dtbase_labels(52)='Te'     ;  dtbase_mass(52)=126.9045
                dtbase_labels(53)='I'      ;  dtbase_mass(53)=127.6
                dtbase_labels(54)='Xe'     ;  dtbase_mass(54)=131.293
                dtbase_labels(55)='Cs'     ;  dtbase_mass(55)=132.9055
                dtbase_labels(56)='Ba'     ;  dtbase_mass(56)=137.327
                dtbase_labels(57)='La'     ;  dtbase_mass(57)=138.9055
                dtbase_labels(58)='Ce'     ;  dtbase_mass(58)=140.12
                dtbase_labels(59)='Pr'     ;  dtbase_mass(59)=140.91
                dtbase_labels(60)='Nd'     ;  dtbase_mass(60)=144.24
                dtbase_labels(61)='Pm'     ;  dtbase_mass(61)=145.00
                dtbase_labels(62)='Sm'     ;  dtbase_mass(62)=150.36
                dtbase_labels(63)='Eu'     ;  dtbase_mass(63)=151.96
                dtbase_labels(64)='Gd'     ;  dtbase_mass(64)=157.25
                dtbase_labels(65)='Tb'     ;  dtbase_mass(65)=158.93
                dtbase_labels(66)='Dy'     ;  dtbase_mass(66)=162.50
                dtbase_labels(67)='Ho'     ;  dtbase_mass(67)=164.93
                dtbase_labels(68)='Er'     ;  dtbase_mass(68)=167.26
                dtbase_labels(69)='Tm'     ;  dtbase_mass(69)=168.93
                dtbase_labels(70)='Yb'     ;  dtbase_mass(70)=173.05
                dtbase_labels(71)='Lu'     ;  dtbase_mass(71)=174.97
                dtbase_labels(72)='Hf'     ;  dtbase_mass(72)=178.49
                dtbase_labels(73)='Ta'     ;  dtbase_mass(73)=180.95
                dtbase_labels(74)='W'      ;  dtbase_mass(74)=183.84
                dtbase_labels(75)='Re'     ;  dtbase_mass(75)=186.027
                dtbase_labels(76)='Os'     ;  dtbase_mass(76)=190.23
                dtbase_labels(77)='Ir'     ;  dtbase_mass(77)=192.22
                dtbase_labels(78)='Pt'     ;  dtbase_mass(78)=195.08
                dtbase_labels(79)='Au'     ;  dtbase_mass(79)=196.97
                dtbase_labels(80)='Hg'     ;  dtbase_mass(80)=200.59
                dtbase_labels(81)='Tl'     ;  dtbase_mass(81)=204.38
                dtbase_labels(82)='Pb'     ;  dtbase_mass(82)=207.20
                dtbase_labels(83)='Bi'     ;  dtbase_mass(83)=208.98
                dtbase_labels(84)='Po'     ;  dtbase_mass(84)=209
                dtbase_labels(85)='At'     ;  dtbase_mass(85)=210
                dtbase_labels(86)='Rn'     ;  dtbase_mass(86)=222
                dtbase_labels(87)='Fr'     ;  dtbase_mass(87)=223
                dtbase_labels(88)='Ra'     ;  dtbase_mass(88)=226
                dtbase_labels(89)='Ac'     ;  dtbase_mass(89)=227
                dtbase_labels(90)='Th'     ;  dtbase_mass(90)=232.04
                dtbase_labels(91)='Pa'     ;  dtbase_mass(91)=231.04
                dtbase_labels(92)='U'      ;  dtbase_mass(92)=238.03
                dtbase_labels(93)='Np'     ;  dtbase_mass(93)=237
                dtbase_labels(94)='Pu'     ;  dtbase_mass(94)=244
                dtbase_labels(95)='Am'     ;  dtbase_mass(95)=243
                dtbase_labels(96)='Cm'     ;  dtbase_mass(96)=247
                dtbase_labels(97)='Bk'     ;  dtbase_mass(97)=247
                dtbase_labels(98)='Cf'     ;  dtbase_mass(98)=251
                dtbase_labels(99)='Es'     ;  dtbase_mass(99)=252
                dtbase_labels(100)='Fm'    ;  dtbase_mass(100)=257
                dtbase_labels(101)='Md'    ;  dtbase_mass(101)=258
                dtbase_labels(102)='No'    ;  dtbase_mass(102)=259
                dtbase_labels(103)='Lr'    ;  dtbase_mass(103)=262
                dtbase_labels(104)='Rf'    ;  dtbase_mass(104)=267
                dtbase_labels(105)='Db'    ;  dtbase_mass(105)=268
                dtbase_labels(106)='Sg'    ;  dtbase_mass(106)=269
                dtbase_labels(107)='Bh'    ;  dtbase_mass(107)=270
                dtbase_labels(108)='Hs'    ;  dtbase_mass(108)=269
                dtbase_labels(109)='Mt'    ;  dtbase_mass(109)=277
                dtbase_labels(110)='Ds'    ;  dtbase_mass(110)=281
                dtbase_labels(111)='Rg'    ;  dtbase_mass(111)=282
                dtbase_labels(112)='Cn'    ;  dtbase_mass(112)=285
                dtbase_labels(113)='Nh'    ;  dtbase_mass(113)=286
                dtbase_labels(114)='Fl'    ;  dtbase_mass(114)=290
                dtbase_labels(115)='Mc'    ;  dtbase_mass(115)=290
                dtbase_labels(116)='Lv'    ;  dtbase_mass(116)=293
                dtbase_labels(117)='Ts'    ;  dtbase_mass(117)=294
                dtbase_labels(118)='Og'    ;  dtbase_mass(118)=294


                allocate(this%atoms_re(this%at))
                ! Initialisation
                this%cntr_mass=0.0d0
                this%mass_tot=0.0d0

                do i=1,this%at ! loop over the number of atoms.

                        if (this%atoms(i)%label == '') then
                               ! if (this%atoms(i)%z_number /= '') then
                                        this%atoms(i)%label=dtbase_labels(this%atoms(i)%z_number)
                               ! end if
                        else

                                k=FINDLOC(dtbase_labels,this%atoms(i)%label,1) ! Find the element of atom i in the database.
                                this%atoms(i)%z_number=k
                                this%atoms(i)%mass=dtbase_mass(k)              ! Assign the mass
                                this%mass_tot=this%mass_tot+this%atoms(i)%mass ! comput the total mass of the molecule
                                do j=1,3     
                                  this%cntr_mass(j)=this%cntr_mass(j)+this%atoms(i)%mass*this%atoms(i)%coord(j) ! Compute the centre of mass
                                end do
                        end if
                end do
                this%cntr_mass=this%cntr_mass/this%mass_tot 
                !! Centre the molecule around the centre of mass.                
                do i=1,this%at
                        this%atoms_re(i)%coord(:)=this%atoms(i)%coord(:)-this%cntr_mass(:)
                        this%atoms_re(i)%dist_com=normvec(this%atoms_re(i)%coord(:))
                        this%atoms(i)%dist_com=normvec(this%atoms(i)%coord(:)-this%cntr_mass(:))
                end do

this%inertia_mat=0.0d0  !! Build the inertia matrice
                 do k=1,3
                    do l=1,3
                        if (l.eq.k) then
                            do j=1,this%at
this%inertia_mat(k,l)=this%inertia_mat(k,l)+this%atoms(j)%mass*(normvec(this%atoms_re(j)%coord(:))**2-&
&this%atoms_re(j)%coord(k)*this%atoms_re(j)%coord(l))
                             end do
                        else
                             do j=1,this%at
this%inertia_mat(k,l)=this%inertia_mat(k,l)-this%atoms(j)%mass*this%atoms_re(j)%coord(k)*this%atoms_re(j)%coord(l)
                             end do
                        end if
                    end do
                end do
                
                call diag(this%inertia_mat,3,this%inertia_rot,this%inertia_moment) ! Diagonalisation
                do i=1,this%at
                this%atoms_re(i)%coord(:)=MATMUL(TRANSPOSE(this%inertia_rot),this%atoms_re(i)%coord( :))
                end do


                 CLOSE(98)
        end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine mol_get_sym(this)                                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                class(mol_class),intent(inout)::this                             
                integer::counter,i,j,k
                integer,dimension(2)::non_null
                double precision,dimension(3,3)::matop
                double precision::angle
                double precision::scalar,normvec
                double precision,dimension(3)::vec,vec_a,vec_b
                integer::equi
                logical::checker
                integer,dimension(3)::c2_axis
                double precision,parameter::pi=4.0d0*atan(1.0d0)

                this%atoms(:)%sym_equi=0
                this%size_groups=0
                counter=1
                do i=1,this%at
                        if (this%atoms(i)%sym_equi.eq.0) then
                        where (ABS(this%atoms(:)%dist_com) < ABS(this%atoms(i)%dist_com) + 0.001 &
                        &.AND. ABS(this%atoms(:)%dist_com) > ABS(this%atoms(i)%dist_com) - 0.001 )
                        this%atoms(:)%sym_equi=counter
                        end where
                        counter=MAXVAL(this%atoms(:)%sym_equi)+1
                        end if
                        this%size_groups=MAX(this%size_groups,count(this%atoms(:)%sym_equi == counter-1))
                end do
                this%n_groups=counter-1
                
                allocate(this%sym_groups(this%size_groups,this%n_groups))
                this%sym_groups=0
                do i=1,this%n_groups
                        this%sym_groups(:,i)=pack([(j, j=1, this%at)] , this%atoms(:)%sym_equi == i)
                end do

                !!! Check if equivalent atoms are not in different planes
                               call chk_plane(this,checker,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!! vv thresh
                do i=1,2
                non_null(i)=count( this%inertia_moment(:) <= this%inertia_moment(i)+0.01 &
                &.AND. this%inertia_moment(:) >= this%inertia_moment(i)-0.01)
                end do

                if (MAXVAL(non_null).eq.3) then
                        write(*,*) 'Molecule is Cubic !'
                        call get_inv(matop)
                        call get_equivalence(this,matop,checker)
                        if (checker.eqv..TRUE.) then
                              write(*,*) 'Oh or Ih'
                              angle=2*pi/4
                              c2_axis=0
                              do i=1,3
                                      vec=0
                                      vec(i)=i
                                      call get_rot(matop,vec,angle)
                                      call get_equivalence(this,matop,checker)
                                      if (checker.eqv..TRUE.) then
                                            c2_axis(i)=2
                                      else
                                           c2_axis(i)=0
                                      end if       
                              end do
                              if ( count(c2_axis(:) == 0 ) < 3 ) then
                                write(*,*) 'Oh'
                              else
                                write(*,*) 'Ih'
                              end if
                        else
                                write(*,*) 'Td'
                        end if
                else if (MAXVAL(non_null) == 2 .AND. MINVAL(this%inertia_moment) <= 0.01 ) then
                        write(*,*) 'Molecule is Linear'
                        call get_inv(matop)
                        call get_equivalence(this,matop,checker)
                        if (checker.eqv..TRUE.) then
                                write(*,*) 'D_{\infty h}' 
                        else
                                write(*,*) 'C_{\infty v}'
                        end if
                        !!!!!!!!!!!!  S Y M M E T R I C   T O P  !!!!!!!!!!!!!!!!!!!
                else if (MAXVAL(non_null) == 2 .AND. MINVAL(this%inertia_moment) > 0.01 ) then
                        write(*,*) 'Molecule is a Symmetric top'
                        ! Need to check for planes (see above) and get larger group.
                        k=this%size_groups
                        angle=2*pi/k
                        c2_axis=0

                                call vecto(this%atoms_re(this%sym_groups(1,this%max_size))%coord(:)-&
                                         &this%atoms_re(this%sym_groups(2,this%max_size))%coord(:),&
                                         &this%atoms_re(this%sym_groups(2,this%max_size))%coord(:)-&
                                         &this%atoms_re(this%sym_groups(3,this%max_size))%coord(:),vec)

                                call get_rot(matop,vec,angle)
                                call get_equivalence(this,matop,checker)
                                do i=1,3 ; if (ABS(vec(i)).lt.1E-4) then ; vec(i)=0.0d0 ; end if ; end do

                                if (checker.eqv..TRUE.) then        !!! If Cn axis then
                                        Write(*,'(a,i0)') 'Rotation of Order: C', k
                                        !! Check for C2 along bissector   <<< This is not correct.
 
                                        !! Need to select the group with the larges number of equivalent atoms
                                        !! Checking for ONE C2 axis perpendicular to the Cn axis. If one exist, there is n of them.

                                        vec_a(:)=(this%atoms_re(this%sym_groups(1,this%max_size))%coord(:)-&
                                        &this%atoms_re(this%sym_groups(2,this%max_size))%coord(:))/2
                                        do i=1,3 ; if (ABS(vec_a(i)).lt.1E-4) then ; vec_a(i)=0.0d0 ; end if ; end do
                                         angle=2*pi/2

                                        call get_rot(matop,vec_a,angle)
                                        call get_equivalence(this,matop,checker)
 
                                                if (checker.eqv..TRUE.) then !!! C2 perpendicular axis -- > D
                                                        write(*,'(a,i0)') ' C2 perpendicular to C',k
                                                         call get_sigma(matop,vec) 
                                                         call get_equivalence(this,matop,checker)
                                                                 if (checker.eqv..TRUE.) then !! Sigma h plane  -- > h
                                                                        write(*,'(a,i0)') 'Sigma plane perpendicular to C',k
                                                                        call vecto(vec_a,vec,vec_b)
                                                                        call get_sigma(matop,vec_b)   
                                                                        call get_equivalence(this,matop,checker)
                                                                        if (checker.eqv..TRUE.) then !! Sigma d plane 
                                                                              write(*,'(a,i0)') 'Sigma planes along C',k
                                                                              write(*,'(a,i0,a)') 'Molecule is D',k,'h'
                                                                     !   else
                                                                      !        write(*,*) 'Molecule is D',k,'h'
                                                                        end if
                                                                else                         !! No Sigma h plane -- > Dx(d/?)
                                                                        write(*,'(a,i0)') 'No Sigma plane perpendicular to C',k
                                                                        call vecto(vec_a,vec,vec_b)
                                                                        call get_sigma(matop,vec_b)   
                                                                        call get_equivalence(this,matop,checker)     
                                                                        if (checker.eqv..TRUE.) then  !! Sigma v >> Dxd
                                                                                write(*,'(a,i0)') 'Sigma planes along C',k
                                                                                write(*,'(a,i0,a)') 'Molecule is D',k,'d'
                                                                        else  !!>> Dn
                                                                                write(*,'(a,i0)') 'No Sigma planes along C',k
                                                                                write(*,'(a,i0)') 'Molecule is D',k

                                                                        end if 
                                                                 end if

                                                else                         !! No C2 perpendicular axis.  -- > C
                                                        c2_axis(i)=-1       
                                                        vec=0
                                                        vec(i)=1
                                                        call get_sigma(matop,vec)
                                                        call get_equivalence(this,matop,checker)
                                                                if (checker.eqv..TRUE.) then !! Sigma h plane
                                                                        write(*,'(a,i0,a)') 'Molecule is C',k,'h'
                                                                else                         !! No Sigma h plane
                                                        vec_a=(this%atoms_re(this%sym_groups(1,this%max_size))%coord(:)-&
                                                        &this%atoms_re(this%sym_groups(2,this%max_size))%coord(:))/2
                                                                        vec_b=0
                                                                        vec_b(i)=1
                                                                        call vecto(vec_a,vec_b,vec)
                                                                        call get_sigma(matop,vec)   
                                                                        call get_equivalence(this,matop,checker)     
                                                                        if (checker.eqv..TRUE.) then  !! Sigma v
                                                                                write(*,'(a,i0,a)') 'Molecule is C',k,'v'
                                                                        else
                                                                                angle=pi/k
                                                                                call get_sn(matop,vec,angle)
                                                                                call get_equivalence(this,matop,checker) 
                                                                                if (checker.eqv..TRUE.) then  !! Sigma v
                                                                                        write(*,'(a,i0)') 'Molecule is S',k*2
                                                                                else
                                                                                        write(*,'(a,i0)') 'Molecule is C',k
                                                                                end if
                                                                        end if 
                                                                end if
                                                end if

                                else   !! If no Cn axis

                                        !! check for a sigma h
                                        c2_axis(i)=0
                                        call get_sigma(matop,vec)
                                        call get_equivalence(this,matop,checker)
                                        if (checker.eqv..TRUE.) then
                                                !c2_axis(i)=-1
                                                write(*,'(a,i0,a)') 'Molecule is Cs'
                                        else
                                                call get_inv(matop)
                                                call get_equivalence(this,matop,checker)
                                                if (checker.eqv..TRUE.) then
                                                        write(*,'(a,i0,a)') 'Molecule is Ci'
                                                else
                                                        write(*,'(a,i0,a)') 'Molecule is C1'
                                                end if
                                        end if
                                end if       
                        !end do
                        !!!!!!!!!!!! A S Y M M E T R I C   T O P  !!!!!!!!!!!!!!!!!!!


                else if (MAXVAL(non_null) == 1 ) then
                        write(*,*) 'Molecule is an asymmetric top'
                        k=2
                        angle=2*pi/k
                        c2_axis=0
                        do i=1,3
                                vec=0
                                vec(i)=i
                                call get_rot(matop,vec,angle)
                                call get_equivalence(this,matop,checker)
                                if (checker.eqv..TRUE.) then
                                        c2_axis(i)=1
                                        call get_sigma(matop,vec)
                                        if (checker.eqv..TRUE.) then
                                                c2_axis(i)=2
                                        end if
                                else
                                        c2_axis(i)=0
                                        call get_sigma(matop,vec)
                                        if (checker.eqv..TRUE.) then
                                                c2_axis(i)=-1
                                        end if
                                end if       
                        end do
                        if (count(c2_axis == 0 ) == 3 ) then ; 
                                call get_inv(matop)
                                call get_equivalence(this,matop,checker)
                                if (checker.eqv..TRUE.) then
                                  write(*,*) 'Molecule is Ci'
                                else
                                  write(*,*) 'Molecule is C1'
                                end if
                        else if (count(c2_axis ==-1) == 1 ) then ; write(*,*) 'Molecule is Cs'
                        else if (count(c2_axis == 2) == 1 ) then ; write(*,*) 'Molecule is C2h'
                        else if (count(c2_axis ==-1) == 2 & 
                        & .AND.  count(c2_axis == 1) == 1 ) then ; write(*,*) 'Molecule is C2v'
                        else if (count(c2_axis == 1) == 3 ) then ; write(*,*) 'Molecule is D2'
                        else if (count(c2_axis == 1) == 1 &
                        & .AND.  count(c2_axis == 2) == 2 ) then ; write(*,*) 'Molecule is D2d'
                        else if (count(c2_axis == 2) == 3 ) then ; write(*,*) 'Molecule is D2h'
                        end if
                end if



        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        SUBROUTINE chk_plane(this,checker,k)                                     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!!  In rare case equivalency between atoms can not be accurately defined        !!!
!!!  based on their distance to the centre of mass alone. "n" equivalent atoms   !!!
!!! can be positioned in two different planes of n/2 atoms each.                 !!!
!!! This subroutines checks for cases like this and modify the this%sym_groups   !!!  
!!! Attention needs to be put to identifying the proper plane ! e.g. CH3-CH3     !!!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                
                implicit none                                                    !!!
                class(mol_class),intent(inout)::this                             !!!
                integer::n,i,j,k                                                 !!!
                integer::tmp                                                     !!! integer storage
                integer,dimension(:,:),allocatable::index_tmp                    !!! matrice storage
                double precision,dimension(:,:),allocatable::coords_pl           !!! coordinates storage
                double precision,dimension(3)::cent                              !!! Centroid of the LSQ plane
                double precision,dimension(4)::param                             !!! definition of the LSQ plane
                double precision,dimension(3,3)::evec                            !!!
                double precision,dimension(3)::eval                              !!!
                double precision,dimension(3)::axis                              !!!
                double precision::angle                                          !!!
                double precision,parameter::pi=4.0d0*atan(1.0d0)                 !!!
                double precision,dimension(3,3)::matop                           !!!
                logical::checker                                                 !!! Logical: .TRUE. if all atoms in the same plane.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                
                this%max_size=0
        do  j=1,this%n_groups!! Loop over the groups of equivalence.

        n=count(this%sym_groups(:,j) /= 0)  !! We loop only over the non-zero indices. 
        
        if (n>=4) then ! There can only be two planes for four atoms and more.
        allocate(coords_pl(3,n))
        
        do i=1,n        ! we copy the coordinates of the atoms in the group of equivalence j.
               coords_pl(:,i)=this%atoms_re(this%sym_groups(i,j))%coord(:)
        end do     
                call lsq(coords_pl,n,cent,evec,eval)   !! See LSQ subroutine in math.f95
      
                checker=.FALSE.

                do i=1,3
                        if ( eval(i).lt.1E-4) then    !!! here is a threshold.
                                checker=.TRUE.
                        end if
                end do

                if (checker.eqv..FALSE.) then !If there is two planes, we proceed and identify them

                        do i=1,3
                                do k=1,3
                                        if (ABS(evec(i,k)).lt.1E-4) then
                                        evec(i,k) = 0.0d0
                                        elseif (evec(i,k).lt.1.0001.AND.evec(i,k).gt.0.9999) then
                                        evec(i,k) = 1.0d0
                                        elseif (evec(i,k).gt.-1.0001.AND.evec(i,k).lt.-0.9999) then
                                                evec(i,k) = -1.0d0
                                        end if
                                end do
                        end do

                        !! We first to identify the good LSQ.
                        i=0
                        do while (checker.eqv..FALSE. &  !! Stops when the proper LSQ plane is identified.
                                & .AND.i.le.3)           !! Prevent infinite loop
                                i=i+1
                                matop=0.0d0
                                axis(:)=evec(:,i)
                                angle=pi/n
                                call get_rot(matop,axis,angle)                  !! There need to be a axis of rotation of order k
                                call get_equivalence(this,matop,checker)        !! k being (n of group)/2 as there are two planes.
                                param(1)=evec(1,i)  !!A
                                param(2)=evec(2,i)  !!B
                                param(3)=evec(3,i)  !!C
                                param(4)=cent(1)*evec(1,i)+cent(2)*evec(2,i)+cent(3)*evec(3,i) !!D
                        end do 
                        
                        allocate(index_tmp(this%size_groups,this%n_groups+1))  !! Allocation of the temporary index matrix with an additional group.
                        do i=1,this%n_groups
                                index_tmp(:,i)=this%sym_groups(:,i)  !! We copy the matrice over
                                if (i == j) then
                                        index_tmp(:,i)=0             !! we initialise the current group of interest with 0s
                                end if
                        end do
                        index_tmp(:,this%n_groups+1)=0               !! we initialise the new group with 0s.
                        
                        do i=1,n                                      !! we loop over the non zero indices // n = nb of atom in the group
                                if (param(1)*coords_pl(1,i)+param(2)*coords_pl(2,i)+param(3)*coords_pl(3,i)-param(4)&
                                &/(param(1)**2+param(2)**2 +param(3)**2) > 0 ) then    !! If the projection of the point to the plane if +1, then it belongs to plane I
                                        index_tmp(i,j)=this%sym_groups(i,j)
                                else
                                        index_tmp(i,1+this%n_groups)=this%sym_groups(i,j) !! If the projection is -1, it belongs to plane II
                                end if
                        end do

                        this%n_groups=this%n_groups+1                           !! We in crease the number of groups by 1.

                        tmp=0
                        !! We look for the larger number of equivalent of atom, and to which group it corresponds
                        do i=1,this%n_groups
                                 tmp=MAX(tmp,count(index_tmp(:,i) /= 0 ))
                                if ( tmp == MAX(tmp,count(index_tmp(:,i) /= 0 ))) then 
                                        this%max_size=i 
                                end if
                        end do
                        this%size_groups=tmp

                        !We reshape the matrice this%sym_groups
                        deallocate(this%sym_groups)
                        allocate(this%sym_groups(this%size_groups,this%n_groups))
                        this%sym_groups=0
                        do i=1,this%n_groups
                        ! And we copy over all the non 0 elements of the matrice.
                        this%sym_groups(:,i)=pack(index_tmp(:,i),index_tmp(:,i) /= 0 )
                        end do
                else !!! If all in the same plane
        
                        !! We look for the larger number of equivalent of atom, and to which group it corresponds
                        do i=1,this%n_groups
                               if ( n == MAX(n,count(this%sym_groups(:,i) /= 0 ))) then 
                                       this%max_size=j
                               end if
                       end do
        
                end if
                deallocate(coords_pl)

        end if
        end do 
              END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
              SUBROUTINE get_equivalence(this,mat,checker)                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!!  This routine check if every atom of the molecule has an equivalent one      !!!
!!!  through the symetry operiation provided in "mat". It output                 !!!
!!! "checker == .TRUE." if it is the case.                                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                class(mol_class),intent(inout)::this
                integer::i,j                                                     !!! Loop Integers
                double precision::normvec                                        !!! see function normvec in math.f95
                double precision::scalar                                         !!! see function scalar in math.f95
                double precision,dimension(3,3)::mat                             !!! transformation matrice
                double precision,dimension(3)::vec                               !!! transformed coordinates
                logical,dimension(:),allocatable::check                          !!! Logical mask
                logical::checker                                                 !!! Logical output.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                allocate(check(this%at))
                check=.FALSE.   ! Initialisation
                do i=1,this%at  ! loop over the atoms
                        vec=matmul(mat,this%atoms_re(i)%coord(:))  ! vec contains the transformed coordinates of atom i.
                             do j=1,this%at   ! loop over all the atoms.
                                if (check(j).eqv..FALSE.) then !! If no equivalence with atom i has yet be found then
                                        if (scalar(vec,this%atoms_re(j)%coord(:))&
                                        &/(normvec(this%atoms_re(j)%coord(:))**2) <=1.1 .AND. &
                                        &scalar(vec,this%atoms_re(j)%coord(:))&
                                        &/(normvec(this%atoms_re(j)%coord(:))**2) >=0.9 &
                                        &.AND. this%atoms(i)%label.eq.this%atoms(j)%label) then  !! If the atom are the same element and their scalar is comprised between 1+-tolerance then
                                                check(j)=.TRUE.                                  !! The two atoms are equivalents.
                                        else if (normvec(this%atoms_re(j)%coord(:)) == 0 .and. normvec(vec) == 0 ) then !! Molecule is centre of mass and should be on all elements of sym.
                                                check(j)=.TRUE.
                                         end if
                                end if
                        end do
                end do
        
                if ( count(check(:) .eqv. .TRUE.) == this%at ) then   !! If all atoms have equivalent atoms through "mat" then
                      checker=.TRUE.                                  !! outputs .TRUE.
                else                                                  !! Otherwise
                      checker=.FALSE.                                 !! outputs .FALSE.
                end if

                                deallocate(check) ! Deallocation to be clean(ish).
            end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        subroutine mol_print_out(this)                                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!! This is a debug priting routine. You wouldn't expect me to code it in one    !!!
!!! shot without checking intermediate results now, would you ?                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                class(mol_class),intent(inout)::this
                integer::i,j
                integer::k,l
                double precision::normvec

                write(*,'(a,a)') 'File Name: ', trim(adjustl(this%file_xyz))
                write(*,'(a,i0)') '#Atoms: ', this%at
                write(*,*)
                write(*,*) 'Input geometry:'
                write(*,*)
                write(*,'(a,a,a,a,a,a,a,a)') 'At',' | ', '     x      ','     y      ','     z      ',' | '&
                &,'    Mass    ','    Dist    '
                write(*,'(a)') '---+--------------------------------------+-------------'
                do i=1,this%at
                write(*,'(a,a,3f12.6,a,2f12.6)') this%atoms(i)%label,' | ', (this%atoms(i)%coord(j),j=1,3) ,' | ', &
                &this%atoms(i)%mass, this%atoms(i)%dist_com
                end do
                write(*,*)
                write(*,*) 'Centred geometry:'
                write(*,*)
                write(*,'(a,a,a,a,a,a,a,a)') 'At',' | ', '     x      ','     y      ','     z      ',' | '&
                &,'    Mass    ','    Dist    '
                write(*,'(a)') '---+--------------------------------------+-------------'
                do i=1,this%at
                write(*,'(a,a,3f12.6,a,2f12.6)') this%atoms(i)%label,' | ', &
                &(this%atoms_re(i)%coord(j),j=1,3) ,' | ', this%atoms(i)%mass, this%atoms_re(i)%dist_com
                end do
                write(*,*)
                write(*,'(a)') 'Inertia matrix;'
                write(*,'(3f12.3)') this%inertia_mat
                write(*,*)
                write(*,'(a)') 'Moment of Inertia'
                write(*,'(3f12.3)') this%inertia_moment
                write(*,*)
                write(*,*) 'Rotated geometry: (to align inertia axis)'
                write(*,*)
                write(*,'(a,a,a,a,a,a,a,a)') 'At',' | ', '     x      ','     y      ','     z      ',' | '&
                &,'    Mass    ','    Dist    '
                write(*,'(a)') '---+--------------------------------------+-------------'
                do i=1,this%at
                write(*,'(a,a,3f12.6,a,2f12.6)') this%atoms(i)%label,' | ', &
                &(this%atoms_re(i)%coord(j),j=1,3) ,' | ', this%atoms(i)%mass, this%atoms_re(i)%dist_com
                end do
                write(*,*)
        end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module mod_atom
