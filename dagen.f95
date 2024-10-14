Program dagen
use mod_atom
implicit none
type(mol_class)::donor
type(mol_class)::acceptor
type(mol_class)::DA
type(mol_class),dimension(:),allocatable::PICT_DA

character(len=666)::arg3,arg4

!Integers for loops
integer::i,j,k,l

!Integers for atoms
integer::link_at_D
integer::link_at_A
integer::H_at_D
integer::H_at_A

integer,dimension(2)::neigh_C

!Integers for counters
integer::counter_at
integer::counter_clash

!! Clash stuff
integer,dimension(:,:),allocatable::H_clash
integer,dimension(:,:),allocatable::C_clash
integer,dimension(:,:),allocatable::H_clash_tmp
integer,dimension(:,:),allocatable::C_clash_tmp

!! Math Stuff
double precision::normvec
double precision::scalar

!! Al-Kashi Theroem
double precision::alpha,beta,gamma              ! Angles
double precision,dimension(3)::a,b,c            ! Distances
double precision::angle,curr_angle
double precision,dimension(3,3)::rot_mat

!! Alignment stuff
double precision,dimension(3)::axe_y,axe_z
double precision,dimension(3)::vec,axe,veca,vecb
double precision,dimension(3)::shift

!! OUT
character(len=123)::file_out

!!! Get Arguments
call getarg(1,donor%file_xyz)
call getarg(2,acceptor%file_xyz)
call getarg(3,arg3) ; read(arg3,*) link_at_D
call getarg(4,arg4) ; read(arg4,*) link_at_A


!!! Read Coordinates
call mol_read_xyz(donor)
call mol_read_xyz(acceptor)


axe_z=0
axe_z(3)=1

!!!!! REORIENTATION TO ALIGN THE PLANE !!!!!

!!!! Loog at the link atom's neighbours and aligne the plane DONOR...
    counter_at=0

do i=1,donor%at
    if (donor%atoms(i)%label.ne.'H'.AND.i.ne.link_at_D) then
      if (normvec(donor%atoms(i)%coord(:) - donor%atoms(link_at_D)%coord(:)).lt.1.9) then
      counter_at=counter_at+1
        neigh_C(counter_at)=i
      end if
    end if
end do
!!! WORKS
vec(:)=donor%atoms(link_at_D)%coord(:)
do i=1,donor%at
donor%atoms(i)%coord(:)=donor%atoms(i)%coord(:)-vec(:)
end do
veca=donor%atoms(neigh_C(1))%coord-donor%atoms(link_at_D)%coord
vecb=donor%atoms(neigh_C(2))%coord-donor%atoms(link_at_D)%coord

call vecto(veca,vecb,vec)

angle=acos(scalar(vec,axe_z)/(normvec(vec)*normvec(axe_z)))
call vecto(vec,axe_z,axe)
call get_rot(rot_mat,axe,angle)
do i=1,donor%at
donor%atoms(i)%coord(:)=MATMUL(rot_mat,donor%atoms(i)%coord(:))
end do

!!!! Loog at the link atom's neighbours and aligne the plane ACCEPTOR...
counter_at=0
do i=1,acceptor%at
    if (acceptor%atoms(i)%label.ne.'H'.AND.i.ne.link_at_A) then
      if (normvec(acceptor%atoms(i)%coord(:) - acceptor%atoms(link_at_A)%coord(:)).lt.1.8) then
      counter_at=counter_at+1
        neigh_C(counter_at)=i
      end if
    end if
end do
!!! WORKS
vec(:)=acceptor%atoms(link_at_A)%coord(:)
do i=1,acceptor%at
acceptor%atoms(i)%coord(:)=acceptor%atoms(i)%coord(:)-vec(:)
end do

veca=acceptor%atoms(neigh_C(1))%coord-acceptor%atoms(link_at_D)%coord
vecb=acceptor%atoms(neigh_C(2))%coord-acceptor%atoms(link_at_D)%coord
call vecto(veca,vecb,vec)

angle=ACOS(scalar(vec,axe_z)/(normvec(vec)*normvec(axe_z)))
call vecto(vec,axe_z,axe)
call get_rot(rot_mat,axe,angle)
do i=1,acceptor%at
acceptor%atoms(i)%coord(:)=MATMUL(rot_mat,acceptor%atoms(i)%coord(:))
end do
!!!!! END - REORIENTATION TO ALIGN THE PLANE !!!!!



!!!! Identification of H atoms bonded to link atom !!!!
do i=1,donor%at
  if (donor%atoms(i)%label.eq.'H') then
    if (normvec(donor%atoms(i)%coord(:) - donor%atoms(link_at_D)%coord(:)).lt.1.2) then
      H_at_D=i
    end if
  end if
end do

do i=1,acceptor%at
  if (acceptor%atoms(i)%label.eq.'H') then
    if (normvec(acceptor%atoms(i)%coord(:) - acceptor%atoms(link_at_A)%coord(:)).lt.1.2) then
      H_at_A=i
    end if
  end if
end do
!!!! END - Identification of H atoms bonded to link atom !!!!

Write(*,'(a,a,i0,a,a,i0)') 'Donor -link- bond: ',donor%atoms(link_at_D)%label,link_at_D,'-' ,&
&donor%atoms(H_at_D)%label,H_at_D
Write(*,'(a,a,i0,a,a,i0)') 'Acceptor -link- bond: ',acceptor%atoms(link_at_A)%label,link_at_A,'-' ,&
&acceptor%atoms(H_at_A)%label,H_at_A

!!!! Centre FRAGMENTS' link atom on 0 !!!!

vec(:)=acceptor%atoms(link_at_A)%coord(:)
do i=1,acceptor%at
acceptor%atoms(i)%coord(:)=acceptor%atoms(i)%coord(:)-vec(:)
end do
vec(:)=donor%atoms(link_at_D)%coord(:)
do i=1,donor%at
donor%atoms(i)%coord(:)=donor%atoms(i)%coord(:)-vec(:)!
end do
!!!! END - Centre FRAGMENTS' link atom on 0 !!!!
axe_z=0
axe_z(3)=1
!!!! Checking the orientation !!!!
axe_y=0.0d0 ; axe_y(2)=1
vec(:)=(donor%atoms(H_at_D)%coord(:) - donor%atoms(link_at_D)%coord(:)) ; vec(3)=0.0d0

if ((scalar((vec/normvec(vec)),axe_y(:)).lt.0).OR.&
&(ABS(scalar((vec/normvec(vec)),axe_y(:))-1).gt.1E-4)) then
angle=acos(scalar(vec,axe_y)/(normvec(vec)*normvec(axe_y)))
call get_rot(rot_mat,axe_z,angle)
do i=1,donor%at
  donor%atoms(i)%coord(:)=MATMUL(rot_mat,donor%atoms(i)%coord(:))
end do
end if

axe_y=0.0d0 ; axe_y(2)=-1
vec(:)=(acceptor%atoms(H_at_A)%coord(:) - acceptor%atoms(link_at_A)%coord(:)) ; vec(3)=0.0d0

if ((scalar((vec/normvec(vec)),axe_y(:)).lt.0).OR.&
& (ABS(scalar((vec/normvec(vec)),axe_y(:))-1).gt.1E-4)) then
 
angle=-acos(scalar(vec,axe_y)/(normvec(vec)*normvec(axe_y)))
call get_rot(rot_mat,axe_z,angle)
do i=1,acceptor%at
  acceptor%atoms(i)%coord(:)=MATMUL(rot_mat,acceptor%atoms(i)%coord(:))
end do
end if
!!!! END - Checking the orientation !!!!

!! Adding Shift
shift=0.0d0
shift(2)=-1.2 ! Angstroem
do i=1,donor%at
donor%atoms(i)%coord(:)=donor%atoms(i)%coord(:)+shift(:)
end do

!! END - Adding Shift



!!!! PRINTING THE D-A molecule !!!!
 write(file_out,'(a,i0,a)') 'out_',0,'.xyz'
  open(69,file=trim(adjustl(file_out)))
DA%at=donor%at+acceptor%at - 2
write(69,'(i0)') DA%at
write(69,*)
counter_at=0
allocate(DA%atoms(DA%at))
do i=1,donor%at
   if (i.ne.H_at_D) then
    counter_at=counter_at+1
    DA%atoms(counter_at)%coord(:)=donor%atoms(i)%coord(:)
    DA%atoms(counter_at)%label=donor%atoms(i)%label
    write(69,'(A2,3F15.8)') DA%atoms(counter_at)%label, DA%atoms(counter_at)%coord(:)
   end if
end do
do i=1,acceptor%at
   if (i.ne.H_at_A) then
    counter_at=counter_at+1
    DA%atoms(counter_at)%coord(:)=acceptor%atoms(i)%coord(:)
    DA%atoms(counter_at)%label=acceptor%atoms(i)%label
    write(69,'(A2,3F15.8)') DA%atoms(counter_at)%label, DA%atoms(counter_at)%coord(:)
   end if
end do
close(69)
!!!! END - PRINTING THE D-A molecule !!!!

!!!! Checking for Clashes !!!!
counter_clash=0
do i=1,donor%at
  if (donor%atoms(i)%label.eq.'H') then
    do j=1,acceptor%at
      if (acceptor%atoms(j)%label.eq.'H') then
        vec(:)=donor%atoms(i)%coord(:)-acceptor%atoms(j)%coord(:)
        if (normvec(vec).lt.1.1.AND.i.ne.H_at_D.AND.j.ne.H_at_A) then

          counter_clash=counter_clash+1

              if (counter_clash.gt.1) then
               if (allocated(H_clash_tmp)) then ; deallocate(H_clash_tmp) ; end if ; allocate(H_clash_tmp(counter_clash-1,2))
               H_clash_tmp=H_clash
               deallocate(H_clash)
              end if

            allocate(H_clash(counter_clash,2))

            if (counter_clash.gt.1) then
              do k=1,counter_clash-1
                H_clash(k,:)=H_clash_tmp(k,:)
              end do
            end if
              H_clash(counter_clash,1)=i  !! DONOR H
              H_clash(counter_clash,2)=j  !! ACCEPTOR H
        end if
      end if
    end do
  end if
end do
!!!! END - Checking for Clashes !!!!

write(*,*) 'Number of pairs of close lying H atoms: ', counter_clash

!!!! Check for Atom bonded to the Clashing H !!!!
allocate(C_clash(counter_clash,2))
do i=1,counter_clash
  do j=1,donor%at
    if (donor%atoms(j)%label.ne.'H') then
      vec(:)=donor%atoms(j)%coord(:)-donor%atoms(H_clash(i,1))%coord(:)
      if (normvec(vec).lt.1.2 ) then
      C_clash(i,1)=j
      end if
    end if
  end do
    do j=1,acceptor%at
    if (acceptor%atoms(j)%label.ne.'H') then
      vec(:)=acceptor%atoms(j)%coord(:)-acceptor%atoms(H_clash(i,2))%coord(:)
      if (normvec(vec).lt.1.2 ) then
      C_clash(i,2)=j
      end if
    end if
  end do
end do
!!!! END - Check for Atom bonded to the Clashing H !!!!

allocate(PICT_DA(counter_clash))

!!!! Al-Kashi To tilt the acceptor. !!!!
do i=1,counter_clash
 a(:)=acceptor%atoms(C_clash(i,2))%coord(:)-acceptor%atoms(link_at_a)%coord(:)
 c(:)=donor%atoms(C_clash(i,1))%coord(:)-acceptor%atoms(link_at_a)%coord(:)
curr_angle=acos(scalar(a,c)/(normvec(a)*normvec(c)))

beta=acos((1.9**2-normvec(a)**2-normvec(c)**2)/(-2*normvec(a)*normvec(c)))

angle=curr_angle-beta
call vecto(a,c,vec)
call get_rot(rot_mat,vec,angle)
counter_at=0
PICT_DA(i)%at=acceptor%at + donor%at - 4
allocate(PICT_DA(i)%atoms(PICT_DA(i)%at))
    do j=1,donor%at
      if (j.ne.H_at_D.AND.j.ne.H_clash(i,1)) then
        counter_at=counter_at+1
        PICT_DA(i)%atoms(counter_at)%coord(:)=donor%atoms(j)%coord(:)
        PICT_DA(i)%atoms(counter_at)%label=donor%atoms(j)%label

      end if
    end do
      do j=1,acceptor%at
      if (j.ne.H_at_A.AND.j.ne.H_clash(i,2)) then
        counter_at=counter_at+1
        PICT_DA(i)%atoms(counter_at)%coord(:)=MATMUL(rot_mat,acceptor%atoms(j)%coord(:))
        PICT_DA(i)%atoms(counter_at)%label=acceptor%atoms(j)%label

      end if
    end do  

!!!! Printing the isomer. !!!!
  write(file_out,'(a,i0,a)') 'out_',i,'.xyz'
  open(69,file=trim(adjustl(file_out)))
   write(69,'(i0)') acceptor%at + donor%at - 4
  write(69,*)
  do j=1,acceptor%at + donor%at - 4
    write(69,'(A2,3F15.8)') PICT_DA(i)%atoms(j)%label , PICT_DA(i)%atoms(j)%coord(:)
    end do


end do
!!!! END - Al-Kashi To tilt the acceptor. !!!!


end program