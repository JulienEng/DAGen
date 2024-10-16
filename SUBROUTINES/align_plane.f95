subroutine align_plane(mol,linkatom)
use mod_atom
implicit none

integer::i
type(mol_class)::mol
integer::linkatom
double precision::normvec
double precision::scalar
double precision,dimension(3)::veca
double precision,dimension(3)::vecb
double precision,dimension(3)::vec
double precision,dimension(3)::axe
double precision::angle

double precision,dimension(3),parameter::axe_y=[ 0 , 1 , 0 ]
double precision,dimension(3),parameter::axe_z=[ 0 , 0 , 1 ]

double precision,dimension(3,3)::rot_mat

integer,dimension(2)::neigh_C

integer::counter_at

!!!! Loog at the link atom's neighbours and aligne the plane mol...
    counter_at=0

do i=1,mol%at
    if (mol%atoms(i)%label.ne.'H'.AND.i.ne.linkatom) then
      if (normvec(mol%atoms(i)%coord(:) - mol%atoms(linkatom)%coord(:)).lt.1.9) then
      counter_at=counter_at+1
        neigh_C(counter_at)=i
      end if
    end if
end do

! centring on linkatom
vec(:)=mol%atoms(linkatom)%coord(:)
do i=1,mol%at
mol%atoms(i)%coord(:)=mol%atoms(i)%coord(:)-vec(:)
end do
veca=mol%atoms(neigh_C(1))%coord-mol%atoms(linkatom)%coord
vecb=mol%atoms(neigh_C(2))%coord-mol%atoms(linkatom)%coord

call vecto(veca,vecb,vec)

angle=acos(scalar(vec,axe_z)/(normvec(vec)*normvec(axe_z)))
call vecto(vec,axe_z,axe)
call get_rot(rot_mat,axe,angle)

do i=1,mol%at
mol%atoms(i)%coord(:)=MATMUL(rot_mat,mol%atoms(i)%coord(:))
end do

veca=mol%atoms(neigh_C(1))%coord-mol%atoms(linkatom)%coord
vecb=mol%atoms(neigh_C(2))%coord-mol%atoms(linkatom)%coord

call vecto(veca,vecb,vec)

end subroutine



subroutine center_on_at(mol,linkatom,shift,axe,Hat)
use mod_atom
implicit none

integer::i
type(mol_class)::mol
integer::linkatom
integer::Hat
double precision,dimension(3)::shift
double precision,dimension(3)::axe
double precision,dimension(3)::vec
double precision::angle
double precision,dimension(3,3)::rot_mat
double precision::normvec,scalar
double precision,dimension(3),parameter::axe_y=[ 0 , 1 , 0 ]
double precision,dimension(3),parameter::axe_z=[ 0 , 0 , 1 ]

!!!! Identification of H atoms bonded to link atom !!!!
do i=1,mol%at
  if (mol%atoms(i)%label.eq.'H') then
    if (normvec(mol%atoms(i)%coord(:) - mol%atoms(linkatom)%coord(:)).lt.1.2) then
      Hat=i
    end if
  end if
end do
!!!! END - Identification of H atoms bonded to link atom !!!!

!!!! Checking the orientation !!!!
vec(:)=(mol%atoms(Hat)%coord(:) - mol%atoms(linkatom)%coord(:)) ; vec(3)=0.0d0

if ((scalar((vec/normvec(vec)),axe(:)).lt.0).OR.&
&(ABS(scalar((vec/normvec(vec)),axe(:))-1).gt.1E-4)) then
angle=-acos(scalar(vec,axe)/(normvec(vec)*normvec(axe)))
call get_rot(rot_mat,axe_z,angle)
do i=1,mol%at
  mol%atoms(i)%coord(:)=MATMUL(rot_mat,mol%atoms(i)%coord(:))
  mol%atoms(i)%coord(:)=mol%atoms(i)%coord(:)+shift(:)
end do
end if


end subroutine