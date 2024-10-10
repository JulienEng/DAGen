subroutine get_inv(mat)
    implicit none
    double precision,dimension(3,3)::mat
    integer::i

    mat=0.0d0

    do i=1,3
        mat(i,i)=-1.0d0
    end do
end subroutine

subroutine get_sigma(mat,normal)
    implicit none
    double precision,dimension(3,3)::mat
    integer::i,j,k
    double precision,dimension(3)::normal
    double precision::normvec
    double precision,dimension(3)::norm_axis

    norm_axis=normal/normvec(normal)

    mat=0.0d0
    do i=1,3 ; mat(i,i) = 1.0d0 ; end do
    
    do i=1,3
        do j=1,3
            mat(i,j)=mat(i,j) - 2*norm_axis(i)*norm_axis(j)
        end do
    end do


    !mat=0.0d0
    !do i=1,3
    !    mat(i,i)=norm_axis(1)**2+norm_axis(2)**2+norm_axis(3)**2
    !    mat(i,i)=mat(i,i)-2*norm_axis(i)**2
    !    if (i.lt.3) then
    !    do j=i+1,3
    !        do k=1,3
    !            if (k.ne.i.and.k.ne.j) then
    !        mat(i,j)=(-2*norm_axis(1)*norm_axis(2)*norm_axis(3))/norm_axis(k)
    !            end if
    !            mat(i,j)=mat(j,i)
    !        end do
    !    end do
    !    end if
    !end do

end subroutine

subroutine get_rot(mat,axis,angle)
    implicit none
    double precision::angle
    double precision,dimension(3)::axis
    double precision,dimension(3,3)::mat
    double precision,dimension(3,3)::matK
    double precision,dimension(3,3)::matI
    integer::i,j,k
    double precision::normvec
    double precision,dimension(3)::norm_axis
    ! Compute the cosine, sine and (1 - cosine) of the angle
    ! and assign components of the normalized axis
    double precision :: c, s, v

    norm_axis=axis/normvec(axis)

    c=cos(angle)
    s=sin(angle)
    v=1.0d0 - c

    ! Explicitly define each element of the rotation matrix
    mat(1,1)=c+v*norm_axis(1)**2
mat(1,1)=c+norm_axis(1)**2*v
    mat(1,2)= norm_axis(1) * norm_axis(2) * v - norm_axis(3) * s
    mat(1,3)= norm_axis(1) * norm_axis(3) * v + norm_axis(2) * s
    mat(2,1)= norm_axis(2) * norm_axis(1) * v + norm_axis(3) * s
    mat(2,2)= c + norm_axis(2)**2 * v
    mat(2,3)= norm_axis(2) * norm_axis(3) * v - norm_axis(1) * s

    mat(3,1)= norm_axis(3) * norm_axis(1) * v - norm_axis(2) * s
    mat(3,2)= norm_axis(3) * norm_axis(2) * v + norm_axis(1) * s
    mat(3,3)= c + norm_axis(3)**2 * v

end subroutine


subroutine get_sn(mat,axis,angle)
    implicit none
    double precision::angle
    double precision,dimension(3)::axis
    double precision,dimension(3,3)::mat,mata,matb

    call get_rot(mata,axis,angle)
    call get_sigma(matb,axis)

    mat=matmul(mata,matb)

end subroutine

subroutine get_id(mat)
    implicit none
    double precision,dimension(3,3)::mat
    integer::i
do i=1,3
mat(i,i)=1.0d0
end do

end subroutine