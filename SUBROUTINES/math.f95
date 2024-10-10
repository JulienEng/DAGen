!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!! Maths         
       function normvec(vec)
        implicit none
        double precision::normvec
        double precision,dimension(3)::vec

normvec=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine vecto(veca,vecb,vectorial)
        implicit none
        double precision,dimension(3),intent(IN)::veca,vecb
        double precision,dimension(3),intent(OUT)::vectorial
vectorial(1)=(veca(2)*vecb(3)-veca(3)*vecb(2))
vectorial(2)=(veca(3)*vecb(1)-veca(1)*vecb(3))
vectorial(3)=(veca(1)*vecb(2)-veca(2)*vecb(1))
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function scalar(veca,vecb)
        implicit none
        double precision::scalar
        double precision,dimension(3)::veca,vecb

        scalar=veca(1)*vecb(1)+veca(2)*vecb(2)+veca(3)*vecb(3)
        return
 end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE diag(init_mat,n,eigenvec,eigenval)
        implicit none
        integer::n
        double precision,dimension(n,n)::eigenvec,init_mat
        double precision,dimension(n)::eigenval
        ! Diagonalization routine
          integer::mw_lwork
          integer::mw_info
          double precision,dimension(:),allocatable::mw_work
  
        eigenvec=init_mat
  
        mw_lwork=-1
        allocate(mw_work(n))
  
        call DSYEV('V','U',n,eigenvec,n,eigenval,&
        &mw_work,mw_lwork,mw_info)
  
        mw_lwork=INT(mw_work(1))
        deallocate(mw_work)
        allocate(mw_work(mw_lwork))
  
        call DSYEV('V','U',n,eigenvec,n,eigenval,&
        &mw_work,mw_lwork,mw_info)
        deallocate(mw_work)
  
        END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        SUBROUTINE lsq(coord,n,cent,evec,eval)
          IMPLICIT NONE
          INTEGER::i,j,k,l
          INTEGER::n
          DOUBLE PRECISION,DIMENSION(3,n)::coord
          DOUBLE PRECISION,DIMENSION(3)::cent
          DOUBLE PRECISION,DIMENSION(4)::param
          DOUBLE PRECISION,DIMENSION(3,3)::mat,evec
          double precision,dimension(3)::dc,eval
  
          !!! Set to 0
          cent=0.0d0
          param=0.0d0
          mat=0.0d0
          dc=0.0d0
  
          !! Determination of the centroid point
          do j=1,3
              do i=1,n
                   cent(j)=cent(j)+coord(j,i)
              end do
          cent(j)=cent(j)/DBLE(n)
          end do
  
         !!! Tensor
         do i=1,n !! Loop over the atoms
           do j=1,3
           dc(j)=coord(j,i)-cent(j)
           end do
           do j=1,3
             do k=1,3
           mat(j,k)=mat(j,k)+dc(j)*dc(k)
             end do
           end do
         end do !! loop over the atoms
  
         CALL diag(mat,3,evec,eval)

  ! RMS points to plane is: RMS=sqrt(eigenvalue(i)/n)
  ! eigenvectors(XX,i) : i represent an intermediate and worst LSQ 

        END SUBROUTINE       
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
