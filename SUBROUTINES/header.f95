        SUBROUTINE header()
        IMPLICIT NONE
        character(len=8)  :: date
        character(len=10) :: time
        character(len=5)  :: zone
        integer,dimension(8) :: values


        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)

        Write(*,*)
        Write(*,*)
        
write(*,*) '██████╗  █████╗  █',&
&'█████╗ ███████╗███╗   ██╗'
write(*,*) '██╔══██╗██╔══██╗██',&
&'╔════╝ ██╔════╝████╗  ██║'
write(*,*) '██║  ██║███████║██',&
&'║  ███╗█████╗  ██╔██╗ ██║'
write(*,*) '██║  ██║██╔══██║██',&
&'║   ██║██╔══╝  ██║╚██╗██║'
write(*,*) '██████╔╝██║  ██║╚█',&
&'█████╔╝███████╗██║ ╚████║'
write(*,*) '╚═════╝ ╚═╝  ╚═╝ ╚',&
&'═════╝ ╚══════╝╚═╝  ╚═══╝'
        
        write(*,*) '                                        v1.0          '
        Write(*,*)
     !   Write(*,*) '+=============================================================+'
     !   write(*,*) '|  J. Eng*                                                    |'             
     !   write(*,*) '|  * Chemistry, School of Natural and Environmental Sciences, |'
     !   write(*,*) '|    Newcastle University, NE1 7RU, Newcastle Upon Tyne, UK   |'
     !   write(*,*) '|                                                             |' 
     !   write(*,*) '| VCMaker website: VCMaker.glitch.me                          |' 
     !   write(*,*) '|                                                             |' 
     !   write(*,*) '| Cite me as (bibtex entry at the end):                       |'
     !   write(*,*) '| T.J. Penfold, J. Eng, Phys. Chem. Chem. Phys.               |'
     !   write(*,*) '| 25, 7195-7204 (2023)                                        |'
     !   write(*,*) '| VCMaker, 2022, github.com/JulienEng/VCMaker                 |'
     !   write(*,*) '|                                                             |' 
     !   write(*,*) '+=============================================================+'
     !   write(*,*)
        write(*,*) 'Started on the ', date(7:8),'/',date(5:6),'/',&
        &date(1:4),' at ', time(1:2),':',time(3:4),':',time(5:6),'.'
        write(*,*)


        END SUBROUTINE

        SUBROUTINE endoftimes()
        IMPLICIT NONE
        character(len=8)  :: date
        character(len=10) :: time
        character(len=5)  :: zone
        integer,dimension(8) :: values

        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)

   !     write(*,*) '*** BibTeX Entry ***'
   !     write(*,*) '@article{GAP_2023, '
   !     write(*,*) 'title={Mind the GAP: quantifying the breakdown of the linear vibronic coupling Hamiltonian}, '
   !     write(*,*) 'volume={25}, '
   !     write(*,*) 'ISSN={1463-9084}, '
   !     write(*,*) 'url={http://dx.doi.org/10.1039/D2CP05576G}, '
   !     write(*,*) 'DOI={10.1039/d2cp05576g}, '
   !     write(*,*) 'number={10}, '
   !     write(*,*) 'journal={Physical Chemistry Chemical Physics}, '
   !     write(*,*) 'publisher={Royal Society of Chemistry (RSC)}, '
   !     write(*,*) 'author={Penfold, Thomas J and Eng, Julien}, '
   !     write(*,*) 'year={2023}, '
   !     write(*,*) 'pages={7195–7204} }'
   !     write(*,*)
   !     write(*,*) '@misc{VCMaker,'
   !     write(*,*) 'author = {Julien Eng},'
   !     write(*,*) 'title = {VCMaker (v1)},'
   !     write(*,*) 'year  = {2022},'
   !     write(*,*) 'url   = {github.com/JulienEng/VCMaker} }'

        write(*,*)
        write(*,*) 'Stopped on the ', date(7:8),'/',date(5:6),&
        &'/',date(1:4),' at ', time(1:2),':',time(3:4),':',time(5:6),'.'
        write(*,*) ' Thank you for using DAGen!'
        END SUBROUTINE
